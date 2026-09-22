"""Preflight validation for EPFinder preprocessing inputs."""

import os
import shutil


def normalize_ensembl_id(identifier):
    return identifier.strip().split('.', 1)[0]


def normalize_chr(c):
    """Normalize chromosome label to no 'chr' prefix: 'chr1'->'1', 'X'->'X'."""
    c = c.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def detect_chr_prefix(path, column=0):
    """True when the first data row's chromosome column carries a 'chr' prefix."""
    for _, fields in _rows(path):
        return fields[column].strip().lower().startswith('chr')
    return False


def format_chr(chromosome, use_prefix):
    """Render a normalized chromosome label in the requested convention."""
    chromosome = normalize_chr(chromosome)
    return 'chr' + chromosome if use_prefix else chromosome


def _rows(path):
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.strip() and not raw.startswith('#'):
                yield number, raw.rstrip('\n').split('\t')


def _file(path, label):
    if not os.path.isfile(path):
        raise ValueError(f"{label} is not a readable file: {path}")


def _int(value, label, minimum=0):
    try:
        value = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be an integer: {value!r}") from None
    if value < minimum:
        raise ValueError(f"{label} must be >= {minimum}: {value}")
    return value


def _float(value, label):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be numeric: {value!r}") from None


def _expression(path, prefix):
    ids = set()
    for number, fields in _rows(path):
        if not fields[0].startswith(prefix):
            continue
        if len(fields) < 2:
            raise ValueError(f"{path}:{number}: expression rows need ID and value")
        ident = normalize_ensembl_id(fields[0])
        _float(fields[-1], f"{path}:{number}: expression")
        ids.add(ident)
    if not ids:
        raise ValueError(f"{path}: no {prefix} expression rows found")
    return ids


def validate_inputs(config):
    """Validate configured files and cross-file identifier compatibility."""
    required = ('input_gwas', 'hic_folder', 'hic_prefix', 'tss_file',
                'tx_expression', 'gene_list', 'gene_expression',
                'feature_list', 'output_dir', 'enhancer_window',
                'promoter_window', 'hic_bin_size')
    missing = [key for key in required if config.get(key) in (None, '')]
    if missing:
        raise ValueError('Missing required config keys: ' + ', '.join(missing))
    for key in ('input_gwas', 'tss_file', 'tx_expression', 'gene_list',
                'gene_expression', 'feature_list'):
        _file(config[key], key)
    for key in ('enhancer_window', 'promoter_window', 'hic_bin_size'):
        _int(config[key], key, 1)
    _int(config.get('step1_nproc', 4), 'step1_nproc', 1)
    bedtools = config.get('bedtools_path', 'bedtools')
    if not shutil.which(bedtools):
        raise ValueError(f"bedtools executable not found: {bedtools}")

    chromosomes = set()
    snps = 0
    gwas_columns = None
    for number, fields in _rows(config['input_gwas']):
        if len(fields) < 2:
            raise ValueError(f"{config['input_gwas']}:{number}: need chromosome and position")
        if gwas_columns is None:
            gwas_columns = len(fields)
        elif len(fields) != gwas_columns:
            raise ValueError(
                f"{config['input_gwas']}:{number}: every row must have the same number of "
                f"columns (found {len(fields)}, expected {gwas_columns}). Downstream steps "
                f"locate the Hi-C and TSS columns by offset from this width."
            )
        chromosome = normalize_chr(fields[0])
        if not chromosome:
            raise ValueError(f"{config['input_gwas']}:{number}: empty chromosome")
        _int(fields[1], f"{config['input_gwas']}:{number}: position", 1)
        chromosomes.add(chromosome)
        snps += 1
    if not snps:
        raise ValueError('input_gwas contains no data rows')
    config['gwas_columns'] = gwas_columns

    if not os.path.isdir(config['hic_folder']):
        raise ValueError(f"hic_folder is not a directory: {config['hic_folder']}")
    for chromosome in chromosomes:
        path = os.path.join(config['hic_folder'], config['hic_prefix'] + 'chr' + chromosome)
        _file(path, f'Hi-C chromosome {chromosome}')
        for number, fields in _rows(path):
            if len(fields) != 3:
                raise ValueError(f"{path}:{number}: Hi-C rows need 3 columns")
            _int(fields[0], f"{path}:{number}: bin1")
            _int(fields[1], f"{path}:{number}: bin2")
            if fields[2].lower() != 'nan':
                _float(fields[2], f"{path}:{number}: contact")
            break

    config['tss_chr_prefix'] = detect_chr_prefix(config['tss_file'])
    tss_tx, tss_genes = set(), set()
    previous = None
    for number, fields in _rows(config['tss_file']):
        if len(fields) < 5:
            raise ValueError(f"{config['tss_file']}:{number}: TSS rows need >=5 columns")
        start = _int(fields[1], f"{config['tss_file']}:{number}: start")
        end = _int(fields[2], f"{config['tss_file']}:{number}: end")
        if end < start:
            raise ValueError(f"{config['tss_file']}:{number}: end before start")
        tx = normalize_ensembl_id(fields[-2])
        if not tx.startswith('ENST') or not fields[-1]:
            raise ValueError(f"{config['tss_file']}:{number}: invalid transcript/gene columns")
        key = (fields[0], start)
        if previous is not None and key < previous:
            raise ValueError(f"{config['tss_file']}:{number}: TSS file is not sorted")
        previous = key
        tss_tx.add(tx); tss_genes.add(fields[-1])
    tx_ids = _expression(config['tx_expression'], 'ENST')
    gene_ids = _expression(config['gene_expression'], 'ENSG')
    list_ids, list_genes = set(), set()
    for number, fields in _rows(config['gene_list']):
        if len(fields) < 6 or not fields[4] or not normalize_ensembl_id(fields[5]).startswith('ENSG'):
            raise ValueError(f"{config['gene_list']}:{number}: expected symbol in col 5 and ENSG in col 6")
        list_genes.add(fields[4]); list_ids.add(normalize_ensembl_id(fields[5]))
    if not tss_tx & tx_ids or not list_ids & gene_ids:
        raise ValueError('Annotation and expression files have no matching stable Ensembl IDs')
    if not tss_genes & list_genes:
        raise ValueError('TSS and gene-list files have no matching gene symbols')

    feature_count = 0
    names = set()
    feature_prefix = None
    for number, fields in _rows(config['feature_list']):
        if len(fields) != 2 or not fields[0] or fields[0] in names:
            raise ValueError(f"{config['feature_list']}:{number}: need unique name and BED path")
        _file(fields[1], f"{config['feature_list']}:{number}: feature BED")
        prefix = detect_chr_prefix(fields[1])
        if feature_prefix is None:
            feature_prefix = prefix
        elif prefix != feature_prefix:
            raise ValueError(
                f"{fields[1]}: feature signal files disagree on the 'chr' prefix. All files "
                f"in feature_list must use the same chromosome naming convention, because "
                f"they are scored against one shared set of enhancer/promoter regions."
            )
        for bed_number, bed in _rows(fields[1]):
            if len(bed) < 4:
                raise ValueError(f"{fields[1]}:{bed_number}: BED rows need 4 columns")
            start = _int(bed[1], f"{fields[1]}:{bed_number}: start")
            end = _int(bed[2], f"{fields[1]}:{bed_number}: end")
            if end < start: raise ValueError(f"{fields[1]}:{bed_number}: end before start")
            _float(bed[3], f"{fields[1]}:{bed_number}: signal")
            break
        names.add(fields[0]); feature_count += 1
    if not feature_count:
        raise ValueError('feature_list contains no features')
    config['feature_chr_prefix'] = feature_prefix

    def style(prefix):
        return "chrN" if prefix else "N"

    print(f"Preflight passed: {snps} SNPs, {len(tss_tx & tx_ids)} transcripts, "
          f"{len(list_ids & gene_ids)} genes, {feature_count} features")
    print(f"  SNP file: {gwas_columns} columns, chromosomes as "
          f"{style(detect_chr_prefix(config['input_gwas']))}")
    print(f"  Chromosome naming: TSS={style(config['tss_chr_prefix'])}, "
          f"feature signals={style(config['feature_chr_prefix'])}")
