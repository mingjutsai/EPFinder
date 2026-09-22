"""Synthetic EPFinder inputs.

Small enough to run the whole 9-step workflow in under a second, and chosen so
every feature value can be worked out by hand:

  chr1  SNP 10000 -> Hi-C bin 9000, contacts bin 21000, which holds TSS 21500
  chr2  SNP 50000 -> Hi-C bin 48000, contacts bin 60000, which holds TSS 60500

With 250 bp windows a signal interval that spans a whole window contributes
500 bp * its value, so MarkX at the chr1 enhancer is 500 * 2.0 = 1000.

The chromosome style of each input is independent, which is the point: the
workflow has to detect and convert between them rather than assume one.
"""

import os


def _write(path, text):
    with open(path, 'w') as handle:
        handle.write(text)


def build(root, snp_style='bare', tss_style='bare', signal_style='chr',
          snp_columns=3):
    """Write one input set under `root` and return its config dict.

    `*_style` is 'bare' ("1") or 'chr' ("chr1"). The default mirrors the
    production eBMD layout: bare SNP and TSS files, chr-prefixed signal files.
    """
    root = str(root)
    hic = os.path.join(root, 'hic')
    os.makedirs(hic, exist_ok=True)

    def snp_c(n):
        return 'chr' + n if snp_style == 'chr' else n

    def tss_c(n):
        return 'chr' + n if tss_style == 'chr' else n

    def sig_c(n):
        return 'chr' + n if signal_style == 'chr' else n

    if snp_columns == 3:
        _write(os.path.join(root, 'snps.tsv'),
               '#Chr\tPos\tRSID\n'
               '%s\t10000\trs111\n' % snp_c('1') +
               '%s\t50000\trs222\n' % snp_c('2'))
    else:
        _write(os.path.join(root, 'snps.tsv'),
               '#Chr\tPos\tRSID\tPvalue\n'
               '%s\t10000\trs111\t1e-9\n' % snp_c('1') +
               '%s\t50000\trs222\t2e-10\n' % snp_c('2'))

    # Row 2 exercises the 'nan' contact path; row 3 a bin with no SNP in it.
    _write(os.path.join(hic, 'TEST.hic.KR.chr1'),
           '9000\t21000\t5.5\n9000\t30000\tnan\n12000\t21000\t1.1\n')
    _write(os.path.join(hic, 'TEST.hic.KR.chr2'), '48000\t60000\t7.25\n')

    _write(os.path.join(root, 'tss.tsv'),
           '%s\t21500\t21500\tENST00000000001.1\tGENEA\n' % tss_c('1') +
           '%s\t60500\t60500\tENST00000000002.1\tGENEB\n' % tss_c('2'))

    _write(os.path.join(root, 'tx_expr.tsv'),
           'ENST00000000001.1\t12.5\nENST00000000002.1\t0.75\n')
    _write(os.path.join(root, 'gene_list.tsv'),
           'x\tx\tx\tx\tGENEA\tENSG00000000001.5\n'
           'x\tx\tx\tx\tGENEB\tENSG00000000002.3\n')
    _write(os.path.join(root, 'gene_expr.tsv'),
           'ENSG00000000001.5\t33.0\nENSG00000000002.3\t4.5\n')

    for mark, chr1_enh, chr1_prom, chr2_enh, chr2_prom in (
            ('markX', 2.0, 3.0, 4.0, 5.0),
            ('markY', 1.0, 1.0, 1.0, 1.0)):
        _write(os.path.join(root, mark + '.bedGraph'),
               '%s\t9000\t11000\t%s\n' % (sig_c('1'), chr1_enh) +
               '%s\t21000\t22000\t%s\n' % (sig_c('1'), chr1_prom) +
               '%s\t49000\t51000\t%s\n' % (sig_c('2'), chr2_enh) +
               '%s\t60000\t61000\t%s\n' % (sig_c('2'), chr2_prom))

    _write(os.path.join(root, 'feature_list'),
           'MarkX\t%s\n' % os.path.join(root, 'markX.bedGraph') +
           'MarkY\t%s\n' % os.path.join(root, 'markY.bedGraph'))

    return {
        'input_gwas': os.path.join(root, 'snps.tsv'),
        'bedtools_path': 'bedtools',
        'hic_folder': hic,
        'hic_prefix': 'TEST.hic.KR.',
        'tss_file': os.path.join(root, 'tss.tsv'),
        'tx_expression': os.path.join(root, 'tx_expr.tsv'),
        'gene_list': os.path.join(root, 'gene_list.tsv'),
        'gene_expression': os.path.join(root, 'gene_expr.tsv'),
        'feature_list': os.path.join(root, 'feature_list'),
        'enhancer_window': 250,
        'promoter_window': 250,
        'hic_bin_size': 3000,
        'output_dir': os.path.join(root, 'out'),
        'output_file': 'test_features.tsv',
        'step1_nproc': 2,
    }


# One row per enhancer-promoter pair, every value hand-derived from the inputs.
EXPECTED = [
    {'#Enh_chr': '1', 'Enh_start': 9750, 'Enh_end': 10250,
     'Prom_start': 21250, 'Prom_end': 21750,
     'Prom_TXID': 'ENST00000000001.1', 'Prom_TSS': 21500, 'Prom_gene': 'GENEA',
     'SNPID_at_Enh': 'rs111', 'HiC_Contact': 5.5,
     'Tx_expression': 12.5, 'Gene_expression': 33.0,
     'MarkX_Enh': 1000.0, 'MarkX_Prom': 1500.0,
     'MarkY_Enh': 500.0, 'MarkY_Prom': 500.0},
    {'#Enh_chr': '2', 'Enh_start': 49750, 'Enh_end': 50250,
     'Prom_start': 60250, 'Prom_end': 60750,
     'Prom_TXID': 'ENST00000000002.1', 'Prom_TSS': 60500, 'Prom_gene': 'GENEB',
     'SNPID_at_Enh': 'rs222', 'HiC_Contact': 7.25,
     'Tx_expression': 0.75, 'Gene_expression': 4.5,
     'MarkX_Enh': 2000.0, 'MarkX_Prom': 2500.0,
     'MarkY_Enh': 500.0, 'MarkY_Prom': 500.0},
]
