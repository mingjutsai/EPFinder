#!/usr/bin/env python3
"""
EPFinder PreProcessing Workflow
"""

import os
import sys
import yaml
from pybedtools import BedTool
import subprocess

def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config["base_name"] = os.path.splitext(os.path.basename(config["input_gwas"]))[0]
    return config

def step1_hic_matrix_mapping(config):
    """Step 1: Map GWAS SNPs to Hi-C contact matrices."""
    input_file = config['input_gwas']
    hic_folder = config['hic_folder']
    hic_prefix = config['hic_prefix']
    output_file = os.path.join(config['tmp_dir'], config['base_name'] + ".hic_contact")

    print("Step 1: Mapping Hi-C contacts...")

    with open(output_file, 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line_num, line in enumerate(in_f, 1):
                line = line.strip()
                if line.startswith('#'):
                    continue
                print(f"No.{line_num}\t{line}")

                fields = line.split('\t')
                snp_chr = "chr" + fields[0]
                snp_start = int(fields[1])

                hic_file = os.path.join(hic_folder, hic_prefix + snp_chr)
                if not os.path.exists(hic_file):
                    print(f"Error: {hic_file} doesn't exist")
                    sys.exit(1)

                match = False
                with open(hic_file, 'r') as hic_f:
                    for hic_line in hic_f:
                        hic_fields = hic_line.strip().split('\t')
                        hic_bin1_start = int(hic_fields[0])
                        hic_bin1_end = hic_bin1_start + config['hic_bin_size']
                        hic_bin2_start = int(hic_fields[1])
                        hic_bin2_end = hic_bin2_start + config['hic_bin_size']
                        contact_freq = hic_fields[2]
                        if contact_freq == "nan":
                            contact_freq = "0"

                        if hic_bin1_start <= snp_start <= hic_bin1_end:
                            out_f.write(f"{line}\t{hic_bin2_start}\t{hic_bin2_end}\t{contact_freq}\n")
                            match = True
                        elif hic_bin2_start <= snp_start <= hic_bin2_end:
                            out_f.write(f"{line}\t{hic_bin1_start}\t{hic_bin1_end}\t{contact_freq}\n")
                            match = True

                if not match:
                    out_f.write(f"{line}\tnan\tnan\t0\n")

    print(f"Step 1 completed. Output: {output_file}")

def step2_hic_prom(config):
    """Step 2: Filter and format Hi-C contact data."""
    input_file = os.path.join(config['tmp_dir'], config['base_name'] + ".hic_contact")
    output_file = input_file + "_hicbin2"

    print("Step 2: Filtering Hi-C data...")

    with open(output_file, 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                hic_bin_start = fields[3]
                hic_bin_end = fields[4]
                if hic_bin_start != "nan" and hic_bin_end != "nan":
                    chr_name = fields[0]
                    out_f.write(f"{chr_name}\t{hic_bin_start}\t{hic_bin_end}\t{line}\n")

    # Sort lexicographically
    sorted_file = output_file + "_sorted_lexicographical"
    subprocess.run(f"sort -k1,1 -k2,2n {output_file} > {sorted_file}", shell=True)

    print(f"Step 2 completed. Output: {sorted_file}")

def step3_hicbin2_tss_overlap(config):
    """Step 3: Find overlaps between Hi-C bins and TSS."""
    input_file = os.path.join(config['tmp_dir'], config['base_name'] + ".hic_contact_hicbin2_sorted_lexicographical")
    tss_file = config['tss_file']
    output_file = os.path.join(config['tmp_dir'], f"{config['base_name']}.hic_contact_hicbin2_tss")

    print("Step 3: Finding TSS overlaps...")

    a = BedTool(input_file)
    b = BedTool(tss_file)
    result = a.intersect(b, wa=True, wb=True, sorted=True, F=1.0)
    result.saveas(output_file)

    print(f"Step 3 completed. Output: {output_file}")

def step4_format(config):
    """Step 4: Format the data with enhancer and promoter regions."""
    input_file = os.path.join(config['tmp_dir'], f"{config['base_name']}.hic_contact_hicbin2_tss")
    output_file = input_file + ".format"
    enhancer_window = config['enhancer_window']
    promoter_window = config['promoter_window']

    print("Step 4: Formatting data...")

    with open(output_file, 'w') as out_f:
        out_f.write("#Enh_chr\tEnh_start\tEnh_end\tProm_start\tProm_end\tProm_TXID\tProm_TSS\tProm_gene\tSNPID_at_Enh\tHiC_Contact\n")

        with open(input_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                enh_chr = fields[3]
                snp_pos = int(fields[4])
                enh_start = snp_pos - enhancer_window
                enh_end = snp_pos + enhancer_window
                tss_pos = int(fields[10])
                prom_start = tss_pos - promoter_window
                prom_end = tss_pos + promoter_window
                prom_txid = fields[12]
                prom_gene = fields[13]
                hic_contact = fields[8]
                snpid = fields[5]

                info = f"{enh_chr}\t{enh_start}\t{enh_end}\t{prom_start}\t{prom_end}\t{prom_txid}\t{tss_pos}\t{prom_gene}\t{snpid}\t{hic_contact}\n"
                out_f.write(info)

    print(f"Step 4 completed. Output: {output_file}")

def step5_mapping_tx_quantifications(config):
    """Step 5: Map transcript quantifications."""
    input_file = os.path.join(config['tmp_dir'], f"{config['base_name']}.hic_contact_hicbin2_tss.format")
    tx_file = config['tx_expression']
    output_file = input_file + ".Tx_expression"

    print("Step 5: Mapping transcript expressions...")

    # Load tx expressions
    tx_dict = {}
    with open(tx_file, 'r') as tx_f:
        for line in tx_f:
            line = line.strip()
            fields = line.split('\t')
            if fields[0].startswith('ENST'):
                tx_id = fields[0].split('.')[0]
                tx_dict[tx_id] = fields[1]

    with open(output_file, 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line.startswith('#'):
                    out_f.write(line + "\tTx_expression\n")
                    continue

                fields = line.split('\t')
                txid_full = fields[5]
                txid = txid_full.split('.')[0]
                expr = tx_dict.get(txid, '0')
                out_f.write(line + f"\t{expr}\n")
                if expr == '0' and txid in tx_dict:
                    print(f"Unmapping ID: {txid}\t{line}")

    print(f"Step 5 completed. Output: {output_file}")

def step6_mapping_gene_quantification(config):
    """Step 6: Map gene quantifications."""
    input_file = os.path.join(config['tmp_dir'], f"{config['base_name']}.hic_contact_hicbin2_tss.format.Tx_expression")
    gene_list_file = config['gene_list']
    gene_expr_file = config['gene_expression']
    output_file = input_file + ".GENEexpression"

    print("Step 6: Mapping gene expressions...")

    # Load gene mappings
    gene_dict = {}
    with open(gene_list_file, 'r') as gl_f:
        for line in gl_f:
            line = line.strip()
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            gene_dict[fields[4]] = fields[5]

    # Load gene expressions
    expr_dict = {}
    with open(gene_expr_file, 'r') as ge_f:
        for line in ge_f:
            line = line.strip()
            if line.startswith('ENSG'):
                fields = line.split('\t')
                expr_dict[fields[0]] = fields[1]

    with open(output_file, 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line.startswith('#'):
                    out_f.write(line + "\tGene_expression\n")
                    continue

                fields = line.split('\t')
                gene = fields[7]
                if gene == "C19orf43":
                    gene = "TRIR"

                if gene in gene_dict:
                    gene_id = gene_dict[gene]
                    if gene_id in expr_dict:
                        expr = expr_dict[gene_id]
                        out_f.write(line + f"\t{expr}\n")
                    else:
                        print(f"Error: No expression for gene_id {gene_id}")
                        sys.exit(1)
                else:
                    print(f"Error: Gene {gene} not found")
                    sys.exit(1)

    print(f"Step 6 completed. Output: {output_file}")

def step7_split_class_enh_prom(config):
    """Step 7: Split into enhancer and promoter BED files."""
    input_file = os.path.join(config['tmp_dir'], f"{config['base_name']}.hic_contact_hicbin2_tss.format.Tx_expression.GENEexpression")
    enh_file = input_file + ".enh"
    prom_file = input_file + ".prom"

    print("Step 7: Splitting into enhancer and promoter files...")

    with open(enh_file, 'w') as enh_f, open(prom_file, 'w') as prom_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                info = ','.join(fields)
                enh_f.write(f"chr{fields[0]}\t{fields[1]}\t{fields[2]}\t{info}\n")
                prom_f.write(f"chr{fields[0]}\t{fields[3]}\t{fields[4]}\t{info}\n")

    print(f"Step 7 completed. Outputs: {enh_file}, {prom_file}")

def step8_bedtool_overlap(config):
    """Step 8: Run bedtools overlap for each feature."""
    input_base = f"{config['base_name']}.hic_contact_hicbin2_tss.format.Tx_expression.GENEexpression"
    enh_file = os.path.join(config['tmp_dir'], input_base + ".enh")
    prom_file = os.path.join(config['tmp_dir'], input_base + ".prom")
    feature_list = config['feature_list']
    # Sort the files
    enh_sorted = enh_file + ".sorted"
    BedTool(enh_file).sort().saveas(enh_sorted)
    enh_file = enh_sorted
    prom_sorted = prom_file + ".sorted"
    BedTool(prom_file).sort().saveas(prom_sorted)
    prom_file = prom_sorted

    print("Step 8: Running bedtools overlap for features...")

    with open(feature_list, 'r') as fl_f:
        for line in fl_f:
            fields = line.strip().split('	')
            mark = fields[0]
            target_bed = fields[1]
            print(f"Processing {mark}...")

            if not os.path.exists(os.path.join(config['tmp_dir'], mark)):
                os.makedirs(os.path.join(config['tmp_dir'], mark))

            enh_re = os.path.join(config['tmp_dir'], mark, os.path.basename(enh_file) + "." + mark)
            if not os.path.exists(enh_re):
                BedTool(enh_file).intersect(BedTool(target_bed), wao=True, sorted=True).saveas(enh_re)

            enh_merge = enh_re + ".merge"
            if not os.path.exists(enh_merge):
                merge_overlaps(enh_re, enh_merge)

            prom_re = os.path.join(config['tmp_dir'], mark, os.path.basename(prom_file) + "." + mark)
            if not os.path.exists(prom_re):
                BedTool(prom_file).intersect(BedTool(target_bed), wao=True, sorted=True).saveas(prom_re)

            prom_merge = prom_re + ".merge"
            if not os.path.exists(prom_merge):
                merge_overlaps(prom_re, prom_merge)

    print("Step 8 completed.")

def merge_overlaps(input_file, output_file):
    """Merge overlap results by summing value * bp."""
    df = pd.read_csv(input_file, sep='\t', header=None, names=['chr_a', 'start_a', 'end_a', 'id', 'chr_b', 'start_b', 'end_b', 'value', 'bp'])
    df['value'] = pd.to_numeric(df['value'], errors='coerce').fillna(0)
    df['bp'] = pd.to_numeric(df['bp'], errors='coerce').fillna(0)
    df['area'] = df['value'] * df['bp']
    merged = df.groupby('id')['area'].sum().reset_index()
    merged.to_csv(output_file, sep='\t', header=False, index=False)

def add_feature(input_file, feature, fix_in, config):
    print(f"input_file:{input_file}\n")
    base_name = os.path.basename(fix_in)
    enh_merge = os.path.join(config["tmp_dir"], feature, base_name + ".enh.sorted." + feature + ".merge")
    prom_merge = os.path.join(config["tmp_dir"], feature, base_name + ".prom.sorted." + feature + ".merge")
    
    enh_values = {}
    if os.path.exists(enh_merge):
        with open(enh_merge, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("	")
                    if len(parts) == 2:
                        #print(f"enh part[0]:{parts[0]}")
                        enh_values[parts[0]] = parts[1]
    
    prom_values = {}
    if os.path.exists(prom_merge):
        with open(prom_merge, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("	")
                    if len(parts) == 2:
                        #print(f"prom parts[0]:{parts[0]}")
                        prom_values[parts[0]] = parts[1]
    
    output_file = input_file + "_" + feature
    with open(output_file, "w") as out_f:
        with open(input_file, "r") as in_f:
            count = 0
            for line in in_f:
                count += 1
                line = line.strip()
                if line.startswith("#"):
                    out_f.write(line + "	" + feature + "_Enh	" + feature + "_Prom\n")
                    continue
                fields = line.split("	")
                info = ",".join(fields[0:12])
                #print(f"info:{info}")
                enh_val = enh_values.get(info)
                if enh_val is None:
                    raise KeyError(f"Missing key in enh_value: {info}")
                prom_val = prom_values.get(info)
                if prom_val is None:
                    raise KeyError(f"Missing key in prom_values: {info}")
                out_f.write(line + "	" + str(enh_val) + "	" + str(prom_val) + "\n")

def step9_add_allfeatures(config):
    """Step 9: Add all features to the main file."""
    base_input = os.path.join(config["tmp_dir"], f"{config['base_name']}.hic_contact_hicbin2_tss.format.Tx_expression.GENEexpression")
    features_sort_file = config["features_sort"]

    print("Step 9: Adding all features...")

    current_input = None
    with open(features_sort_file, "r") as fs_f:
        for line in fs_f:
            feature = line.strip()
            print(f"Adding {feature}...")

            if current_input is None:
                add_feature(base_input, feature, base_input, config)
                current_input = f"{base_input}_{feature}"
            else:
                add_feature(current_input, feature, base_input, config)
                current_input += f"_{feature}"

    print("Step 9 completed.")
    return current_input
    

def main():
    """Main function to run the entire workflow."""
    if len(sys.argv) != 2:
        print("Usage: python EPFinder_preprocessing.py config.yaml")
        sys.exit(1)

    config_file = sys.argv[1]
    config = load_config(config_file)

    # Change to output directory
    os.makedirs(config["output_dir"], exist_ok=True)
    os.chdir(config['output_dir'])
    config['tmp_dir'] = os.path.abspath("tmp/") + "/"
    os.makedirs(config['tmp_dir'], exist_ok=True)

    config["config['tmp_dir']"] = config['tmp_dir']
    #Run all steps
    step1_hic_matrix_mapping(config)
    step2_hic_prom(config)
    step3_hicbin2_tss_overlap(config)
    step4_format(config)
    step5_mapping_tx_quantifications(config)
    step6_mapping_gene_quantification(config)
    step7_split_class_enh_prom(config)
    step8_bedtool_overlap(config)
    final_file = step9_add_allfeatures(config)
    if "output_file" in config:
        os.rename(final_file, config["output_file"])
        print(f"Renamed final output to {config["output_file"]}")

    print("Workflow completed successfully!")

if __name__ == "__main__":
    main()
