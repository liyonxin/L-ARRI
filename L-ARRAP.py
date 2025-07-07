#!/usr/bin/env python3
"""
L-ARRI (Long-read Antibiotic Resistance Risk Index) Analysis Pipeline
Author: Yongxin Li
"""

import os
import sys
import argparse
import itertools
import subprocess

def create_directories():
    """Create all necessary directories"""
    dirs = [
        'arg-mge', 'Nanofilt_result', 'all_centrifuge', 'raw_mges',
        'extract_mges', 'mges_abundance', '3_hbp_abundance', 'ARRI',
        'ARG_raw_result', 'ARG_abundance', 'ARG_centrifuge', 'ARG_reads',
        'extract_ARG', 'arg_hbp'
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def get_sample_list(rawdata_path):
    """Get sample list and calculate sequence lengths"""
    samples = []
    for f in os.listdir(rawdata_path):
        if f.endswith('.fastq.gz'):
            sample_name = f.replace('.fastq.gz', '')
            samples.append(sample_name)
    
    sample_lengths = {}
    with open('sum_length.txt', 'w') as outfile0:
        for sample in samples:
            # Calculate sequence length
            cmd = f'seqkit stats {rawdata_path}/{sample}.fastq.gz -T -o {sample}_length.txt'
            subprocess.run(cmd, shell=True, check=True)
            
            with open(f'{sample}_length.txt', 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    length = lines[1].split('\t')[4].strip()
                    outfile0.write(f'{sample}\t{length}\n')
                    sample_lengths[sample] = float(length)
            
            os.remove(f'{sample}_length.txt')
    
    return samples, sample_lengths

def quality_control(sample, nanofilt_path, min_quality, min_length):
    """Perform quality control on raw reads"""
    raw_path = f"rawdata/{sample}.fastq.gz"
    output_gz = f"{nanofilt_path}/{sample}_nanofilt.gz"
    output_fa = f"{nanofilt_path}/{sample}.fa"
    
    # Quality control with Chopper
    cmd = f"gunzip -c {raw_path} | chopper -q {min_quality} -l {min_length} | gzip > {output_gz}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Convert to FASTA format
    cmd = f"seqtk seq -A {output_gz} > {output_fa}"
    subprocess.run(cmd, shell=True, check=True)

def run_centrifuge(sample, nanofilt_path, centrifuge_path, threads):
    """Run Centrifuge classification"""
    # Get database path
    base_dir = os.getcwd()
    centrifuge_db = os.path.join(base_dir, "database", "centrifuge", "p_compressed+h+v")
    
    input_fa = f"{nanofilt_path}/{sample}.fa"
    report_file = f"{centrifuge_path}/{sample}_all_report.tsv"
    result_file = f"{centrifuge_path}/{sample}_all_result.tsv"
    
    cmd = (f"centrifuge -f -x {centrifuge_db} -U {input_fa} "
           f"--report-file {report_file} -S {result_file} -p {threads}")
    subprocess.run(cmd, shell=True, check=True)

def run_minimap2_ARG(sample, nanofilt_path, arg_result_path, threads, preset):
    """Run Minimap2 alignment against ARG database"""
    # Get database path
    base_dir = os.getcwd()
    arg_db = os.path.join(base_dir, "database", "minimap2", "SARG_20211207_14210_filter.ffn")
    
    input_fa = f"{nanofilt_path}/{sample}.fa"
    output_file = f"{arg_result_path}/{sample}_ARG"
    
    # Run Minimap2 with specified preset
    cmd = (f"minimap2 -x {preset} -t {threads} "
           f"{arg_db} {input_fa} > {output_file}")
    subprocess.run(cmd, shell=True, check=True)

def run_last_MGE(sample, nanofilt_path, raw_mges_path, threads):
    """Run LAST alignment against MGE database"""
    # Get database path
    base_dir = os.getcwd()
    mge_db = os.path.join(base_dir, "database", "mobileOG", "mobileOG-db_beatrix-1.6.All.faa")
    
    input_fa = f"{nanofilt_path}/{sample}.fa"
    
    # Create LAST database
    cmd = f"lastdb -P{threads} -q -c trandb {mge_db}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Training and alignment
    train_file = f"{raw_mges_path}/{sample}_codon.train"
    maf_file = f"{raw_mges_path}/{sample}_out.maf"
    psl_file = f"{raw_mges_path}/{sample}_my-alignments.psl"
    
    cmd = f"last-train -P{threads} --codon trandb {input_fa} > {train_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    cmd = (f"lastal -P{threads} -p {train_file} -m100 -D1e9 -K1 trandb "
           f"{input_fa} > {maf_file}")
    subprocess.run(cmd, shell=True, check=True)
    
    cmd = f"maf-convert psl {maf_file} > {psl_file}"
    subprocess.run(cmd, shell=True, check=True)

def extract_ARG(sample, arg_result_path, extract_arg_path, min_coverage, min_identity):
    """Extract ARG results based on coverage and identity thresholds"""
    with open(f"{arg_result_path}/{sample}_ARG", 'r') as f:
        lines = f.readlines()
    
    with open(f"{extract_arg_path}/{sample}_ARG_cov_identity.txt", 'w') as out:
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 18:
                try:
                    coverage = (float(parts[8]) - float(parts[7])) / float(parts[6])
                    identity = float(parts[9]) / (float(parts[8]) - float(parts[7]))
                    
                    if coverage > min_coverage and identity > min_identity:
                        out.write(line)
                except (ValueError, ZeroDivisionError):
                    continue

def calculate_ARG_abundance(sample, extract_arg_path, abundance_path, sample_lengths):
    """Calculate ARG abundance"""
    k = float(sample_lengths[sample]) / 1e9
    
    with open(f"{extract_arg_path}/{sample}_ARG_cov_identity.txt", 'r') as f:
        lines = [line.strip().split('\t') for line in f.readlines()]
    
    gene_ids = set(parts[5] for parts in lines if len(parts) > 5)
    
    with open(f"{abundance_path}/{sample}_abundance_ARG.txt", 'w') as out:
        for gene_id in gene_ids:
            total_coverage = sum(
                (float(parts[8]) - float(parts[7])) / float(parts[6])
                for parts in lines if parts[5] == gene_id
            )
            abundance = total_coverage / k
            out.write(f"{gene_id}\t{abundance}\n")

def extract_MGE(sample, raw_mges_path, extract_mges_path, min_coverage, min_identity):
    """Extract MGE results with complex deduplication logic"""
    input_file = f"{raw_mges_path}/{sample}_my-alignments.psl"
    with open(input_file, 'r') as f:
        raw_lines = f.readlines()
    
    # Step 1: Basic filtering (coverage and identity)
    filtered_lines = []
    for line in raw_lines:
        parts = line.strip().split('\t')
        if len(parts) < 18:
            continue
            
        try:
            # Calculate coverage and identity
            target_length = float(parts[10])
            aln_length = float(parts[0])
            target_start = float(parts[11])
            target_end = float(parts[12])
            
            coverage = (target_end - target_start) / target_length
            identity = aln_length / (target_end - target_start)
            
            # Apply filtering thresholds
            if coverage > min_coverage and identity > min_identity:
                # Save complete line and parsed parts
                filtered_lines.append({
                    'raw': line,
                    'parts': parts,
                    'coverage': coverage,
                    'identity': identity,
                    'target_length': target_length,
                    'aln_length': aln_length,
                    'target_start': target_start,
                    'target_end': target_end,
                    'read_id': parts[13].strip()
                })
        except (ValueError, IndexError, ZeroDivisionError):
            continue
    
    # Step 2: Group by read
    read_groups = {}
    for item in filtered_lines:
        read_id = item['read_id']
        if read_id not in read_groups:
            read_groups[read_id] = []
        read_groups[read_id].append(item)
    
    # Step 3: Complex deduplication for multiple alignments per read
    final_results = []
    
    for read_id, alignments in read_groups.items():
        # If only one alignment, keep directly
        if len(alignments) == 1:
            final_results.append(alignments[0]['raw'])
            continue
        
        # Step 3.1: Identify highly overlapping alignment pairs
        highly_overlapping = set()
        combinations = list(itertools.combinations(alignments, 2))
        
        for align1, align2 in combinations:
            # Calculate overlap ratio
            start1, end1 = align1['target_start'], align1['target_end']
            start2, end2 = align2['target_start'], align2['target_end']
            
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            overlap_length = max(0, overlap_end - overlap_start)
            
            # Calculate overlap ratio relative to first alignment
            align_length = align1['target_end'] - align1['target_start']
            overlap_ratio = overlap_length / align_length if align_length > 0 else 0
            
            # Mark as highly overlapping if ratio > 80%
            if overlap_ratio > 0.8:
                highly_overlapping.add(tuple(align1['parts']))
                highly_overlapping.add(tuple(align2['parts']))
        
        # Step 3.2: Separate overlapping and non-overlapping alignments
        non_overlapping = []
        overlapping_group = []
        
        for align in alignments:
            if tuple(align['parts']) in highly_overlapping:
                overlapping_group.append(align)
            else:
                non_overlapping.append(align)
        
        # Step 3.3: Process non-overlapping alignments - keep directly
        for align in non_overlapping:
            final_results.append(align['raw'])
        
        # Step 3.4: Process overlapping group
        if not overlapping_group:
            continue
            
        # Sort by alignment start position
        overlapping_group.sort(key=lambda x: x['target_start'])
        
        # Step 3.5: Cluster nearby alignments
        clusters = []
        current_cluster = [overlapping_group[0]]
        
        for i in range(1, len(overlapping_group)):
            prev = current_cluster[-1]
            curr = overlapping_group[i]
            
            # Calculate distance between alignments
            distance = curr['target_start'] - prev['target_start']
            max_distance = prev['target_length'] * 0.4  # 40% of target length
            
            # If distance within threshold, add to current cluster
            if distance < max_distance:
                current_cluster.append(curr)
            else:
                # Save current cluster and start new one
                clusters.append(current_cluster)
                current_cluster = [curr]
        
        clusters.append(current_cluster)
        
        # Step 3.6: Process each cluster
        for cluster in clusters:
            # If cluster has only one alignment, keep directly
            if len(cluster) == 1:
                final_results.append(cluster[0]['raw'])
                continue
                
            # Find best alignment in cluster (highest identity)
            best_alignment = max(cluster, key=lambda x: x['identity'])
            final_results.append(best_alignment['raw'])
    
    # Save final results
    output_file = f"{extract_mges_path}/{sample}_same-reads_result.txt"
    with open(output_file, 'w') as out:
        for line in final_results:
            parts = line.strip().split('\t')
            # Keep only first 11 columns
            out.write('\t'.join(parts[:11]) + '\n')
    
    # Save intermediate filtered results (optional)
    filtered_file = f"{extract_mges_path}/{sample}_cov-idty.txt"
    with open(filtered_file, 'w') as out:
        for align in filtered_lines:
            out.write(align['raw'])
    
    print(f"Processed {len(filtered_lines)} alignments for {sample}, kept {len(final_results)} after deduplication")

def calculate_MGE_abundance(sample, extract_mges_path, abundance_path, sample_lengths):
    """Calculate MGE abundance"""
    k = float(sample_lengths[sample]) / 1e9
    
    with open(f"{extract_mges_path}/{sample}_same-reads_result.txt", 'r') as f:
        lines = [line.strip().split('\t') for line in f.readlines()]
    
    gene_ids = set(parts[2] for parts in lines if len(parts) > 2)
    
    with open(f"{abundance_path}/{sample}_abundance_mges.txt", 'w') as out:
        for gene_id in gene_ids:
            total_coverage = sum(
                (float(parts[5]) - float(parts[4])) / float(parts[3])
                for parts in lines if parts[2] == gene_id
            )
            abundance = total_coverage / k
            out.write(f"{gene_id}\t{abundance}\n")

def calculate_arg_mge(sample, extract_arg_path, extract_mges_path, arg_mge_path, sample_lengths):
    """Calculate ARG-MGE associations and abundances"""
    k = float(sample_lengths[sample]) / 1e9
    
    # Read ARG results
    with open(f"{extract_arg_path}/{sample}_ARG_cov_identity.txt", 'r') as f:
        arg_lines = [line.strip().split('\t') for line in f.readlines()]
    
    # Read MGE results
    with open(f"{extract_mges_path}/{sample}_same-reads_result.txt", 'r') as f:
        mge_lines = [line.strip().split('\t') for line in f.readlines()]
    
    # Find associations
    associations = []
    for arg in arg_lines:
        for mge in mge_lines:
            if arg[0] == mge[6]:
                associations.append((mge[2], arg[5], mge[6]))
    
    # Calculate abundances
    unique_reads = set(assoc[2] for assoc in associations)
    abundances = {}
    
    for read_id in unique_reads:
        coverage_sum = sum(
            (float(arg[8]) - float(arg[7])) / float(arg[6])
            for arg in arg_lines if arg[0] == read_id
        )
        abundances[read_id] = coverage_sum / k
    
    # Save results
    with open(f"{arg_mge_path}/{sample}_arg-mge-result.txt", 'w') as assoc_out:
        for assoc in associations:
            assoc_out.write(f"{assoc[0]}\t{assoc[1]}\t{assoc[2]}\n")
    
    with open(f"{arg_mge_path}/{sample}_arg-mge_abundance.txt", 'w') as abund_out:
        for read_id, abund in abundances.items():
            abund_out.write(f"{read_id}\t{abund}\n")

def calculate_HBP_abundance(sample, centrifuge_path, extract_arg_path, extract_mges_path, 
                            hbp_path, sample_lengths):
    """Calculate HBP-related abundances"""
    k = float(sample_lengths[sample]) / 1e9
    
    # Load HBP database
    base_dir = os.getcwd()
    hbp_file = os.path.join(base_dir, "database", "HBP", "new-species-all.txt")
    
    hbp_data = {}
    with open(hbp_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                hbp_data[parts[1]] = parts[0]  # taxid -> name
    
    # Load Centrifuge results
    with open(f"{centrifuge_path}/{sample}_all_result.tsv", 'r') as f:
        centrifuge_results = [line.strip().split('\t') for line in f if line.split('\t')[7] == '1']
    
    with open(f"{centrifuge_path}/{sample}_all_report.tsv", 'r') as f:
        centrifuge_report = {line.split('\t')[1]: line.split('\t') for line in f}
    
    # Calculate HBP abundance
    with open(f"{hbp_path}/{sample}_HBP_abundance.txt", 'w') as hbp_out:
        for taxid, name in hbp_data.items():
            if taxid in centrifuge_report:
                total_reads = float(centrifuge_report[taxid][3])
                matched_reads = sum(1 for res in centrifuge_results if res[2] == taxid)
                abundance = (matched_reads / total_reads) / k if total_reads > 0 else 0
                hbp_out.write(f"{name}\t{abundance}\n")
    
    # Calculate ARG-HBP abundance
    with open(f"{extract_arg_path}/{sample}_ARG_cov_identity.txt", 'r') as f:
        arg_data = [line.strip().split('\t') for line in f]
    
    with open(f"{hbp_path}/{sample}_ARG-HBP_abundance.txt", 'w') as arg_hbp_out:
        for res in centrifuge_results:
            read_id = res[0]
            coverage_sum = sum(
                (float(arg[8]) - float(arg[7])) / float(arg[6])
                for arg in arg_data if arg[0] == read_id
            )
            if coverage_sum > 0:
                abundance = coverage_sum / k
                arg_hbp_out.write(f"{read_id}\t{abundance}\n")
    
    # Calculate ARG-MGE-HBP abundance
    with open(f"{extract_mges_path}/{sample}_same-reads_result.txt", 'r') as f:
        mge_data = [line.strip().split('\t') for line in f]
    
    associated_reads = set()
    for arg in arg_data:
        for mge in mge_data:
            if arg[0] == mge[6]:
                associated_reads.add(arg[0])
    
    with open(f"{hbp_path}/{sample}_ARG-MGE-HBP_abundance.txt", 'w') as full_out:
        for read_id in associated_reads:
            coverage_sum = sum(
                (float(arg[8]) - float(arg[7])) / float(arg[6])
                for arg in arg_data if arg[0] == read_id
            )
            if coverage_sum > 0:
                abundance = coverage_sum / k
                full_out.write(f"{read_id}\t{abundance}\n")

def calculate_ARRI(sample, arg_abundance_path, mge_abundance_path, hbp_path, 
                   arg_mge_path, arri_path, sample_lengths):
    """Calculate Antibiotic Resistance Risk Index (ARRI)"""
    # Load abundance data
    def load_abundance(file_path):
        total = 0.0
        with open(file_path, 'r') as f:
            for line in f:
                total += float(line.split('\t')[1])
        return total
    
    arg_abund = load_abundance(f"{arg_abundance_path}/{sample}_abundance_ARG.txt")
    mge_abund = load_abundance(f"{mge_abundance_path}/{sample}_abundance_mges.txt")
    hbp_abund = load_abundance(f"{hbp_path}/{sample}_HBP_abundance.txt")
    arg_hbp_abund = load_abundance(f"{hbp_path}/{sample}_ARG-HBP_abundance.txt")
    arg_mge_hbp_abund = load_abundance(f"{hbp_path}/{sample}_ARG-MGE-HBP_abundance.txt")
    arg_mge_abund = load_abundance(f"{arg_mge_path}/{sample}_arg-mge_abundance.txt")
    
    # Calculate ARRI
    numerator = arg_mge_abund + arg_hbp_abund + arg_mge_hbp_abund
    denominator = arg_abund + mge_abund + hbp_abund
    arri = (numerator / denominator) * 10000 if denominator > 0 else 0
    
    # Save results
    with open(f"{arri_path}/{sample}_ARRI.txt", 'w') as f:
        f.write(f"ARG-MGE\t{arg_mge_abund}\n")
        f.write(f"ARG-HBP\t{arg_hbp_abund}\n")
        f.write(f"ARG-MGE-HBP\t{arg_mge_hbp_abund}\n")
        f.write(f"Abundance(ARGs)\t{arg_abund}\n")
        f.write(f"Abundance(MGEs)\t{mge_abund}\n")
        f.write(f"Abundance(HBPs)\t{hbp_abund}\n")
        f.write(f"ARRI\t{arri}\n")

def calculate_subARG_abundance(sample, extract_arg_path, arg_abundance_path, sample_lengths):
    """Calculate subARG/ARG category abundances"""
    k = float(sample_lengths[sample]) / 1e9
    
    # Load ARG structure database
    base_dir = os.getcwd()
    structure_file = os.path.join(base_dir, "database", "structure", "structure_20181107.LIST")
    
    arg_structure = {}
    with open(structure_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                key = parts[0]
                values = parts[1].replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')
                arg_structure[key] = values
    
    # Process ARG data
    with open(f"{extract_arg_path}/{sample}_ARG_cov_identity.txt", 'r') as f:
        arg_data = [line.strip().split('\t') for line in f]
    
    # Calculate subARG abundance
    subarg_counts = {}
    for arg in arg_data:
        read_id = arg[0]
        gene_id = arg[13].split(':')[0]
        coverage = (float(arg[8]) - float(arg[7])) / float(arg[6])
        
        for category, genes in arg_structure.items():
            if gene_id in genes:
                if category not in subarg_counts:
                    subarg_counts[category] = 0
                subarg_counts[category] += coverage
    
    # Calculate ARG category abundance
    arg_category_counts = {}
    for category in subarg_counts.keys():
        main_category = category.split('_')[0]
        if main_category not in arg_category_counts:
            arg_category_counts[main_category] = 0
        arg_category_counts[main_category] += subarg_counts[category]
    
    # Save results
    with open(f"{arg_abundance_path}/{sample}_subARG_abundance2.txt", 'w') as sub_out:
        for category, count in subarg_counts.items():
            abundance = count / k
            sub_out.write(f"{category}\t{abundance}\n")
    
    with open(f"{arg_abundance_path}/{sample}_ARG_abundance3.txt", 'w') as cat_out:
        for category, count in arg_category_counts.items():
            abundance = count / k
            cat_out.write(f"{category}\t{abundance}\n")

def main():
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description='L-ARRI Analysis Pipeline')
    parser.add_argument('-i', '--input', default='rawdata', 
                        help='Input directory with FASTQ files (default: rawdata)')
    parser.add_argument('-t', '--threads', type=int, default=56, 
                        help='Number of threads (default: 56)')
    parser.add_argument('-a_id', '--arg_identity', type=float, default=0.75, 
                        help='Minimum identity for ARG extraction (default: 0.75)')
    parser.add_argument('-a_cov', '--arg_coverage', type=float, default=0.9, 
                        help='Minimum coverage for ARG extraction (default: 0.9)')
    parser.add_argument('-m_id', '--mge_identity', type=float, default=0.75, 
                        help='Minimum identity for MGE extraction (default: 0.75)')
    parser.add_argument('-m_cov', '--mge_coverage', type=float, default=0.9, 
                        help='Minimum coverage for MGE extraction (default: 0.9)')
    parser.add_argument('-p', '--preset', choices=['map-ont', 'map-pb'], default='map-ont',
                        help='Minimap2 preset parameter (map-ont for Nanopore, map-pb for PacBio) (default: map-ont)')
    
    args = parser.parse_args()
    
    # Create directory structure
    create_directories()
    
    # Get sample list and sequence lengths
    samples, sample_lengths = get_sample_list(args.input)
    
    # Process each sample
    for sample in samples:
        print(f"\n{'='*50}")
        print(f"Processing sample: {sample}")
        print(f"{'='*50}")
        
        # Quality control
        print("[1/12] Performing quality control...")
        quality_control(sample, "Nanofilt_result", min_quality=10, min_length=500)
        
        # Centrifuge classification
        print("[2/12] Running Centrifuge classification...")
        run_centrifuge(sample, "Nanofilt_result", "all_centrifuge", args.threads)
        
        # ARG analysis
        print("[3/12] Identifying ARGs with Minimap2...")
        run_minimap2_ARG(sample, "Nanofilt_result", "ARG_raw_result", 
                         args.threads, args.preset)
        
        print("[4/12] Extracting ARG results...")
        extract_ARG(sample, "ARG_raw_result", "extract_ARG", 
                   min_coverage=args.arg_coverage, min_identity=args.arg_identity)
        
        print("[5/12] Calculating ARG abundance...")
        calculate_ARG_abundance(sample, "extract_ARG", "ARG_abundance", sample_lengths)
        
        # MGE analysis
        print("[6/12] Identifying MGEs with LAST...")
        run_last_MGE(sample, "Nanofilt_result", "raw_mges", args.threads)
        
        print("[7/12] Extracting and deduplicating MGE results...")
        extract_MGE(sample, "raw_mges", "extract_mges", 
                    min_coverage=args.mge_coverage, min_identity=args.mge_identity)
        
        print("[8/12] Calculating MGE abundance...")
        calculate_MGE_abundance(sample, "extract_mges", "mges_abundance", sample_lengths)
        
        # ARG-MGE association analysis
        print("[9/12] Analyzing ARG-MGE associations...")
        calculate_arg_mge(sample, "extract_ARG", "extract_mges", "arg-mge", sample_lengths)
        
        # HBP analysis
        print("[10/12] Calculating HBP abundances...")
        calculate_HBP_abundance(
            sample, "all_centrifuge", "extract_ARG", "extract_mges", 
            "3_hbp_abundance", sample_lengths
        )
        
        # ARRI calculation
        print("[11/12] Calculating ARRI...")
        calculate_ARRI(
            sample, "ARG_abundance", "mges_abundance", "3_hbp_abundance", 
            "arg-mge", "ARRI", sample_lengths
        )
        
        # subARG analysis
        print("[12/12] Calculating subARG abundances...")
        calculate_subARG_abundance(
            sample, "extract_ARG", "ARG_abundance", sample_lengths
        )
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main()