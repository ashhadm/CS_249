# task1_3_part2.py
#!/usr/bin/env python3
# Run DBG on reads_r.fastq with k=35 and k=45

import os
import sys
import pandas as pd

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from genome_assembly_src.dbg_assembler import DBGAssembler
from genome_assembly_src.evaluation import run_quast, parse_quast_report, compare_assemblies

def calculate_n50(contigs):
    """Calculate N50 of contigs."""
    contig_lengths = sorted([len(c) for c in contigs], reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length / 2
    
    current_length = 0
    for length in contig_lengths:
        current_length += length
        if current_length >= target_length:
            return length
    
    return 0

def main():
    # Set up paths
    input_file = "cs249_hw2/toy_dataset/reads_r.fastq"
    reference_file = "cs249_hw2/toy_dataset/reference_r.fasta"
    output_dir = "genome_assembly_output/dbg_results"
    quast_dir = "genome_assembly_output/quast_results"
    
    # Create output directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(quast_dir, exist_ok=True)
    
    # Run DBG assembly with k=35
    print("Running De Bruijn Graph assembly on reads_r.fastq with k=35...")
    assembler_k35 = DBGAssembler(k=35)
    fasta_k35, gfa_k35, contigs_k35 = assembler_k35.assemble(input_file, output_dir)
    
    # Run DBG assembly with k=45
    print("\nRunning De Bruijn Graph assembly on reads_r.fastq with k=45...")
    assembler_k45 = DBGAssembler(k=45)
    fasta_k45, gfa_k45, contigs_k45 = assembler_k45.assemble(input_file, output_dir)
    
    # Run QUAST for both assemblies
    print("\nEvaluating assemblies with QUAST...")
    quast_k35_dir = os.path.join(quast_dir, "reads_r_k35")
    quast_k45_dir = os.path.join(quast_dir, "reads_r_k45")
    
    quast_k35_report = run_quast(fasta_k35, reference_file, quast_k35_dir)
    quast_k45_report = run_quast(fasta_k45, reference_file, quast_k45_dir)
    
    # Compare assemblies
    print("\nComparing assemblies...")
    reports = [
        ("k=35", quast_k35_report),
        ("k=45", quast_k45_report)
    ]
    
    comparison_file = os.path.join(quast_dir, "k35_vs_k45_comparison.tsv")
    comparison = compare_assemblies(reports, comparison_file)
    
    # Print basic statistics
    print("\nAssembly Comparison:")
    print(f"K=35: {len(contigs_k35)} contigs, total length: {sum(len(c) for c in contigs_k35)}, N50: {calculate_n50(contigs_k35)}")
    print(f"K=45: {len(contigs_k45)} contigs, total length: {sum(len(c) for c in contigs_k45)}, N50: {calculate_n50(contigs_k45)}")
    
    print("\nTo visualize the assembly graphs, run:")
    print(f"Bandage load {gfa_k35}")
    print(f"Bandage load {gfa_k45}")
    
    print(f"\nFull comparison saved to {comparison_file}")
    
    # Create a simple report file
    report_file = os.path.join(output_dir, "task1_part2_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 2: Effect of k-mer Size on De Bruijn Graph Assembly\n\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"Reference file: {reference_file}\n\n")
        
        f.write("## Assembly Statistics\n")
        f.write(f"K=35: {len(contigs_k35)} contigs, total length: {sum(len(c) for c in contigs_k35)}, N50: {calculate_n50(contigs_k35)}\n")
        f.write(f"K=45: {len(contigs_k45)} contigs, total length: {sum(len(c) for c in contigs_k45)}, N50: {calculate_n50(contigs_k45)}\n\n")
        
        f.write("## Analysis\n")
        f.write("The choice of k-mer size affects the assembly in several ways:\n\n")
        f.write("1. Smaller k (k=35): Potentially more connections between reads but more ambiguity in repeat regions\n")
        f.write("2. Larger k (k=45): Better resolution of repeats but may lead to fragmentation due to sequencing errors\n\n")
        
        f.write("To visualize the assembly graphs, run:\n")
        f.write(f"Bandage load {gfa_k35}\n")
        f.write(f"Bandage load {gfa_k45}\n\n")
        
        f.write("For the full comparison of quality metrics, see the QUAST report at:\n")
        f.write(f"{comparison_file}\n")
    
    print(f"\nReport written to {report_file}")

if __name__ == "__main__":
    main()
