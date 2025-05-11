# task1_3_part1.py
#!/usr/bin/env python3
# Run DBG on reads_b.fastq with k=40

import os
import sys

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from genome_assembly_src.dbg_assembler import DBGAssembler

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
    input_file = "cs249_hw2/toy_dataset/reads_b.fastq"
    output_dir = "genome_assembly_output/dbg_results"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Run DBG assembly with k=40
    print("Running De Bruijn Graph assembly on reads_b.fastq with k=40...")
    assembler = DBGAssembler(k=40)
    fasta_file, gfa_file, contigs = assembler.assemble(input_file, output_dir)
    
    print("\nAssembly Statistics:")
    print(f"Number of contigs: {len(contigs)}")
    print(f"Total assembly length: {sum(len(c) for c in contigs)}")
    print(f"Largest contig: {max(len(c) for c in contigs) if contigs else 0}")
    print(f"N50: {calculate_n50(contigs)}")
    
    print("\nTo visualize the assembly graph, run:")
    print(f"Bandage load {gfa_file}")
    
    # Create a simple report file
    report_file = os.path.join(output_dir, "task1_part1_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 1: De Bruijn Graph Analysis\n\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"k-mer size: 40\n\n")
        f.write("## Assembly Statistics\n")
        f.write(f"Number of contigs: {len(contigs)}\n")
        f.write(f"Total assembly length: {sum(len(c) for c in contigs)}\n")
        f.write(f"Largest contig: {max(len(c) for c in contigs) if contigs else 0}\n")
        f.write(f"N50: {calculate_n50(contigs)}\n\n")
        f.write("## Assembly Graph Analysis\n")
        f.write("The assembly graph can be visualized using Bandage with the command:\n")
        f.write(f"Bandage load {gfa_file}\n\n")
        f.write("The graph structure shows how the reads connect to each other.\n")
        f.write("Look for any cycles, bubbles, or dead-ends in the graph,\n")
        f.write("which can indicate repetitive regions or sequencing errors.\n")
    
    print(f"\nReport written to {report_file}")

if __name__ == "__main__":
    main()
