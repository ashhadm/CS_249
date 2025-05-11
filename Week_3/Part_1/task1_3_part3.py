# task1_3_part3.py
#!/usr/bin/env python3
# Run DBG and OLC on MERS datasets

import os
import sys
import pandas as pd

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from genome_assembly_src.dbg_assembler import DBGAssembler
from genome_assembly_src.olc_assembler import OLCAssembler
from genome_assembly_src.evaluation import run_quast, parse_quast_report, compare_assemblies

def main():
    # Set up paths
    reference_file = "cs249_hw2/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna"
    output_dir_dbg = "genome_assembly_output/dbg_results"
    output_dir_olc = "genome_assembly_output/olc_results"
    quast_dir = "genome_assembly_output/quast_results"
    
    # Create output directories if they don't exist
    os.makedirs(output_dir_dbg, exist_ok=True)
    os.makedirs(output_dir_olc, exist_ok=True)
    os.makedirs(quast_dir, exist_ok=True)
    
    # Dataset paths
    datasets = [
        "cs249_hw2/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq",
        "cs249_hw2/synthetic_dataset/reads/no_error_ont_hq_50x.fastq",
        "cs249_hw2/synthetic_dataset/reads/reads_hiseq_5k.fastq",
        "cs249_hw2/synthetic_dataset/reads/ont_hq_50x.fastq"
    ]
    
    # Parameters
    dbg_k = 31  # k-mer size for DBG
    olc_min_overlap = 20  # minimum overlap length for OLC
    
    reports = []
    results = []
    
    # Process each dataset
    for dataset in datasets:
        basename = os.path.basename(dataset).split('.')[0]
        print(f"\nProcessing dataset: {basename}")
        
        # Run DBG assembly
        print(f"Running De Bruijn Graph assembly with k={dbg_k}...")
        assembler_dbg = DBGAssembler(k=dbg_k)
        fasta_dbg, _, contigs_dbg = assembler_dbg.assemble(dataset, output_dir_dbg)
        
        # Run OLC assembly
        print(f"Running Overlap-Layout-Consensus assembly with min_overlap={olc_min_overlap}...")
        assembler_olc = OLCAssembler(min_overlap=olc_min_overlap)
        fasta_olc, contigs_olc = assembler_olc.assemble(dataset, output_dir_olc)
        
        # Run QUAST for both assemblies
        print("Evaluating assemblies with QUAST...")
        quast_dbg_dir = os.path.join(quast_dir, f"{basename}_dbg")
        quast_olc_dir = os.path.join(quast_dir, f"{basename}_olc")
        
        quast_dbg_report = run_quast(fasta_dbg, reference_file, quast_dbg_dir)
        quast_olc_report = run_quast(fasta_olc, reference_file, quast_olc_dir)
        
        # Add to reports list
        reports.append((f"{basename}_dbg", quast_dbg_report))
        reports.append((f"{basename}_olc", quast_olc_report))
        
        # Collect basic statistics
        results.append({
            "Dataset": basename,
            "Algorithm": "DBG",
            "Contigs": len(contigs_dbg),
            "Total Length": sum(len(c) for c in contigs_dbg),
            "Largest Contig": max(len(c) for c in contigs_dbg) if contigs_dbg else 0
        })
        
        results.append({
            "Dataset": basename,
            "Algorithm": "OLC",
            "Contigs": len(contigs_olc),
            "Total Length": sum(len(c) for c in contigs_olc),
            "Largest Contig": max(len(c) for c in contigs_olc) if contigs_olc else 0
        })
    
    # Compare all assemblies
    print("\nComparing all assemblies...")
    comparison_file = os.path.join(quast_dir, "all_assemblies_comparison.tsv")
    comparison = compare_assemblies(reports, comparison_file)
    
    print(f"\nFull comparison saved to {comparison_file}")
    
    # Create a summary report
    results_df = pd.DataFrame(results)
    results_file = os.path.join(quast_dir, "assembly_results_summary.tsv")
    results_df.to_csv(results_file, sep='\t', index=False)
    
    print(f"\nSummary report saved to {results_file}")
    
    # Create a simple report file
    report_file = os.path.join(quast_dir, "task1_part3_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 3: Comparing DBG and OLC Algorithms on MERS Datasets\n\n")
        f.write(f"Reference genome: {reference_file}\n\n")
        
        f.write("## Assembly Statistics\n")
        f.write(results_df.to_string(index=False))
        f.write("\n\n")
        
        f.write("## Analysis\n")
        f.write("This analysis compares the De Bruijn Graph (DBG) and Overlap-Layout-Consensus (OLC)\n")
        f.write("algorithms on both error-free and error-containing datasets.\n\n")
        
        f.write("### Key Observations:\n")
        f.write("1. Error Impact: Both algorithms perform better on error-free reads\n")
        f.write("2. Read Type: Performance varies between short (Illumina HiSeq) and long (ONT) reads\n")
        f.write("3. Algorithm Strengths:\n")
        f.write("   - DBG typically handles short reads well but is sensitive to errors\n")
        f.write("   - OLC can work better with longer reads and may be more resilient to some errors\n\n")
        
        f.write("For detailed metrics, refer to the QUAST reports in each respective directory.\n")
        f.write(f"The full comparison is available at: {comparison_file}\n")
    
    print(f"\nReport written to {report_file}")

if __name__ == "__main__":
    main()
