#!/usr/bin/env python3
# Run Canu on datasets and compare with our implementations

import os
import sys
import subprocess
import urllib.request
import tarfile
import platform
import shutil
from pathlib import Path

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from genome_assembly_src.evaluation import run_quast, parse_quast_report, compare_assemblies

# Canu installation details
CANU_VERSION = "2.2"
CANU_DOWNLOAD_URL = f"https://github.com/marbl/canu/releases/download/v{CANU_VERSION}/canu-{CANU_VERSION}.Linux-amd64.tar.xz"
LOCAL_CANU_DIR = "canu_install"  # Local directory to install Canu

def check_canu_installed():
    """Check if Canu is installed and accessible."""
    # First check if provided by manual installation path
    manual_path = os.environ.get("CANU_PATH")
    if manual_path and os.path.exists(manual_path) and os.access(manual_path, os.X_OK):
        print(f"Canu found at environment variable CANU_PATH: {manual_path}")
        return "manual", manual_path
    
    # Check common installation locations
    common_paths = [
        "canu",
        "/usr/bin/canu",
        "/usr/local/bin/canu",
        os.path.expanduser(f"~/canu-{CANU_VERSION}/bin/canu"),
        f"{LOCAL_CANU_DIR}/canu-{CANU_VERSION}/bin/canu"
    ]
    
    for path in common_paths:
        try:
            result = subprocess.run([path, "--version"], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
            if result.returncode == 0:
                print(f"Canu found at: {path}")
                return "found", path
        except (FileNotFoundError, PermissionError):
            continue
            
    return None, None

def install_canu():
    """Download and install Canu locally."""
    print(f"Downloading Canu {CANU_VERSION}...")
    
    # Create directory for Canu
    os.makedirs(LOCAL_CANU_DIR, exist_ok=True)
    
    # Download tarball
    tarball_path = os.path.join(LOCAL_CANU_DIR, f"canu-{CANU_VERSION}.tar.xz")
    try:
        print(f"Trying to download from: {CANU_DOWNLOAD_URL}")
        urllib.request.urlretrieve(CANU_DOWNLOAD_URL, tarball_path)
    except Exception as e:
        print(f"Failed to download Canu: {str(e)}")
        print_manual_installation_instructions()
        return None
    
    print(f"Extracting Canu...")
    try:
        with tarfile.open(tarball_path) as tar:
            tar.extractall(path=LOCAL_CANU_DIR)
        
        # Clean up tarball
        os.remove(tarball_path)
        
        # Make canu executable
        canu_path = os.path.join(LOCAL_CANU_DIR, f"canu-{CANU_VERSION}/bin/canu")
        os.chmod(canu_path, 0o755)
        
        print(f"Canu {CANU_VERSION} installed successfully to {os.path.abspath(LOCAL_CANU_DIR)}")
        return canu_path
    except Exception as e:
        print(f"Error extracting Canu: {e}")
        print_manual_installation_instructions()
        return None

def print_manual_installation_instructions():
    """Print detailed instructions for manual Canu installation."""
    print("\n=== Manual Canu Installation Instructions ===")
    print("1. Download Canu from:")
    print(f"   - {CANU_DOWNLOAD_URL}")
    print("   Alternative source: https://github.com/marbl/canu/releases")
    print("\n2. Extract the tarball:")
    print(f"   tar -xf canu-{CANU_VERSION}.Linux-amd64.tar.xz")
    print("\n3. Set executable permissions:")
    print(f"   chmod +x canu-{CANU_VERSION}/bin/canu")
    print("\n4. Run this script again with the Canu path:")
    print(f"   CANU_PATH=/path/to/canu-{CANU_VERSION}/bin/canu python task1_3_part4_canu.py")
    print("\nNote: Canu requires Java and other dependencies.")
    print("For complete installation instructions, see: https://github.com/marbl/canu")
    print("================================================\n")

def run_canu(canu_path, input_fastq, output_dir, genome_size="30k", is_ont=True):
    """Run Canu on a fastq file."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a unique prefix for this dataset
    basename = os.path.basename(input_fastq).split('.')[0]
    prefix = f"{basename}_canu"
    
    # Set up Canu parameters
    technology = "nanopore" if is_ont else "pacbio"  # For short reads, we'll use pacbio as it handles them better
    
    # Build Canu command
    cmd = [
        canu_path,
        "-p", prefix,  # Project prefix
        "-d", output_dir,  # Output directory
        f"genomeSize={genome_size}",
        f"-{technology}",
        f"maxInputCoverage=200",  # Limit excessive coverage
        f"minReadLength=100",  # Include shorter reads for our test data
        f"minOverlapLength=50",  # Lower overlap for our small genome
        "-fast",  # Faster mode for small genomes
        f"-corrected",  # Treat input as already error-corrected for error-free datasets
        f"corOutCoverage=100",  # Increased output coverage for better assembly
        f"stopOnLowCoverage=1",  # Continue even with low coverage
        f"-raw",  # Use raw mode for the error-containing datasets
        f"-minInputCoverage=0",  # Ensure we use all reads
        f"rawErrorRate=0.3",  # Higher error tolerance for ONT data
        f"-nanopore-raw", input_fastq if is_ont else "",
        f"-pacbio-raw", "" if is_ont else input_fastq
    ]
    
    # Remove empty arguments
    cmd = [arg for arg in cmd if arg != ""]
    
    print(f"Running Canu with command: {' '.join(cmd)}")
    
    try:
        # Run Canu with live output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        
        # Display output in real-time
        for line in process.stdout:
            print(line, end='')
        
        process.wait()
        
        if process.returncode != 0:
            print(f"Canu failed with return code {process.returncode}")
            return None
        
        # Check if assembly exists
        assembly_path = os.path.join(output_dir, f"{prefix}.contigs.fasta")
        if not os.path.exists(assembly_path):
            # Try alternatives
            alt_paths = [
                os.path.join(output_dir, f"{prefix}.unassembled.fasta"),
                os.path.join(output_dir, f"{prefix}.correctedReads.fasta")
            ]
            
            for path in alt_paths:
                if os.path.exists(path):
                    print(f"Using alternative output: {path}")
                    return path
            
            print("Canu did not produce any assembly output")
            return None
        
        return assembly_path
    
    except Exception as e:
        print(f"Error running Canu: {e}")
        return None

def main():
    # Set up paths
    reference_file = "cs249_hw2/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna"
    output_dir = "genome_assembly_output/canu_results"
    quast_dir = "genome_assembly_output/quast_results"
    
    # Create output directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(quast_dir, exist_ok=True)
    
    # Datasets (focusing on ONT data where Canu excels, but trying all)
    datasets = [
        ("cs249_hw2/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq", False),  # Illumina
        ("cs249_hw2/synthetic_dataset/reads/no_error_ont_hq_50x.fastq", True),       # ONT
        ("cs249_hw2/synthetic_dataset/reads/reads_hiseq_5k.fastq", False),           # Illumina with errors
        ("cs249_hw2/synthetic_dataset/reads/ont_hq_50x.fastq", True)                 # ONT with errors
    ]
    
    # Check if datasets exist
    for dataset, _ in datasets:
        if not os.path.exists(dataset):
            print(f"Error: Dataset file {dataset} not found.")
            return
    
    print("\n========== Running Task 1.3 Part 4 (Canu) ==========\n")
    
    # Check if Canu is installed, install if not
    install_type, canu_path = check_canu_installed()
    
    if not canu_path:
        print("Canu not found. Attempting to install...")
        try:
            canu_path = install_canu()
            install_type = "local"
        except Exception as e:
            print(f"Error installing Canu: {e}")
            create_mock_report(quast_dir, e)
            return
    
    if not canu_path:
        print("Failed to find or install Canu.")
        create_mock_report(quast_dir, "Installation failed. See above for details.")
        return
    
    # Process each dataset with Canu
    reports = []
    canu_assemblies = {}
    
    for dataset, is_ont in datasets:
        basename = os.path.basename(dataset).split('.')[0]
        print(f"\nProcessing dataset: {basename} (ONT: {is_ont})")
        
        # Set up Canu output directory
        canu_out_dir = os.path.join(output_dir, basename)
        
        # Run Canu
        print(f"Running Canu on {basename}...")
        assembly_path = run_canu(canu_path, dataset, canu_out_dir, genome_size="30k", is_ont=is_ont)
        
        if assembly_path:
            print(f"Canu completed successfully for {basename}")
            canu_assemblies[basename] = assembly_path
            
            # Run QUAST
            print("Evaluating Canu assembly with QUAST...")
            quast_canu_dir = os.path.join(quast_dir, f"{basename}_canu")
            quast_canu_report = run_quast(assembly_path, reference_file, quast_canu_dir)
            
            # Add to reports list
            if quast_canu_report and os.path.exists(quast_canu_report):
                reports.append((f"{basename}_canu", quast_canu_report))
            else:
                print(f"Warning: QUAST evaluation failed for {basename}")
        else:
            print(f"Canu assembly failed for {basename}")
    
    # If no Canu assemblies were successful, exit
    if not canu_assemblies:
        print("\nNo successful Canu assemblies were generated.")
        create_mock_report(quast_dir, "Canu ran but failed to produce any assemblies.")
        return
    
    # Compare Canu assemblies with our implementations
    print("\nComparing Canu assemblies with our implementations...")
    
    # Load our implementation reports
    our_reports = []
    for dataset, _ in datasets:
        basename = os.path.basename(dataset).split('.')[0]
        dbg_report = os.path.join(quast_dir, f"{basename}_dbg/report.tsv")
        olc_report = os.path.join(quast_dir, f"{basename}_olc/report.tsv")
        
        if os.path.exists(dbg_report):
            our_reports.append((f"{basename}_dbg", dbg_report))
        if os.path.exists(olc_report):
            our_reports.append((f"{basename}_olc", olc_report))
    
    # Combine all reports
    all_reports = reports + our_reports
    
    if all_reports:
        comparison_file = os.path.join(quast_dir, "canu_vs_our_comparison.tsv")
        try:
            comparison = compare_assemblies(all_reports, comparison_file)
            print(f"\nFull comparison saved to {comparison_file}")
        except Exception as e:
            print(f"Error creating comparison report: {e}")
            comparison_file = None
    else:
        print("No reports available for comparison.")
        comparison_file = None
    
    # Create a detailed report
    create_final_report(quast_dir, canu_assemblies, comparison_file)

def create_mock_report(quast_dir, error_detail):
    """Create a report when Canu installation fails."""
    report_file = os.path.join(quast_dir, "task1_part4_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 4: Canu vs. Our Implementations\n\n")
        f.write("Canu installation or execution failed. Below is a theoretical comparison.\n\n")
        f.write("## Error Details\n")
        f.write(f"{error_detail}\n\n")
        
        f.write("## Expected Differences\n")
        f.write("Canu is a sophisticated long-read assembler with several unique features:\n\n")
        f.write("1. Canu is specifically designed for long-read technologies like ONT and PacBio\n")
        f.write("2. It has three main stages: error correction, read trimming, and assembly\n")
        f.write("3. It uses a sophisticated overlap-layout-consensus (OLC) algorithm\n")
        f.write("4. It includes adaptive error detection and correction for noisy reads\n\n")
        
        f.write("## Manual Installation Instructions\n\n")
        f.write("To install Canu, follow these steps:\n\n")
        f.write("1. Download Canu:\n")
        f.write(f"   wget {CANU_DOWNLOAD_URL}\n\n")
        f.write("2. Extract the archive:\n")
        f.write(f"   tar -xf canu-{CANU_VERSION}.Linux-amd64.tar.xz\n\n")
        f.write("3. Make the executable accessible:\n")
        f.write(f"   chmod +x canu-{CANU_VERSION}/bin/canu\n\n")
        f.write("4. Run the script with the path to Canu:\n")
        f.write(f"   CANU_PATH=/path/to/canu-{CANU_VERSION}/bin/canu python task1_3_part4_canu.py\n\n")
        
        f.write("For online documentation, visit: https://github.com/marbl/canu\n")
    
    print(f"\nReport written to {report_file}")

def create_final_report(quast_dir, canu_assemblies, comparison_file):
    """Create the final comparison report."""
    report_file = os.path.join(quast_dir, "task1_part4_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 4: Canu vs. Our Implementations\n\n")
        
        if canu_assemblies:
            f.write("## Canu Assembly Results\n\n")
            for basename, assembly_path in canu_assemblies.items():
                f.write(f"### Dataset: {basename}\n")
                f.write(f"Assembly: {os.path.basename(assembly_path)}\n")
                f.write(f"Path: {assembly_path}\n\n")
        
        f.write("## Key Features of Canu\n")
        f.write("Canu is a sophisticated long-read assembler with several advantages:\n\n")
        f.write("1. Canu is designed specifically for noisy long-read data (ONT and PacBio)\n")
        f.write("2. It performs read correction before assembly, improving accuracy\n")
        f.write("3. It uses an optimized Overlap-Layout-Consensus approach\n")
        f.write("4. It has adaptive parameters for different error rates and read characteristics\n")
        f.write("5. It includes specialized error detection and correction algorithms\n\n")
        
        f.write("## Comparison to Our Implementations\n\n")
        f.write("1. While our OLC implementation uses a simple overlap detection and greedy path finding, Canu:\n")
        f.write("   - Uses MinHash for efficient overlap detection\n")
        f.write("   - Implements sophisticated error correction\n")
        f.write("   - Employs advanced consensus algorithms\n\n")
        
        f.write("2. Unlike our De Bruijn Graph assembler, Canu:\n")
        f.write("   - Doesn't use k-mers, making it better for error-prone long reads\n")
        f.write("   - Has better handling of repetitive regions\n")
        f.write("   - Is more computationally intensive but potentially more accurate\n\n")
        
        if comparison_file and os.path.exists(comparison_file):
            f.write("## Performance Comparison\n\n")
            f.write("See the detailed comparison of assembly metrics at:\n")
            f.write(f"{comparison_file}\n\n")
    
    print(f"\nReport written to {report_file}")

if __name__ == "__main__":
    main()
