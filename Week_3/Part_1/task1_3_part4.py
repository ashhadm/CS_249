#!/usr/bin/env python3
# Run SPAdes on MERS datasets and compare with our implementations

import os
import sys
import subprocess
import pandas as pd
import urllib.request
import tarfile
import platform
import shutil
from pathlib import Path
import socket

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from genome_assembly_src.evaluation import run_quast, parse_quast_report, compare_assemblies

# SPAdes installation details
SPADES_VERSION = "3.15.5"
SPADES_DOWNLOAD_URLS = [
    f"http://cab.spbu.ru/files/release{SPADES_VERSION}/SPAdes-{SPADES_VERSION}-Linux.tar.gz",
    f"https://github.com/ablab/spades/releases/download/v{SPADES_VERSION}/SPAdes-{SPADES_VERSION}-Linux.tar.gz"
]
LOCAL_SPADES_DIR = "spades_install"  # Local directory to install SPAdes

def check_internet_connection():
    """Check if there's an active internet connection by attempting to connect to Google's DNS."""
    try:
        # Attempt to connect to Google's DNS (8.8.8.8)
        socket.create_connection(("8.8.8.8", 53), timeout=3)
        return True
    except OSError:
        return False

def check_spades_installed():
    """Check if SPAdes is installed and accessible."""
    # First check if provided by manual installation path
    manual_path = os.environ.get("SPADES_PATH")
    if manual_path and os.path.exists(manual_path) and os.access(manual_path, os.X_OK):
        print(f"SPAdes found at environment variable SPADES_PATH: {manual_path}")
        return "manual", manual_path
    
    # Check common installation locations
    common_paths = [
        "spades.py",
        "/usr/bin/spades.py",
        "/usr/local/bin/spades.py",
        os.path.expanduser("~/SPAdes-3.15.5-Linux/bin/spades.py"),
        "SPAdes-3.15.5-Linux/bin/spades.py",
        f"{LOCAL_SPADES_DIR}/SPAdes-{SPADES_VERSION}-Linux/bin/spades.py"
    ]
    
    for path in common_paths:
        try:
            result = subprocess.run([path, "--version"], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
            if result.returncode == 0:
                print(f"SPAdes found at: {path}")
                return "found", path
        except (FileNotFoundError, PermissionError):
            continue
    
    # Check if downloaded but not executable
    local_spades = os.path.join(LOCAL_SPADES_DIR, f"SPAdes-{SPADES_VERSION}-Linux/bin/spades.py")
    if os.path.exists(local_spades):
        try:
            os.chmod(local_spades, 0o755)
            print(f"Found SPAdes at {local_spades} and made it executable")
            return "local", local_spades
        except:
            pass
            
    return None, None

def install_spades():
    """Download and install SPAdes locally."""
    if not check_internet_connection():
        print("No internet connection detected. Manual installation required.")
        print_manual_installation_instructions()
        return None
    
    print(f"Downloading SPAdes {SPADES_VERSION}...")
    
    # Create directory for SPAdes
    os.makedirs(LOCAL_SPADES_DIR, exist_ok=True)
    
    # Try different download URLs
    downloaded = False
    error_messages = []
    tarball_path = os.path.join(LOCAL_SPADES_DIR, f"SPAdes-{SPADES_VERSION}-Linux.tar.gz")
    
    for url in SPADES_DOWNLOAD_URLS:
        try:
            print(f"Trying to download from: {url}")
            urllib.request.urlretrieve(url, tarball_path)
            downloaded = True
            break
        except Exception as e:
            error_messages.append(f"Failed to download from {url}: {str(e)}")
            continue
    
    if not downloaded:
        print("Failed to download SPAdes from all sources:")
        for error in error_messages:
            print(f" - {error}")
        print_manual_installation_instructions()
        return None
    
    print(f"Extracting SPAdes...")
    try:
        with tarfile.open(tarball_path) as tar:
            tar.extractall(path=LOCAL_SPADES_DIR)
        
        # Clean up tarball
        os.remove(tarball_path)
        
        # Make spades.py executable
        spades_path = os.path.join(LOCAL_SPADES_DIR, f"SPAdes-{SPADES_VERSION}-Linux/bin/spades.py")
        os.chmod(spades_path, 0o755)
        
        print(f"SPAdes {SPADES_VERSION} installed successfully to {os.path.abspath(LOCAL_SPADES_DIR)}")
        return spades_path
    except Exception as e:
        print(f"Error extracting SPAdes: {e}")
        print_manual_installation_instructions()
        return None

def print_manual_installation_instructions():
    """Print detailed instructions for manual SPAdes installation."""
    print("\n=== Manual SPAdes Installation Instructions ===")
    print("1. Download SPAdes from one of these URLs:")
    for url in SPADES_DOWNLOAD_URLS:
        print(f"   - {url}")
    print("   Alternative source: https://github.com/ablab/spades/releases")
    print("\n2. Extract the tarball:")
    print("   tar -xzf SPAdes-3.15.5-Linux.tar.gz")
    print("\n3. Set executable permissions:")
    print("   chmod +x SPAdes-3.15.5-Linux/bin/spades.py")
    print("\n4. Run this script again with the SPAdes path:")
    print("   SPADES_PATH=/path/to/SPAdes-3.15.5-Linux/bin/spades.py python task1_3_part4.py")
    print("\nNote: SPAdes also requires Python 3.7+ and several dependencies.")
    print("For complete installation instructions, see: http://cab.spbu.ru/software/spades/")
    print("================================================\n")

def run_spades(spades_path, input_fastq, output_dir):
    """Run SPAdes on a single fastq file."""
    os.makedirs(output_dir, exist_ok=True)
    
    # SPAdes parameters: using --s 1 for single-end reads
    cmd = [
        spades_path,
        "-s", input_fastq,
        "-o", output_dir,
        "--isolate",  # Optimized for isolate genomes (like viral genomes)
        "--phred-offset", "33",  # Standard FASTQ quality encoding
        "-t", "4"  # Use 4 threads
    ]
    
    print(f"Running SPAdes with command: {' '.join(cmd)}")
    
    try:
        # Run SPAdes with live output
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
            print(f"SPAdes failed with return code {process.returncode}")
            return None
        
        # Check if assembly exists
        assembly_path = os.path.join(output_dir, "contigs.fasta")
        if not os.path.exists(assembly_path):
            # Try scaffolds instead
            assembly_path = os.path.join(output_dir, "scaffolds.fasta")
            if not os.path.exists(assembly_path):
                print("SPAdes did not produce any assembly output")
                return None
        
        return assembly_path
    
    except Exception as e:
        print(f"Error running SPAdes: {e}")
        return None

def main():
    # Set up paths
    reference_file = "cs249_hw2/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna"
    output_dir = "genome_assembly_output/spades_results"
    quast_dir = "genome_assembly_output/quast_results"
    
    # Create output directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(quast_dir, exist_ok=True)
    
    # Dataset paths (only use Illumina-like reads for SPAdes)
    datasets = [
        "cs249_hw2/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq",
        "cs249_hw2/synthetic_dataset/reads/reads_hiseq_5k.fastq"
    ]
    
    # Check if datasets exist
    for dataset in datasets:
        if not os.path.exists(dataset):
            print(f"Error: Dataset file {dataset} not found.")
            return
    
    print("\n========== Running Task 1.3 Part 4 ==========\n")
    
    # Check if SPAdes is installed, install if not
    install_type, spades_path = check_spades_installed()
    
    if not spades_path:
        print("SPAdes not found. Attempting to install...")
        try:
            spades_path = install_spades()
            install_type = "local"
        except Exception as e:
            print(f"Error installing SPAdes: {e}")
            create_mock_report(quast_dir, e)
            return
    
    if not spades_path:
        print("Failed to find or install SPAdes.")
        create_mock_report(quast_dir, "Installation failed. See above for details.")
        return
    
    # Process each dataset with SPAdes
    reports = []
    spades_assemblies = {}
    
    for dataset in datasets:
        basename = os.path.basename(dataset).split('.')[0]
        print(f"\nProcessing dataset: {basename}")
        
        # Set up SPAdes output directory
        spades_out_dir = os.path.join(output_dir, basename)
        
        # Run SPAdes
        print(f"Running SPAdes on {basename}...")
        assembly_path = run_spades(spades_path, dataset, spades_out_dir)
        
        if assembly_path:
            print(f"SPAdes completed successfully for {basename}")
            spades_assemblies[basename] = assembly_path
            
            # Run QUAST
            print("Evaluating SPAdes assembly with QUAST...")
            quast_spades_dir = os.path.join(quast_dir, f"{basename}_spades")
            quast_spades_report = run_quast(assembly_path, reference_file, quast_spades_dir)
            
            # Add to reports list
            if quast_spades_report and os.path.exists(quast_spades_report):
                reports.append((f"{basename}_spades", quast_spades_report))
            else:
                print(f"Warning: QUAST evaluation failed for {basename}")
        else:
            print(f"SPAdes assembly failed for {basename}")
    
    # If no SPAdes assemblies were successful, exit
    if not spades_assemblies:
        print("\nNo successful SPAdes assemblies were generated.")
        create_mock_report(quast_dir, "SPAdes ran but failed to produce any assemblies.")
        return
    
    # Compare SPAdes assemblies with our implementations
    print("\nComparing SPAdes assemblies with our implementations...")
    
    # Load our implementation reports
    our_reports = []
    for dataset in datasets:
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
        comparison_file = os.path.join(quast_dir, "spades_vs_our_comparison.tsv")
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
    create_final_report(quast_dir, spades_assemblies, comparison_file)

def create_mock_report(quast_dir, error_detail):
    """Create a report when SPAdes installation fails."""
    report_file = os.path.join(quast_dir, "task1_part4_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 4: SPAdes vs. Our Implementations\n\n")
        f.write("SPAdes installation or execution failed. Below is a theoretical comparison.\n\n")
        f.write("## Error Details\n")
        f.write(f"{error_detail}\n\n")
        
        f.write("## Expected Differences\n")
        f.write("SPAdes is a sophisticated assembler that would likely outperform our implementations:\n\n")
        f.write("1. SPAdes uses multiple k-mer sizes to build a combined de Bruijn graph\n")
        f.write("2. It incorporates error correction using BayesHammer\n")
        f.write("3. It uses sophisticated algorithms for bubble removal and repeat resolution\n")
        f.write("4. It includes advanced graph simplification techniques\n\n")
        
        f.write("## Manual Installation Instructions\n\n")
        f.write("To install SPAdes, follow these steps:\n\n")
        f.write("1. Download SPAdes:\n")
        f.write("   wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz\n\n")
        f.write("2. Extract the archive:\n")
        f.write("   tar -xzf SPAdes-3.15.5-Linux.tar.gz\n\n")
        f.write("3. Make the executable accessible:\n")
        f.write("   chmod +x SPAdes-3.15.5-Linux/bin/spades.py\n\n")
        f.write("4. Run the script with the path to SPAdes:\n")
        f.write("   SPADES_PATH=/path/to/SPAdes-3.15.5-Linux/bin/spades.py python task1_3_part4.py\n\n")
        
        f.write("For online documentation, visit: http://cab.spbu.ru/software/spades/\n")
    
    print(f"\nReport written to {report_file}")

def create_final_report(quast_dir, spades_assemblies, comparison_file):
    """Create the final comparison report."""
    report_file = os.path.join(quast_dir, "task1_part4_report.txt")
    with open(report_file, 'w') as f:
        f.write("# Task 1.3 Part 4: SPAdes vs. Our Implementations\n\n")
        
        if spades_assemblies:
            f.write("## SPAdes Assembly Results\n\n")
            for basename, assembly_path in spades_assemblies.items():
                f.write(f"### Dataset: {basename}\n")
                f.write(f"Assembly: {os.path.basename(assembly_path)}\n")
                f.write(f"Path: {assembly_path}\n\n")
        
        f.write("## Key Features of SPAdes\n")
        f.write("SPAdes is a sophisticated assembler with several advantages over our implementations:\n\n")
        f.write("1. SPAdes uses multiple k-mer sizes to build a combined de Bruijn graph\n")
        f.write("2. It incorporates error correction using BayesHammer\n")
        f.write("3. It uses sophisticated algorithms for bubble removal and repeat resolution\n")
        f.write("4. It includes advanced graph simplification techniques\n")
        f.write("5. It can utilize paired-end and mate-pair information\n\n")
        
        if comparison_file and os.path.exists(comparison_file):
            f.write("## Performance Comparison\n\n")
            f.write("See the detailed comparison of assembly metrics at:\n")
            f.write(f"{comparison_file}\n\n")
            
            # Try to highlight key differences
            f.write("### Key Metrics Comparison\n\n")
            try:
                comparison_data = pd.read_csv(comparison_file, sep='\t')
                # Extract a simplified comparison here
                f.write("Key metrics comparison available in the TSV file.\n")
            except:
                f.write("Metrics comparison unavailable. Please check the TSV file directly.\n")
    
    print(f"\nReport written to {report_file}")

if __name__ == "__main__":
    main()
