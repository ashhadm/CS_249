# genome_assembly_src/evaluation.py
import os
import subprocess
import pandas as pd

def run_quast(assembly_file, reference_file=None, output_dir="./quast_results"):
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct QUAST command with lower min-contig
    cmd = ["quast.py", assembly_file, "-o", output_dir, "--min-contig", "100"]
    if reference_file:
        cmd.extend(["-r", reference_file])
    
    # Run QUAST
    print(f"Running QUAST: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return os.path.join(output_dir, "report.tsv")

def parse_quast_report(report_file):
    """Parse QUAST report and return metrics as a dictionary."""
    try:
        df = pd.read_csv(report_file, sep='\t')
        metrics = {}
        
        # Convert DataFrame to dictionary
        for _, row in df.iterrows():
            if not pd.isna(row.iloc[0]):
                metrics[row.iloc[0]] = row.iloc[1]
            
        return metrics
    except Exception as e:
        print(f"Error parsing QUAST report: {e}")
        return {}

def compare_assemblies(reports, output_file="assembly_comparison.tsv"):
    """
    Compare multiple QUAST reports and save results to a file.
    
    Args:
        reports: List of tuples (assembly_name, report_path)
        output_file: Path to output comparison file
        
    Returns:
        DataFrame with comparison results
    """
    all_metrics = {}
    
    for name, report in reports:
        metrics = parse_quast_report(report)
        all_metrics[name] = metrics
    
    # Convert to DataFrame
    df = pd.DataFrame(all_metrics)
    
    # Save to file
    df.to_csv(output_file, sep='\t')
    
    return df
