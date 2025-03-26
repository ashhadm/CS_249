#!/bin/bash
# Script to download SRA samples and run Kraken2 analysis
# For CS249 Assignment 1, Task 3.2
#sra.sh
# Set up directories
BASE_DIR="$(pwd)/metagenomic_analysis"
DB_DIR="${BASE_DIR}/kraken2_db"
SAMPLES_DIR="${BASE_DIR}/samples"
RESULTS_DIR="${BASE_DIR}/results"

# Create directories
mkdir -p ${DB_DIR} ${SAMPLES_DIR} ${RESULTS_DIR}

# List of SRA accessions to analyze
declare -a SRA_ACCESSIONS=(
  "SRR11412973" "SRR11412976" "SRR11412979" "SRR11412980" "SRR11412984"
  "SRR21907296" "SRR21907303" "SRR21907307" "SRR21907332" "SRR21907330"
)

# Step 1: Download Kraken2 Standard-8 database from the link provided in the assignment
echo "=== Downloading Kraken2 Standard-8 database ==="
cd ${DB_DIR}

# URL from the assignment (benlangmead.github.io/aws-indexes/k2)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz

# Extract the database
tar -xzf k2_standard_08gb_20230314.tar.gz
rm k2_standard_08gb_20230314.tar.gz
echo "Database downloaded and extracted to ${DB_DIR}"

# Check if database files are present
if [ ! -f "${DB_DIR}/hash.k2d" ] || [ ! -f "${DB_DIR}/opts.k2d" ] || [ ! -f "${DB_DIR}/taxo.k2d" ]; then
    echo "Warning: Some database files may be missing. Continuing anyway..."
fi

# Step 2: Download and process samples
echo "=== Downloading and processing SRA samples ==="

# Check if SRA toolkit is installed
if ! command -v prefetch &> /dev/null || ! command -v fastq-dump &> /dev/null; then
    echo "SRA toolkit not found. Please install it first."
    echo "Visit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit"
    exit 1
fi

# Check if Kraken2 is installed
if ! command -v kraken2 &> /dev/null; then
    echo "Kraken2 not found. Please install it first."
    exit 1
fi

# Download and process each sample
for accession in "${SRA_ACCESSIONS[@]}"; do
    echo "Processing sample ${accession}..."
    
    # Create sample directory
    sample_dir="${SAMPLES_DIR}/${accession}"
    mkdir -p ${sample_dir}
    
    # Download SRA data
    cd ${sample_dir}
    echo "  Downloading ${accession}..."
    prefetch ${accession}
    
    # Convert to FASTQ
    echo "  Converting to FASTQ format..."
    fastq-dump --split-files --gzip ${accession}/${accession}.sra
    
    # Run Kraken2 on the sample
    echo "  Running Kraken2 analysis..."
    
    # Create results directory for this sample
    result_dir="${RESULTS_DIR}/${accession}"
    mkdir -p ${result_dir}
    
    # Check if we have paired-end or single-end reads
    if [ -f "${sample_dir}/${accession}_2.fastq.gz" ]; then
        # Paired-end analysis
        kraken2 --db ${DB_DIR} \
                --paired \
                --output ${result_dir}/kraken2_output.txt \
                --report ${result_dir}/kraken2_report.txt \
                --use-names \
                --threads 8 \
                ${sample_dir}/${accession}_1.fastq.gz \
                ${sample_dir}/${accession}_2.fastq.gz
    else
        # Single-end analysis
        kraken2 --db ${DB_DIR} \
                --output ${result_dir}/kraken2_output.txt \
                --report ${result_dir}/kraken2_report.txt \
                --use-names \
                --threads 8 \
                ${sample_dir}/${accession}_1.fastq.gz
    fi
    
    echo "  Analysis complete for ${accession}"
done

# Create a summary file with sample info
echo "=== Creating summary of samples ==="
summary_file="${RESULTS_DIR}/sample_summary.txt"
echo "Sample ID,Total Reads,Classified Reads,Unclassified Reads,Classification Rate" > ${summary_file}

for accession in "${SRA_ACCESSIONS[@]}"; do
    result_file="${RESULTS_DIR}/${accession}/kraken2_report.txt"
    if [ -f "$result_file" ]; then
        # Extract classification statistics from the report
        total=$(grep -m1 "^100.00" ${result_file} | awk '{print $2}')
        classified=$(grep -m1 "^  classified" ${result_file} | awk '{print $2}')
        unclassified=$(grep -m1 "^  unclassified" ${result_file} | awk '{print $2}')
        rate=$(grep -m1 "^  classified" ${result_file} | awk '{print $3}' | tr -d "()")
        
        echo "${accession},${total},${classified},${unclassified},${rate}" >> ${summary_file}
    else
        echo "${accession},Error: No results found" >> ${summary_file}
    fi
done

echo "=== Analysis complete ==="
echo "Results are stored in ${RESULTS_DIR}"
echo "A summary is available at ${summary_file}"
