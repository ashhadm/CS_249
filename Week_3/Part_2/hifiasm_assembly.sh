#!/bin/bash
#SBATCH --account=cs249
#SBATCH --job-name=hifiasm_basic
#SBATCH --output=hifiasm_basic_%j.out
#SBATCH --error=hifiasm_basic_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G

# Print start time
echo "====== Script started at: $(date) ======"

# Create working directory
WORK_DIR="$PWD/lizard_assembly_hifiasm"
mkdir -p $WORK_DIR
cd $WORK_DIR

# Define input and output paths
HIFI_READS="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"
OUTPUT_PREFIX="$WORK_DIR/lizard_hifiasm"

# Log file paths
echo "HIFI_READS: $HIFI_READS"
echo "OUTPUT_PREFIX: $OUTPUT_PREFIX"

# Load necessary modules
echo "Loading modules..."
module purge
module load hifiasm

# Check if Hifiasm is available
if [ $? -ne 0 ]; then
    echo "ERROR: Hifiasm module not available. Checking alternative modules..."
    module avail 2>&1 | grep -i hifiasm
    
    # Try to find hifiasm in other modules or bioconda
    module load bio/conda-apps 2>/dev/null || module load bioconda/cf201901/python2.7.18 2>/dev/null
    
    if ! command -v hifiasm &>/dev/null; then
        echo "ERROR: Failed to load Hifiasm. Exiting."
        exit 1
    fi
fi

# Print Hifiasm version
hifiasm --version || echo "Hifiasm version command unavailable"

echo "Starting Hifiasm basic assembly at: $(date)"

# Run HiFi-only assembly with maximum threads
hifiasm -o ${OUTPUT_PREFIX} -t $SLURM_CPUS_PER_TASK $HIFI_READS

# Check if assembly completed successfully
if [ $? -eq 0 ]; then
    echo "Hifiasm assembly completed successfully at: $(date)"
    
    # List output files
    echo "Assembly output files:"
    ls -lh ${OUTPUT_PREFIX}*
    
    # Convert the primary contig GFA file to FASTA format
    echo "Converting primary contigs GFA to FASTA format..."
    
    if [ -f "${OUTPUT_PREFIX}.p_ctg.gfa" ]; then
        awk '/^S/{print ">"$2;print $3}' ${OUTPUT_PREFIX}.p_ctg.gfa > ${OUTPUT_PREFIX}.p_ctg.fasta
        echo "Created FASTA file: ${OUTPUT_PREFIX}.p_ctg.fasta"
        
        # Create symlink to the final assembly for easy access
        ln -sf ${OUTPUT_PREFIX}.p_ctg.fasta $WORK_DIR/final_assembly.fasta
        echo "Created symlink: $WORK_DIR/final_assembly.fasta"
        
        # Get basic file stats
        echo "Assembly file size: $(ls -lh ${OUTPUT_PREFIX}.p_ctg.fasta | awk '{print $5}')"
        echo "Number of contigs: $(grep -c "^>" ${OUTPUT_PREFIX}.p_ctg.fasta)"
    else
        echo "WARNING: Primary contig GFA file not found: ${OUTPUT_PREFIX}.p_ctg.gfa"
    fi
    
    # Also convert the alternative (haplotype-resolved) contigs if available
    if [ -f "${OUTPUT_PREFIX}.a_ctg.gfa" ]; then
        awk '/^S/{print ">"$2;print $3}' ${OUTPUT_PREFIX}.a_ctg.gfa > ${OUTPUT_PREFIX}.a_ctg.fasta
        echo "Created alternative haplotype FASTA file: ${OUTPUT_PREFIX}.a_ctg.fasta"
    fi
else
    echo "ERROR: Hifiasm assembly failed at: $(date)"
fi

echo "====== Script completed at: $(date) ======"
