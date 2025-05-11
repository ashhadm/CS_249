#!/bin/bash
#SBATCH --account=cs249
#SBATCH --job-name=busco_misasm
#SBATCH --output=busco_misasm_%j.out
#SBATCH --error=busco_misasm_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Print start time
echo "====== Script started at: $(date) ======"

# Set working directory
WORK_DIR="$PWD/lizard_assembly_hifiasm"
cd $WORK_DIR

# Define paths
ASSEMBLY_FILE="$WORK_DIR/final_assembly.fasta"
EVAL_DIR="$WORK_DIR/evaluation"
BUSCO_DIR="$EVAL_DIR/busco_new"
MISASM_DIR="$EVAL_DIR/misasm_new"

mkdir -p $BUSCO_DIR $MISASM_DIR

# Check if assembly file exists
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "ERROR: Assembly file not found: $ASSEMBLY_FILE"
    echo "Please run the conversion script first."
    exit 1
fi

# Define paths for raw reads
HIFI_READS="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"

# Create a log file for the evaluation
PARALLEL_LOG="$EVAL_DIR/busco_misasm.log"
touch $PARALLEL_LOG

# Helper function to log messages
log_message() {
    local message="$1"
    echo "$message" | tee -a $PARALLEL_LOG
}

log_message "====== BUSCO and Mis-assembly Detection Started: $(date) ======"
log_message "Assembly file: $ASSEMBLY_FILE"

###########################################
# Run BUSCO with explicit lineage setting
###########################################
log_message "====== Running BUSCO with explicit lineage at: $(date) ======"

module purge
module load busco

if [ $? -eq 0 ]; then
    log_message "BUSCO module loaded successfully."
    
    # List available lineages
    log_message "Available BUSCO lineages:"
    busco --list-datasets | tee -a $PARALLEL_LOG
    
    # Set explicit lineage - using vertebrata_odb10 which should be available
    LINEAGE="vertebrata_odb10"
    log_message "Using lineage: $LINEAGE"
    
    # Run BUSCO with explicit parameters
    log_message "Running BUSCO on primary assembly..."
    cd $BUSCO_DIR
    busco -i $ASSEMBLY_FILE -o primary -m genome -l $LINEAGE --download_path $BUSCO_DIR/downloads --cpu $SLURM_CPUS_PER_TASK --force
    
    if [ $? -eq 0 ]; then
        log_message "BUSCO completed successfully on primary assembly."
        log_message "BUSCO summary for primary assembly:"
        cat $BUSCO_DIR/primary/short_summary*.txt | tee -a $PARALLEL_LOG
    else
        log_message "ERROR: BUSCO failed on primary assembly."
        
        # Try with a different lineage
        log_message "Trying with tetrapoda_odb10 lineage..."
        busco -i $ASSEMBLY_FILE -o primary_tetrapoda -m genome -l tetrapoda_odb10 --download_path $BUSCO_DIR/downloads --cpu $SLURM_CPUS_PER_TASK --force
        
        if [ $? -eq 0 ]; then
            log_message "BUSCO with tetrapoda_odb10 completed successfully."
            log_message "BUSCO summary for primary assembly (tetrapoda_odb10):"
            cat $BUSCO_DIR/primary_tetrapoda/short_summary*.txt | tee -a $PARALLEL_LOG
        else
            log_message "ERROR: BUSCO with tetrapoda_odb10 also failed."
        fi
    fi
    
    # Return to working directory
    cd $WORK_DIR
else
    log_message "ERROR: Failed to load BUSCO module. Skipping BUSCO evaluation."
fi

###########################################
# Run mis-assembly detection
###########################################
log_message "====== Running mis-assembly detection at: $(date) ======"

module purge
module load flagger 2>/dev/null || module load inspector 2>/dev/null || module load minimap2

# Check which tool is available
if command -v flagger &>/dev/null; then
    log_message "Using Flagger for mis-assembly detection..."
    
    # Run Flagger on primary assembly
    cd $MISASM_DIR
    flagger $ASSEMBLY_FILE $HIFI_READS -o $MISASM_DIR/flagger_output -t $SLURM_CPUS_PER_TASK
    
    if [ $? -eq 0 ]; then
        log_message "Flagger completed successfully."
        log_message "Flagger results:"
        ls -lh $MISASM_DIR/flagger_output | tee -a $PARALLEL_LOG
    else
        log_message "ERROR: Flagger failed."
    fi
    
    # Return to working directory
    cd $WORK_DIR
elif command -v inspector.py &>/dev/null; then
    log_message "Using Inspector for mis-assembly detection..."
    
    # Run Inspector on primary assembly
    cd $MISASM_DIR
    inspector.py -c $ASSEMBLY_FILE -r $HIFI_READS -o $MISASM_DIR/inspector_output -t $SLURM_CPUS_PER_TASK
    
    if [ $? -eq 0 ]; then
        log_message "Inspector completed successfully."
        log_message "Inspector results:"
        ls -lh $MISASM_DIR/inspector_output | tee -a $PARALLEL_LOG
    else
        log_message "ERROR: Inspector failed."
    fi
    
    # Return to working directory
    cd $WORK_DIR
elif command -v minimap2 &>/dev/null; then
    log_message "Using minimap2 for basic mis-assembly detection..."
    
    # Run minimap2 for read-to-assembly alignment
    log_message "Aligning HiFi reads to assembly with minimap2..."
    
    # Create output directory
    mkdir -p $MISASM_DIR/minimap2_output
    
    # Run alignment
    minimap2 -ax map-hifi $ASSEMBLY_FILE $HIFI_READS -t $SLURM_CPUS_PER_TASK > $MISASM_DIR/minimap2_output/hifi_to_assembly.sam
    
    if [ $? -eq 0 ]; then
        log_message "Minimap2 alignment completed successfully."
        log_message "SAM file: $MISASM_DIR/minimap2_output/hifi_to_assembly.sam"
        
        # Convert SAM to BAM and sort for easier analysis
        if command -v samtools &>/dev/null; then
            log_message "Converting SAM to sorted BAM..."
            samtools view -bS $MISASM_DIR/minimap2_output/hifi_to_assembly.sam | samtools sort -o $MISASM_DIR/minimap2_output/hifi_to_assembly.sorted.bam -
            
            if [ $? -eq 0 ]; then
                log_message "BAM conversion and sorting completed successfully."
                
                # Calculate basic mapping statistics
                log_message "Calculating mapping statistics..."
                samtools flagstat $MISASM_DIR/minimap2_output/hifi_to_assembly.sorted.bam > $MISASM_DIR/minimap2_output/mapping_stats.txt
                
                log_message "Mapping statistics:"
                cat $MISASM_DIR/minimap2_output/mapping_stats.txt | tee -a $PARALLEL_LOG
                
                # Index the BAM file
                samtools index $MISASM_DIR/minimap2_output/hifi_to_assembly.sorted.bam
                
                # Identify potential misassembly regions (low coverage or many clipped reads)
                log_message "Identifying potential misassembly regions..."
                
                # Look for regions with many clipped reads
                samtools view $MISASM_DIR/minimap2_output/hifi_to_assembly.sorted.bam | \
                  awk '{if($6~/S/ && length($6)>=4) print $3"\t"$4"\t"$6}' | \
                  sort -k1,1 -k2,2n > $MISASM_DIR/minimap2_output/clipped_regions.txt
                
                log_message "Clipped regions saved to: $MISASM_DIR/minimap2_output/clipped_regions.txt"
            fi
        fi
    else
        log_message "ERROR: Minimap2 alignment failed."
    fi
else
    log_message "ERROR: No tools available for mis-assembly detection. Skipping this step."
fi

log_message "====== BUSCO and Mis-assembly Detection Completed: $(date) ======"
echo "====== Script completed at: $(date) ======"
