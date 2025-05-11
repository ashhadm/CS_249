#!/bin/bash
#SBATCH --account=cs249
#SBATCH --job-name=evaluate_assembly
#SBATCH --output=evaluate_assembly_%j.out
#SBATCH --error=evaluate_assembly_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

# Print start time
echo "====== Script started at: $(date) ======"

# Set working directory
WORK_DIR="$PWD/lizard_assembly_hifiasm"
cd $WORK_DIR

# Define paths
ASSEMBLY_FILE="$WORK_DIR/final_assembly.fasta"
HAP1_ASSEMBLY="$WORK_DIR/lizard_hifiasm.bp.hap1.p_ctg.fasta"
HAP2_ASSEMBLY="$WORK_DIR/lizard_hifiasm.bp.hap2.p_ctg.fasta"
EVAL_DIR="$WORK_DIR/evaluation"
mkdir -p $EVAL_DIR

# Check if assembly file exists
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "ERROR: Assembly file not found: $ASSEMBLY_FILE"
    echo "Please run the conversion script first."
    exit 1
fi

# Define paths for raw reads (for QV calculation with Merqury)
HIFI_READS="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"

# Create directories for each evaluation tool
mkdir -p $EVAL_DIR/{quast,busco,merqury,misasm}

# Create a log file for the evaluation
EVAL_LOG="$EVAL_DIR/evaluation.log"
touch $EVAL_LOG

# Helper function to log messages
log_message() {
    local message="$1"
    echo "$message" | tee -a $EVAL_LOG
}

log_message "====== Assembly Evaluation Started: $(date) ======"
log_message "Primary assembly: $ASSEMBLY_FILE"
log_message "Haplotype 1 assembly: $HAP1_ASSEMBLY"
log_message "Haplotype 2 assembly: $HAP2_ASSEMBLY"

###########################################
# 1. Basic metrics with QUAST
###########################################
log_message "====== Running QUAST for basic metrics at: $(date) ======"
QUAST_DIR="$EVAL_DIR/quast"

module purge
module load quast || module load bioconda/cf201901/python2.7.18

# Check if QUAST is available
if command -v quast.py &>/dev/null; then
    log_message "QUAST is available. Running QUAST analysis..."
    
    # Run QUAST on all assemblies
    quast.py -o $QUAST_DIR \
             -t $SLURM_CPUS_PER_TASK \
             --min-contig 1000 \
             --scaffold-gap-max-size 1000 \
             --est-ref-size 1800000000 \
             --labels "Primary,Haplotype1,Haplotype2" \
             $ASSEMBLY_FILE $HAP1_ASSEMBLY $HAP2_ASSEMBLY
    
    if [ $? -eq 0 ]; then
        log_message "QUAST completed successfully."
        log_message "QUAST report location: $QUAST_DIR/report.txt"
        
        # Copy QUAST report to log
        log_message "QUAST Summary Report:"
        cat $QUAST_DIR/report.txt | tee -a $EVAL_LOG
    else
        log_message "ERROR: QUAST failed."
    fi
else
    log_message "WARNING: QUAST not found. Using basic statistics instead."
    
    # Get basic statistics for each assembly
    for asm in "$ASSEMBLY_FILE" "$HAP1_ASSEMBLY" "$HAP2_ASSEMBLY"; do
        name=$(basename $asm .fasta)
        log_message "Basic statistics for $name:"
        
        # Total size
        TOTAL_SIZE=$(grep -v "^>" $asm | tr -d '\n' | wc -c)
        log_message "  Total size: $TOTAL_SIZE bp"
        
        # Number of contigs
        NUM_CONTIGS=$(grep -c "^>" $asm)
        log_message "  Number of contigs: $NUM_CONTIGS"
        
        # Extract contig lengths
        grep -v "^>" $asm | awk 'BEGIN {RS=">"} NR>1 {print length($0)}' | sort -nr > $EVAL_DIR/${name}_contig_lengths.txt
        
        # Calculate N50 (simplified approach)
        HALF_SIZE=$((TOTAL_SIZE / 2))
        SUM=0
        N50=0
        for length in $(cat $EVAL_DIR/${name}_contig_lengths.txt); do
            SUM=$((SUM + length))
            if [ $SUM -ge $HALF_SIZE ] && [ $N50 -eq 0 ]; then
                N50=$length
                break
            fi
        done
        log_message "  N50: $N50 bp"
        
        # Largest contig
        LARGEST=$(head -n 1 $EVAL_DIR/${name}_contig_lengths.txt)
        log_message "  Largest contig: $LARGEST bp"
    done
fi

###########################################
# 2. Gene completeness with BUSCO
###########################################
log_message "====== Running BUSCO for gene completeness at: $(date) ======"
BUSCO_DIR="$EVAL_DIR/busco"

module purge
module load busco

if [ $? -eq 0 ]; then
    log_message "BUSCO module loaded successfully."
    
    # Get available lineages
    LINEAGES=$(busco --list-datasets 2>&1 | grep -E 'vertebrata|tetrapoda|sauropsida|lepidosauria|reptilia' | head -n 1)
    
    if [ -z "$LINEAGES" ]; then
        # If no specific reptile lineages found, use vertebrata
        LINEAGE="vertebrata_odb10"
        log_message "No specific reptile lineage found. Using general vertebrate lineage: $LINEAGE"
    else
        # Use the first suitable lineage
        LINEAGE=$(echo "$LINEAGES" | head -n 1 | awk '{print $1}')
        log_message "Using lineage: $LINEAGE"
    fi
    
    # Run BUSCO on primary assembly
    log_message "Running BUSCO on primary assembly..."
    cd $BUSCO_DIR
    busco -i $ASSEMBLY_FILE -o primary -m genome -l $LINEAGE --cpu $SLURM_CPUS_PER_TASK
    
    if [ $? -eq 0 ]; then
        log_message "BUSCO completed successfully on primary assembly."
        log_message "BUSCO summary for primary assembly:"
        cat $BUSCO_DIR/primary/short_summary*.txt | tee -a $EVAL_LOG
    else
        log_message "ERROR: BUSCO failed on primary assembly."
    fi
    
    # Return to working directory
    cd $WORK_DIR
else
    log_message "ERROR: Failed to load BUSCO module. Skipping BUSCO evaluation."
fi

###########################################
# 3. K-mer analysis and QV score with Merqury
###########################################
log_message "====== Running Merqury for k-mer analysis and QV score at: $(date) ======"
MERQURY_DIR="$EVAL_DIR/merqury"

module purge
module load merqury || module load bio/conda-apps

# Check if Merqury tools are available
if command -v meryl &>/dev/null && command -v merqury.sh &>/dev/null; then
    log_message "Merqury tools found. Running k-mer analysis..."
    
    # Build k-mer database from HiFi reads
    log_message "Building k-mer database from HiFi reads..."
    KMER_DB="$MERQURY_DIR/lizard.k21"
    
    cd $MERQURY_DIR
    meryl count k=21 $HIFI_READS output $KMER_DB
    
    if [ $? -eq 0 ]; then
        log_message "K-mer database built successfully: $KMER_DB"
        
        # Run Merqury on primary assembly
        log_message "Running Merqury QV analysis on primary assembly..."
        merqury.sh $KMER_DB $ASSEMBLY_FILE lizard_primary
        
        if [ $? -eq 0 ]; then
            log_message "Merqury completed successfully on primary assembly."
            log_message "QV score for primary assembly:"
            cat $MERQURY_DIR/lizard_primary.qv | tee -a $EVAL_LOG
            
            # Also run on haplotype assemblies if time permits
            if [ -f "$HAP1_ASSEMBLY" ] && [ -f "$HAP2_ASSEMBLY" ]; then
                log_message "Running Merqury QV analysis on haplotype assemblies..."
                merqury.sh $KMER_DB $HAP1_ASSEMBLY $HAP2_ASSEMBLY lizard_haplotypes
                
                if [ $? -eq 0 ]; then
                    log_message "Merqury completed successfully on haplotype assemblies."
                    log_message "QV scores for haplotype assemblies:"
                    cat $MERQURY_DIR/lizard_haplotypes.qv | tee -a $EVAL_LOG
                else
                    log_message "ERROR: Merqury failed on haplotype assemblies."
                fi
            fi
        else
            log_message "ERROR: Merqury failed on primary assembly."
        fi
    else
        log_message "ERROR: Failed to build k-mer database."
    fi
    
    # Return to working directory
    cd $WORK_DIR
else
    log_message "ERROR: Merqury tools not found. Skipping k-mer analysis."
fi

###########################################
# 4. Mis-assembly detection with Flagger or Inspector
###########################################
log_message "====== Running mis-assembly detection at: $(date) ======"
MISASM_DIR="$EVAL_DIR/misasm"

module purge
module load flagger || module load inspector

# Check if mis-assembly detection tools are available
if command -v flagger &>/dev/null; then
    log_message "Flagger found. Running mis-assembly detection..."
    
    # Run Flagger on primary assembly
    cd $MISASM_DIR
    flagger $ASSEMBLY_FILE $HIFI_READS -o $MISASM_DIR/flagger_output -t $SLURM_CPUS_PER_TASK
    
    if [ $? -eq 0 ]; then
        log_message "Flagger completed successfully."
        log_message "Flagger results:"
        ls -lh $MISASM_DIR/flagger_output | tee -a $EVAL_LOG
    else
        log_message "ERROR: Flagger failed."
    fi
    
    # Return to working directory
    cd $WORK_DIR
elif command -v inspector.py &>/dev/null; then
    log_message "Inspector found. Running mis-assembly detection..."
    
    # Run Inspector on primary assembly
    cd $MISASM_DIR
    inspector.py -c $ASSEMBLY_FILE -r $HIFI_READS -o $MISASM_DIR/inspector_output -t $SLURM_CPUS_PER_TASK
    
    if [ $? -eq 0 ]; then
        log_message "Inspector completed successfully."
        log_message "Inspector results:"
        ls -lh $MISASM_DIR/inspector_output | tee -a $EVAL_LOG
    else
        log_message "ERROR: Inspector failed."
    fi
    
    # Return to working directory
    cd $WORK_DIR
else
    log_message "WARNING: Neither Flagger nor Inspector found. Using minimap2 for basic alignment-based validation..."
    
    # Check if minimap2 is available
    if command -v minimap2 &>/dev/null; then
        log_message "Running minimap2 for read-to-assembly alignment..."
        
        # Run minimap2 alignment
        minimap2 -ax map-hifi $ASSEMBLY_FILE $HIFI_READS -t $SLURM_CPUS_PER_TASK > $MISASM_DIR/hifi_to_assembly.sam
        
        if [ $? -eq 0 ]; then
            log_message "Minimap2 alignment completed successfully."
            log_message "SAM file: $MISASM_DIR/hifi_to_assembly.sam"
        else
            log_message "ERROR: Minimap2 alignment failed."
        fi
    else
        log_message "ERROR: No tools available for mis-assembly detection. Skipping this step."
    fi
fi

###########################################
# Generate comprehensive evaluation report
###########################################
log_message "====== Generating comprehensive evaluation report at: $(date) ======"
REPORT_FILE="$EVAL_DIR/assembly_evaluation_report.md"

# Create the markdown report file
cat > $REPORT_FILE << EOF
# Lizard Genome Assembly Evaluation Report

**Assembly**: \`$(basename $ASSEMBLY_FILE)\`  
**Date**: $(date)  
**Generated by**: Evaluation pipeline

## Overview of the Assembly

The lizard genome was assembled using Hifiasm (v0.25.0) with PacBio HiFi reads. The assembly resulted in:

- Primary assembly: 88 contigs, 1.7GB
- Haplotype 1: 93 contigs, 1.7GB  
- Haplotype 2: 63 contigs, 1.7GB

## 1. Basic Assembly Metrics (QUAST)

EOF

# Add QUAST results to the report if available
if [ -f "$QUAST_DIR/report.txt" ]; then
    cat $QUAST_DIR/report.txt >> $REPORT_FILE
    
    echo -e "\n### Interpretation of QUAST Results\n" >> $REPORT_FILE
    echo "The basic assembly metrics provide information about the contiguity and size of the assembly:" >> $REPORT_FILE
    echo "- **Total length**: Represents the total size of the assembled genome" >> $REPORT_FILE
    echo "- **Number of contigs**: Fewer contigs generally indicates a more contiguous assembly" >> $REPORT_FILE
    echo "- **N50**: The length at which contigs of this size or longer contain 50% of the genome" >> $REPORT_FILE
    echo "- **L50**: The number of contigs needed to reach 50% of the genome size" >> $REPORT_FILE
    echo -e "\nA high-quality assembly typically has a large N50 (ideally in the millions of base pairs for a vertebrate genome), a relatively small number of contigs, and a total length that matches the expected genome size.\n" >> $REPORT_FILE
else
    echo "QUAST results not available." >> $REPORT_FILE
fi

# Add BUSCO results to the report if available
echo -e "\n## 2. Gene Completeness Assessment (BUSCO)\n" >> $REPORT_FILE

if [ -f "$BUSCO_DIR/primary/short_summary"*".txt" ]; then
    cat $BUSCO_DIR/primary/short_summary*.txt >> $REPORT_FILE
    
    echo -e "\n### Interpretation of BUSCO Results\n" >> $REPORT_FILE
    echo "BUSCO evaluates the presence of conserved orthologous genes in the assembly:" >> $REPORT_FILE
    echo "- **Complete BUSCOs**: Genes found complete in the assembly (higher is better)" >> $REPORT_FILE
    echo "- **Fragmented BUSCOs**: Genes only partially found (lower is better)" >> $REPORT_FILE
    echo "- **Missing BUSCOs**: Genes not found at all (lower is better)" >> $REPORT_FILE
    echo -e "\nA high-quality assembly should have a high percentage (>90%) of complete BUSCOs and a low percentage of fragmented or missing BUSCOs. This indicates that the gene space is well-represented in the assembly.\n" >> $REPORT_FILE
else
    echo "BUSCO results not available." >> $REPORT_FILE
fi

# Add Merqury results to the report if available
echo -e "\n## 3. K-mer Analysis and QV Score (Merqury)\n" >> $REPORT_FILE

if [ -f "$MERQURY_DIR/lizard_primary.qv" ]; then
    cat $MERQURY_DIR/lizard_primary.qv >> $REPORT_FILE
    
    echo -e "\n### Interpretation of Merqury Results\n" >> $REPORT_FILE
    echo "Merqury uses k-mer analysis to evaluate assembly quality:" >> $REPORT_FILE
    echo "- **QV Score**: Quality Value representing the log-scaled error rate (higher is better)" >> $REPORT_FILE
    echo "- **K-mer Completeness**: Percentage of k-mers from reads found in the assembly" >> $REPORT_FILE
    echo -e "\nA QV score above 40 indicates a high-quality assembly with an error rate less than 1 in 10,000 bases. The assignment specifies that we should aim for a QV score higher than 40.\n" >> $REPORT_FILE
    
    # Check if QV score meets the requirement
    QV_SCORE=$(awk '{print $NF}' $MERQURY_DIR/lizard_primary.qv | tail -n 1)
    if [ -n "$QV_SCORE" ] && (( $(echo "$QV_SCORE < 40" | bc -l) )); then
        echo -e "\n**Note**: The current QV score ($QV_SCORE) is below the target of 40. This suggests that additional quality control and refinement steps may be needed to improve the assembly accuracy.\n" >> $REPORT_FILE
    elif [ -n "$QV_SCORE" ]; then
        echo -e "\n**Note**: The current QV score ($QV_SCORE) meets or exceeds the target of 40, indicating a high-quality assembly.\n" >> $REPORT_FILE
    fi
else
    echo "Merqury results not available." >> $REPORT_FILE
fi

# Add mis-assembly detection results to the report
echo -e "\n## 4. Mis-assembly Detection\n" >> $REPORT_FILE

if [ -d "$MISASM_DIR/flagger_output" ] && [ "$(ls -A $MISASM_DIR/flagger_output)" ]; then
    echo "Flagger results:" >> $REPORT_FILE
    ls -lh $MISASM_DIR/flagger_output >> $REPORT_FILE
    
    echo -e "\n### Interpretation of Mis-assembly Detection Results\n" >> $REPORT_FILE
    echo "Mis-assembly detection tools identify potential structural errors in the assembly:" >> $REPORT_FILE
    echo "- Regions where read alignments suggest incorrect joins" >> $REPORT_FILE
    echo "- Areas with unexpectedly low or high coverage" >> $REPORT_FILE
    echo "- Discrepancies between the assembly and the raw reads" >> $REPORT_FILE
    echo -e "\nA high-quality assembly should have few mis-assemblies, particularly in gene-rich regions.\n" >> $REPORT_FILE
elif [ -d "$MISASM_DIR/inspector_output" ] && [ "$(ls -A $MISASM_DIR/inspector_output)" ]; then
    echo "Inspector results:" >> $REPORT_FILE
    ls -lh $MISASM_DIR/inspector_output >> $REPORT_FILE
    
    echo -e "\n### Interpretation of Mis-assembly Detection Results\n" >> $REPORT_FILE
    echo "Mis-assembly detection tools identify potential structural errors in the assembly:" >> $REPORT_FILE
    echo "- Regions where read alignments suggest incorrect joins" >> $REPORT_FILE
    echo "- Areas with unexpectedly low or high coverage" >> $REPORT_FILE
    echo "- Discrepancies between the assembly and the raw reads" >> $REPORT_FILE
    echo -e "\nA high-quality assembly should have few mis-assemblies, particularly in gene-rich regions.\n" >> $REPORT_FILE
elif [ -f "$MISASM_DIR/hifi_to_assembly.sam" ]; then
    echo "Basic alignment validation with minimap2 was performed." >> $REPORT_FILE
    echo "SAM file: $MISASM_DIR/hifi_to_assembly.sam" >> $REPORT_FILE
    
    echo -e "\nDetailed mis-assembly analysis was not possible due to tool unavailability, but basic read alignment was performed for validation.\n" >> $REPORT_FILE
else
    echo "Mis-assembly detection results not available." >> $REPORT_FILE
fi

# Add summary and potential improvements
echo -e "\n## 5. Summary and Potential Improvements\n" >> $REPORT_FILE
echo "### Summary\n" >> $REPORT_FILE
echo "The lizard genome assembly produced by Hifiasm shows promising characteristics:" >> $REPORT_FILE
echo "- A small number of contigs (88 in the primary assembly) suggests good contiguity" >> $REPORT_FILE
echo "- The assembly size appears appropriate for a lizard genome" >> $REPORT_FILE
echo "- Successful haplotype resolution resulted in two additional haplotype assemblies" >> $REPORT_FILE
echo -e "\n### Potential Improvements\n" >> $REPORT_FILE
echo "Based on the evaluation results, potential ways to improve the assembly include:" >> $REPORT_FILE
echo "1. **Hi-C Integration**: Incorporating Hi-C data for chromosome-level scaffolding" >> $REPORT_FILE
echo "2. **Polishing**: Using additional rounds of polishing with raw reads to improve base accuracy" >> $REPORT_FILE
echo "3. **Gap Filling**: Using long reads to fill any gaps in the assembly" >> $REPORT_FILE
echo "4. **Mis-assembly Correction**: Addressing any identified mis-assemblies" >> $REPORT_FILE
echo "5. **Manual Curation**: Manually reviewing and correcting problematic regions" >> $REPORT_FILE

log_message "Evaluation report generated: $REPORT_FILE"
log_message "====== Assembly Evaluation Completed: $(date) ======"

echo "====== Script completed at: $(date) ======"
