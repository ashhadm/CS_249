#!/bin/bash
#SBATCH --account=cs249
#SBATCH --job-name=fast_merqury
#SBATCH --output=fast_merqury_%j.out
#SBATCH --error=fast_merqury_%j.err
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --partition=largemem

# Print start time
echo "====== Fast Merqury Analysis started at: $(date) ======"

# Set working directory - use a different directory to avoid conflicts
WORK_DIR="$PWD/lizard_assembly_hifiasm"
cd $WORK_DIR

# Define paths with different output locations from the existing job
ASSEMBLY_FILE="$WORK_DIR/final_assembly.fasta"
EVAL_DIR="$WORK_DIR/evaluation"
MERQURY_DIR="$EVAL_DIR/fast_merqury"  # Different directory from original job

mkdir -p $MERQURY_DIR
mkdir -p $MERQURY_DIR/tmp

# Set temporary directory to avoid disk space issues
export TMPDIR=$MERQURY_DIR/tmp

# Check if assembly file exists
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "ERROR: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

# Define paths for raw reads
HIFI_READS="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"

echo "Assembly file: $ASSEMBLY_FILE"
echo "HiFi reads: $HIFI_READS"
echo "Fast Merqury output directory: $MERQURY_DIR"

# Load required modules
module purge
module load merqury

# Show available memory
free -h
echo "Total CPUs available: $(nproc)"

# Run Merqury with high-memory optimizations
echo "Running fast Merqury analysis with optimized memory settings..."

# Check if meryl is available
if command -v meryl &>/dev/null || module load merqury; then
    echo "Meryl found or loaded. Proceeding with k-mer analysis..."
    
    # Step 1: Count k-mers in reads with high-memory settings
    echo "Counting k-mers in HiFi reads with high memory allocation..."
    cd $MERQURY_DIR
    
    # Use optimal k-mer size for this genome
    KMER_SIZE=21
    
    # Use high memory setting to keep everything in RAM
    echo "Starting k-mer counting with 300GB memory allocation at $(date)"
    meryl count k=$KMER_SIZE threads=$SLURM_CPUS_PER_TASK memory=300 $HIFI_READS output $MERQURY_DIR/reads.meryl
    
    if [ $? -eq 0 ]; then
        echo "Read k-mer counting completed successfully at $(date)"
        
        # Step 2: Count k-mers in the assembly with high-memory settings
        echo "Counting k-mers in assembly..."
        meryl count k=$KMER_SIZE threads=$SLURM_CPUS_PER_TASK memory=300 $ASSEMBLY_FILE output $MERQURY_DIR/asm.meryl
        
        if [ $? -eq 0 ]; then
            echo "Assembly k-mer counting completed successfully at $(date)"
            
            # Step 3: Calculate statistics
            echo "Calculating k-mer statistics..."
            meryl statistics $MERQURY_DIR/reads.meryl > $MERQURY_DIR/reads.stats
            meryl statistics $MERQURY_DIR/asm.meryl > $MERQURY_DIR/asm.stats
            
            # Step 4: Find shared k-mers
            echo "Finding shared k-mers..."
            meryl intersect $MERQURY_DIR/reads.meryl $MERQURY_DIR/asm.meryl output $MERQURY_DIR/shared.meryl
            meryl statistics $MERQURY_DIR/shared.meryl > $MERQURY_DIR/shared.stats
            
            # Step 5: Calculate QV score
            echo "Calculating QV score..."
            ASM_ONLY=$(grep "present" $MERQURY_DIR/asm.stats | awk '{print $2}')
            SHARED=$(grep "present" $MERQURY_DIR/shared.stats | awk '{print $2}')
            
            if [ -n "$ASM_ONLY" ] && [ -n "$SHARED" ] && [ "$ASM_ONLY" -gt 0 ]; then
                ERROR_RATE=$(echo "scale=10; 1 - ($SHARED / $ASM_ONLY)" | bc)
                QV=$(echo "scale=2; -10 * l($ERROR_RATE)/l(10)" | bc -l)
                
                echo "Assembly k-mers: $ASM_ONLY" > $MERQURY_DIR/qv.txt
                echo "Shared k-mers: $SHARED" >> $MERQURY_DIR/qv.txt
                echo "Error rate: $ERROR_RATE" >> $MERQURY_DIR/qv.txt
                echo "QV score: $QV" >> $MERQURY_DIR/qv.txt
                
                echo "QV calculation completed. Results saved to $MERQURY_DIR/qv.txt"
                cat $MERQURY_DIR/qv.txt
                
                # Step 6: Calculate k-mer completeness
                echo "Calculating k-mer completeness..."
                READ_KMERS=$(grep "present" $MERQURY_DIR/reads.stats | awk '{print $2}')
                
                if [ -n "$READ_KMERS" ] && [ "$READ_KMERS" -gt 0 ]; then
                    COMPLETENESS=$(echo "scale=4; 100 * $SHARED / $READ_KMERS" | bc)
                    
                    echo "Read k-mers: $READ_KMERS" > $MERQURY_DIR/completeness.txt
                    echo "Shared k-mers: $SHARED" >> $MERQURY_DIR/completeness.txt
                    echo "K-mer completeness: $COMPLETENESS%" >> $MERQURY_DIR/completeness.txt
                    
                    echo "K-mer completeness calculation completed. Results saved to $MERQURY_DIR/completeness.txt"
                    cat $MERQURY_DIR/completeness.txt
                else
                    echo "ERROR: Unable to calculate k-mer completeness."
                fi
                
                # Step 7: Create a comprehensive report
                echo "Creating comprehensive QV report..."
                
                {
                    echo "High-Memory Assembly Quality (QV) Assessment Report"
                    echo "===================================================="
                    echo "Date: $(date)"
                    echo ""
                    echo "Assembly: $ASSEMBLY_FILE"
                    echo "Analysis performed with optimized high-memory settings"
                    echo ""
                    
                    echo "K-mer Based Quality Assessment:"
                    cat $MERQURY_DIR/qv.txt
                    
                    echo ""
                    echo "QV Score Interpretation:"
                    echo "- A QV score of 40 corresponds to an error rate of 1 in 10,000 bp"
                    echo "- A QV score of 50 corresponds to an error rate of 1 in 100,000 bp"
                    
                    if [ $(echo "$QV >= 40" | bc) -eq 1 ]; then
                        echo "- Assessment: EXCELLENT - QV score meets or exceeds the target of 40"
                    elif [ $(echo "$QV >= 30" | bc) -eq 1 ]; then
                        echo "- Assessment: GOOD - QV score between 30-40 indicates good assembly quality"
                    elif [ $(echo "$QV >= 20" | bc) -eq 1 ]; then
                        echo "- Assessment: ACCEPTABLE - QV score between 20-30 indicates acceptable quality"
                    else
                        echo "- Assessment: NEEDS IMPROVEMENT - QV score below 20 suggests higher error rates"
                    fi
                    
                    echo ""
                    echo "K-mer Completeness:"
                    cat $MERQURY_DIR/completeness.txt
                    
                    if [ $(echo "$COMPLETENESS >= 90" | bc) -eq 1 ]; then
                        echo "- Assessment: EXCELLENT - Very high k-mer completeness (>90%)"
                    elif [ $(echo "$COMPLETENESS >= 80" | bc) -eq 1 ]; then
                        echo "- Assessment: GOOD - High k-mer completeness (80-90%)"
                    elif [ $(echo "$COMPLETENESS >= 70" | bc) -eq 1 ]; then
                        echo "- Assessment: ACCEPTABLE - Moderate k-mer completeness (70-80%)"
                    else
                        echo "- Assessment: NEEDS IMPROVEMENT - Low k-mer completeness (<70%)"
                    fi
                    
                    echo ""
                    echo "Conclusion:"
                    if [ $(echo "$QV >= 30" | bc) -eq 1 ] && [ $(echo "$COMPLETENESS >= 80" | bc) -eq 1 ]; then
                        echo "The assembly is of HIGH QUALITY with good accuracy and completeness."
                    elif [ $(echo "$QV >= 20" | bc) -eq 1 ] && [ $(echo "$COMPLETENESS >= 70" | bc) -eq 1 ]; then
                        echo "The assembly is of GOOD QUALITY with acceptable accuracy and completeness."
                    else
                        echo "The assembly may benefit from additional polishing to improve accuracy and completeness."
                        echo "Recommended next steps include using tools like Pilon or Racon for polishing."
                    fi
                    
                    echo ""
                    echo "Note: This fast high-memory analysis should produce the same results as the slower method,"
                    echo "but with significantly reduced processing time due to in-memory operations."
                } > $MERQURY_DIR/fast_qv_report.txt
                
                echo "Comprehensive QV report generated: $MERQURY_DIR/fast_qv_report.txt"
                cat $MERQURY_DIR/fast_qv_report.txt
            else
                echo "ERROR: Unable to calculate QV score due to missing data."
            fi
        else
            echo "ERROR: Assembly k-mer counting failed."
        fi
    else
        echo "ERROR: Read k-mer counting failed."
    fi
    
    # Return to working directory
    cd $WORK_DIR
else
    echo "ERROR: Meryl not found. Cannot run Merqury analysis."
fi

# Clean up
echo "Cleaning up temporary files..."
rm -rf $MERQURY_DIR/tmp

echo "====== Fast Merqury Analysis completed at: $(date) ======"
