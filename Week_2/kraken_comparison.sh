#!/bin/bash
# kraken_comparison.sh - Script to compare custom implementations with Kraken2

# Get the path to kraken2 tools
KRAKEN2_BIN=$(which kraken2)
KRAKEN_BUILD=$(which kraken2-build)

# Set paths
DB_DIR="./kraken2_custom_db"
OUTPUT_DIR="./kraken2_results"

# Create output directories
mkdir -p $DB_DIR
mkdir -p $OUTPUT_DIR

echo "================ KRAKEN2 COMPARISON TASK ================"
echo "Using Kraken2 at: $KRAKEN2_BIN"

# 1. Download taxonomy (required before building)
echo "Downloading NCBI taxonomy..."
/usr/bin/time -v $KRAKEN_BUILD --download-taxonomy --db $DB_DIR 2>&1 | tee $OUTPUT_DIR/taxonomy_download.log

# 2. Add each genome to the library (one at a time)
echo "Adding genomes to library..."
/usr/bin/time -v $KRAKEN_BUILD --add-to-library ./e_coli.fna --db $DB_DIR 2>&1 | tee -a $OUTPUT_DIR/add_library.log
/usr/bin/time -v $KRAKEN_BUILD --add-to-library ./b_subtilis.fna --db $DB_DIR 2>&1 | tee -a $OUTPUT_DIR/add_library.log
/usr/bin/time -v $KRAKEN_BUILD --add-to-library ./m_tuberculosis.fna --db $DB_DIR 2>&1 | tee -a $OUTPUT_DIR/add_library.log
/usr/bin/time -v $KRAKEN_BUILD --add-to-library ./p_aeruginosa.fna --db $DB_DIR 2>&1 | tee -a $OUTPUT_DIR/add_library.log
/usr/bin/time -v $KRAKEN_BUILD --add-to-library ./s_aureus.fna --db $DB_DIR 2>&1 | tee -a $OUTPUT_DIR/add_library.log

# 3. Build the database with longer k-mers and minimizers to make it more exact-matching like
# Using maximum k-mer size and minimizer size for more exact-matching behavior
echo "Building Kraken2 database with stricter matching parameters..."
/usr/bin/time -v $KRAKEN_BUILD --build --threads 8 \
  --kmer-len 35 --minimizer-len 35 --minimizer-spaces 0 \
  --db $DB_DIR 2>&1 | tee $OUTPUT_DIR/db_build.log

echo "Database building complete!"

# 4. Process each FASTQ file separately (not paired)
# For error-free reads, R1
echo "Running Kraken2 on error-free reads (R1)..."
/usr/bin/time -v $KRAKEN2_BIN \
  --db $DB_DIR \
  --threads 8 \
  --output $OUTPUT_DIR/error_free_R1_results.txt \
  --report $OUTPUT_DIR/error_free_R1_report.txt \
  --use-names \
  --confidence 1.0 \
  --quick \
  --no-unclassified-out \
  ./simulated_reads_no_errors_10k_R1.fastq \
  2>&1 | tee $OUTPUT_DIR/error_free_R1_classification.log

# For error-free reads, R2
echo "Running Kraken2 on error-free reads (R2)..."
/usr/bin/time -v $KRAKEN2_BIN \
  --db $DB_DIR \
  --threads 8 \
  --output $OUTPUT_DIR/error_free_R2_results.txt \
  --report $OUTPUT_DIR/error_free_R2_report.txt \
  --use-names \
  --confidence 1.0 \
  --quick \
  --no-unclassified-out \
  ./simulated_reads_no_errors_10k_R2.fastq \
  2>&1 | tee $OUTPUT_DIR/error_free_R2_classification.log

# For reads with errors, R1
echo "Running Kraken2 on reads with errors (R1)..."
/usr/bin/time -v $KRAKEN2_BIN \
  --db $DB_DIR \
  --threads 8 \
  --output $OUTPUT_DIR/error_reads_R1_results.txt \
  --report $OUTPUT_DIR/error_reads_R1_report.txt \
  --use-names \
  --confidence 1.0 \
  --quick \
  --no-unclassified-out \
  ./simulated_reads_miseq_10k_R1.fastq \
  2>&1 | tee $OUTPUT_DIR/error_reads_R1_classification.log

# For reads with errors, R2
echo "Running Kraken2 on reads with errors (R2)..."
/usr/bin/time -v $KRAKEN2_BIN \
  --db $DB_DIR \
  --threads 8 \
  --output $OUTPUT_DIR/error_reads_R2_results.txt \
  --report $OUTPUT_DIR/error_reads_R2_report.txt \
  --use-names \
  --confidence 1.0 \
  --quick \
  --no-unclassified-out \
  ./simulated_reads_miseq_10k_R2.fastq \
  2>&1 | tee $OUTPUT_DIR/error_reads_R2_classification.log

echo "Classification complete!"

# 5. Combine results for comparison
echo "Combining results from all files..."
echo "Error-free reads (both R1 and R2):"
cat $OUTPUT_DIR/error_free_R1_results.txt $OUTPUT_DIR/error_free_R2_results.txt > $OUTPUT_DIR/error_free_combined_results.txt
cat $OUTPUT_DIR/error_free_R1_report.txt $OUTPUT_DIR/error_free_R2_report.txt > $OUTPUT_DIR/error_free_combined_report.txt

echo "Error-containing reads (both R1 and R2):"
cat $OUTPUT_DIR/error_reads_R1_results.txt $OUTPUT_DIR/error_reads_R2_results.txt > $OUTPUT_DIR/error_reads_combined_results.txt
cat $OUTPUT_DIR/error_reads_R1_report.txt $OUTPUT_DIR/error_reads_R2_report.txt > $OUTPUT_DIR/error_reads_combined_report.txt

# 6. Extract and display key metrics
echo "================ PERFORMANCE SUMMARY ================"
echo "Database Building:"
grep "Maximum resident set size" $OUTPUT_DIR/db_build.log | awk '{printf "Peak Memory: %.2f MB\n", $6/1024}'
grep "Elapsed (wall clock) time" $OUTPUT_DIR/db_build.log | awk -F: '{print "Wall time:", $0}'

echo "Error-free Reads Classification (total):"
grep "Maximum resident set size" $OUTPUT_DIR/error_free_R1_classification.log $OUTPUT_DIR/error_free_R2_classification.log | \
  awk '{sum += $6} END {printf "Peak Memory: %.2f MB\n", sum/1024}'
grep "Elapsed (wall clock) time" $OUTPUT_DIR/error_free_R1_classification.log $OUTPUT_DIR/error_free_R2_classification.log | \
  awk -F: '{print "Wall time:", $0}'

echo "Error-containing Reads Classification (total):"
grep "Maximum resident set size" $OUTPUT_DIR/error_reads_R1_classification.log $OUTPUT_DIR/error_reads_R2_classification.log | \
  awk '{sum += $6} END {printf "Peak Memory: %.2f MB\n", sum/1024}'
grep "Elapsed (wall clock) time" $OUTPUT_DIR/error_reads_R1_classification.log $OUTPUT_DIR/error_reads_R2_classification.log | \
  awk -F: '{print "Wall time:", $0}'

# 7. Parse and display combined classification results
echo "================ CLASSIFICATION SUMMARY ================"

# Count classified and unclassified reads across all files
error_free_total=$(cat $OUTPUT_DIR/error_free_combined_results.txt | wc -l)
error_free_classified=$(grep -c "^C" $OUTPUT_DIR/error_free_combined_results.txt)
error_free_unclassified=$(grep -c "^U" $OUTPUT_DIR/error_free_combined_results.txt)

error_reads_total=$(cat $OUTPUT_DIR/error_reads_combined_results.txt | wc -l)
error_reads_classified=$(grep -c "^C" $OUTPUT_DIR/error_reads_combined_results.txt)
error_reads_unclassified=$(grep -c "^U" $OUTPUT_DIR/error_reads_combined_results.txt)

# Display statistics
echo "Error-free Reads Results (Combined R1+R2):"
echo "Total reads processed: $error_free_total"
echo "Classified reads: $error_free_classified ($(echo "scale=2; $error_free_classified*100/$error_free_total" | bc -l)%)"
echo "Unclassified reads: $error_free_unclassified ($(echo "scale=2; $error_free_unclassified*100/$error_free_total" | bc -l)%)"

echo "Error-containing Reads Results (Combined R1+R2):"
echo "Total reads processed: $error_reads_total"
echo "Classified reads: $error_reads_classified ($(echo "scale=2; $error_reads_classified*100/$error_reads_total" | bc -l)%)"
echo "Unclassified reads: $error_reads_unclassified ($(echo "scale=2; $error_reads_unclassified*100/$error_reads_total" | bc -l)%)"

echo "=================================================="
echo "See detailed results in the $OUTPUT_DIR directory"
