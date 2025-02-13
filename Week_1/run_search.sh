#!/bin/bash

# Ensure we exit the script on error
set -e

# Get inputs from command line arguments
GENOME_ZIP_PATH=$1
PATTERN_PATH=$2
NUM_CORES=$3
OUTPUT_FILE="output.log"

# Measure time and memory usage using the `time` command
# Running the Python script
echo "Starting pattern matching with $NUM_CORES cores..."

# Time and memory profiling
/usr/bin/time -v python3 genome_search.py "$GENOME_ZIP_PATH" "$PATTERN_PATH" "$NUM_CORES" > "$OUTPUT_FILE" 2>&1

# Output the results
echo "Pattern matching finished. Results are stored in $OUTPUT_FILE."

