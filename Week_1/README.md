# Week 1 (Assignment 0)
# Genome Pattern Matching Tool

## Overview
This tool searches for exact and approximate matches of a given pattern in a genomic dataset. It supports large genome files in compressed formats and utilizes parallel processing for efficiency.

## Features
- **Exact pattern matching** using the Knuth-Morris-Pratt (KMP) algorithm.
- **Approximate pattern matching** using an optimized edit distance function (`quick_distance`).
  - Allows up to one mismatch, including **substitution, insertion, or deletion**.
  - Uses early termination to improve efficiency, stopping when more than one mismatch is found.
  - Efficiently determines the edit type: **'S' for substitution, 'I' for insertion, and 'D' for deletion**.
- **Multi-threaded processing** for large genome datasets.
- **Support for compressed and zipped genome files**.
- **Performance tracking** with runtime and memory usage statistics.

## Requirements
Make sure you have the following dependencies installed:

- Python 3.x
- Biopython (`pip install biopython`)
- GNU time (for profiling)

## Installation
Clone this repository and navigate into the project folder:
```bash
git clone https://github.com/ashhadm/CS_249.git
cd Week_1

```

## Usage
Run the script using the provided shell script and provide the path to the genome and the pattern to be searched:
```bash
./run_search.sh <genome_file.zip> <pattern.fna> <num_cores>
```

### Example:
```bash
./run_search.sh genome_data.zip pattern.fna 4
```

This will:
- Search for the pattern in `genome_data.zip`.
- Use `pattern.fna` as the query pattern.
- Utilize `4` CPU cores for faster processing.
- Store the results in `output.log`.

## Output
- **Exact and mismatched matches** will be stored in `all_matches.txt`.
- **Profiling data** (runtime, memory usage) will be stored in `output.log`.

## Script Details
### `genome_search.py`
This script performs the following operations:
1. Reads the genome file (supports ZIP and GZ formats).
2. Reads the pattern from the `.fna` file.
3. Uses the KMP algorithm for exact matches.
4. Identifies mismatched matches (allowing up to one edit distance) using the `quick_distance` function.
5. Outputs detailed results in `all_matches.txt`.

### `run_search.sh`
This shell script automates execution and profiling:
- Takes three arguments: `genome_file.zip`, `pattern.fna`, and `num_cores`.
- Runs the Python script and captures performance metrics.
- Logs results into `output.log`.

## License
This project is licensed under the MIT License.

## Author
Mohammad Ashhad

For any issues or contributions, feel free to submit a pull request or open an issue.
