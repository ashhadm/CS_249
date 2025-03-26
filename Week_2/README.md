# CS249 - Assignment 1: k-mer Index with Minimizer

This repository contains the implementation and analysis of various metagenomic classification approaches for CS249 Assignment 1. The project examines different string matching and indexing techniques for classifying DNA sequence reads against a reference database of bacterial genomes.

## Project Overview

In this assignment, I implemented and compared four different approaches for metagenomic classification:

1. **Suffix Array-Based Exact Matching**: Implementation of a memory-efficient suffix array approach for exact string matching
2. **K-mer Index-Based Classification**: Implementation of a k-mer index (k=31) for fast lookup of sequence reads
3. **Minimizer-Based Classification**: Optimization of the k-mer approach using a minimizer scheme to reduce memory usage
4. **BLAST Integration**: Comparison with BLAST (Basic Local Alignment Search Tool) configured for exact matching

Additionally, I performed a real-world use case analysis comparing human gut and wastewater metagenomic samples using Kraken2.

## Repository Structure

```
.
├── suffix_array_classification.ipynb  # Task 1.2: Suffix array implementation
├── blast_classification.ipynb         # Task 1.4: BLAST comparison
├── kmer_index_classification.ipynb    # Task 2.1-2.2: K-mer index implementation
├── minimizer_classification.ipynb     # Task 2.3: Minimizer implementation
├── kraken2_comparison.sh              # Task 3.1: Kraken2 comparison script
├── kraken2_realworld.sh               # Task 3.2: Real-world samples script
├── metagenomic_analysis.ipynb         # Task 3.2: Analysis of real-world data
├── report.pdf                         # Full assignment report
└── README.md                          # This readme file
```

## Implementation Details

### Task 1: Metagenome Classification by String Matching

#### Suffix Array Implementation (`suffix_array_classification.ipynb`)
- Memory-efficient implementation of suffix arrays for each reference genome
- Binary search algorithm for exact pattern matching
- LCA (Lowest Common Ancestor) approach for handling multi-matching reads
- Performance metrics tracking (time, memory, classification rate)

#### BLAST Comparison (`blast_classification.ipynb`)
- Integration with BLASTN for exact matching
- Configuration for forward-strand search with 100% identity requirement
- Performance comparison with the suffix array implementation
- Analysis of classification discrepancies

### Task 2: Metagenomic Classification by k-mer Index

#### K-mer Index Implementation (`kmer_index_classification.ipynb`)
- Dictionary-based k-mer index using 2-bit encoding for DNA bases
- Multi-process k-mer extraction for improved performance
- Analysis of k-mer distribution across reference genomes
- Memory usage and performance benchmarking

#### Minimizer Implementation (`minimizer_classification.ipynb`)
- Window-based minimizer selection (window size w=10)
- Lexicographic ordering for minimizer selection
- Memory reduction analysis compared to full k-mer index
- Classification performance evaluation with error-containing reads

### Task 3: Real-world Data and Tools

#### Kraken2 Comparison (`kraken2_comparison.sh`)
- Custom Kraken2 database with the 5 reference genomes
- Optimized parameters for more exact-like matching
- Performance measurements (time, memory, classification rate)
- Comparison with custom implementations

#### Real-world Use Case (`kraken2_realworld.sh` and `metagenomic_analysis.ipynb`)
- Analysis of human gut (SRR11412*) and wastewater (SRR21907*) samples
- Taxonomic profiling at genus level using Kraken2
- Statistical analysis (PCA, hierarchical clustering, discriminating taxa)
- Environmental microbial signature identification

## How to Run the Code

### Python Notebooks
All Python implementations can be run as Jupyter notebooks:

```bash
# Install dependencies
pip install numpy pandas matplotlib seaborn scipy scikit-learn psutil

# Run notebooks
jupyter notebook suffix_array_classification.ipynb
jupyter notebook blast_classification.ipynb
jupyter notebook kmer_index_classification.ipynb
jupyter notebook minimizer_classification.ipynb
jupyter notebook metagenomic_analysis.ipynb
```

### Shell Scripts
The Kraken2 comparison and real-world analysis scripts can be run as follows:

```bash
# Make scripts executable
chmod +x kraken2_comparison.sh
chmod +x kraken2_realworld.sh

# Run scripts
./kraken2_comparison.sh
./kraken2_realworld.sh
```

### Data Requirements
- Reference genomes for 5 bacterial species (E. coli, B. subtilis, P. aeruginosa, S. aureus, M. tuberculosis)
- Simulated read datasets (error-free and error-containing)
- For Task 3.2: Access to SRA data for the specified accessions

## Results Summary

### Classification Performance
| Method | Error-free Reads (%) | Error Reads (%) | Memory (MB) | Time (s) |
|--------|----------------------|-----------------|-------------|----------|
| Suffix Array | 50.68% | 3.76% | 1,567.41 | 41.00 |
| BLAST | 50.90% | 5.51% | 255.20 | 5.83 |
| K-mer Index | 50.71% | 3.77% | 8,685.21 | 154.50 |
| Minimizer | 50.71% | 5.07% | 1,917.81 | 103.78 |
| Kraken2 | 98.08% | 92.23% | 7,806.00* | 12.70 |

*For Kraken2, peak memory is during database building; classification used only 55-70 MB.

### Real-world Analysis
The analysis of human gut and wastewater samples revealed:
- Clear separation of samples by environment in both PCA and hierarchical clustering
- Betacoronavirus as the most discriminating taxon (abundant in wastewater, absent in gut)
- Human gut samples dominated by obligate anaerobes (Bacteroides, Prevotella, etc.)
- Wastewater samples characterized by environmentally resilient bacteria

## Dependencies

### Python Dependencies
- numpy
- pandas
- matplotlib
- seaborn
- scipy
- scikit-learn
- psutil
- multiprocessing

### External Tools
- NCBI BLAST+ (blastn, makeblastdb)
- Kraken2
- SRA Toolkit (for downloading SRA data)

## References

- BLAST: Altschul SF, et al. (1990). Basic local alignment search tool. J Mol Biol, 215(3):403-10.
- Kraken2: Wood DE, et al. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257.
- Minimizers: Roberts M, et al. (2004). A preprocessor for shotgun assembly of large genomes. J Comput Biol, 11(4):734-52.

---

For a detailed analysis and discussion of results, please refer to the [full report](report.pdf).
