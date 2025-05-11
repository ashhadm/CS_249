# Lizard Genome Assembly and Evaluation (Task 2)

This repository contains scripts and results for the de novo assembly and evaluation of the Scincus mitranus (sand fish lizard) genome using PacBio HiFi reads. The assembly was performed using Hifiasm, and the quality was evaluated using multiple metrics.

## Overview

The Scincus mitranus genome is approximately 1.8 Gb in size. This project demonstrates the application of modern long-read assembly techniques to generate a high-quality reference genome. The assembly process involves:

1. PacBio HiFi-based assembly using Hifiasm
2. Evaluation using multiple quality assessment tools
3. Analysis of assembly metrics including contiguity, completeness, and base-level accuracy

## Data Sources

The data used for this assembly is from NCBI SRA under study accession PRJNA1221355:
- PacBio HiFi reads from lizard liver tissue (SRR32302813)

File location on Ibex cluster:
```
/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz
```

## Scripts and Workflow

### 1. Genome Assembly with Hifiasm

The script `hifiasm_assembly.sh` performs the assembly using Hifiasm:

```bash
sbatch hifiasm_assembly.sh
```

This script:
- Loads the Hifiasm module
- Runs assembly on PacBio HiFi reads
- Converts output GFA to FASTA format
- Creates a symbolic link to the final assembly

### 2. Assembly Evaluation

Multiple scripts were used for evaluation due to issues with the all-in-one approach:

#### Basic Metrics with QUAST

The script `evaluate_assembly.sh` runs QUAST for basic assembly metrics:

```bash
sbatch evaluate_assembly.sh
```

#### K-mer Analysis and QV Score with Merqury

The script `fast_merqury.sh` performs optimized Merqury analysis for QV score:

```bash
sbatch fast_merqury.sh
```

This script:
- Builds k-mer databases from HiFi reads
- Calculates k-mer completeness
- Computes the QV score (quality value)
- Generates a comprehensive report

#### BUSCO and Mis-assembly Detection Attempts

The script `busco_misasm.sh` attempts to run BUSCO and mis-assembly detection:

```bash
sbatch busco_misasm.sh
```

Note: We experienced challenges with BUSCO and Inspector/Flagger due to computational limitations and compatibility issues.

## Key Results

### Assembly Statistics (QUAST)

| Metric | Primary Assembly | Haplotype 1 | Haplotype 2 |
|--------|------------------|-------------|-------------|
| # contigs | 88 | 93 | 63 |
| Total length | 1,806,549,621 | 1,807,993,030 | 1,801,103,764 |
| Largest contig | 341,991,134 | 341,991,134 | 267,965,059 |
| GC (%) | 45.47 | 45.47 | 45.47 |
| N50 | 138,412,277 | 138,412,277 | 124,578,456 |
| L50 | 4 | 4 | 6 |
| # N's per 100 kbp | 0.00 | 0.00 | 0.00 |

### Merqury Evaluation

| Metric | Primary | Haplotype 1 | Haplotype 2 |
|--------|---------|-------------|-------------|
| K-mer completeness (%) | 97.3 | 97.1 | 96.8 |
| QV score | 41.2 | 40.9 | 40.7 |

## Analysis and Discussion

The Scincus mitranus genome assembly shows excellent quality metrics:

- **Exceptional contiguity**: The N50 of 138 Mb is remarkable, indicating that chromosomes are nearly fully assembled. Only 4 contigs are needed to reach 50% of the genome (L50).

- **High base-level accuracy**: The QV score of 41.2 exceeds the assignment target of 40, indicating approximately 1 error per 13,200 bases.

- **Completeness**: The k-mer completeness of 97.3% suggests that the assembly captures most of the genome.

- **Chromosome-scale contigs**: The largest contig spans an impressive 342 Mb, likely representing entire chromosomes or chromosome arms.

- **No gaps**: The assembly contains no artificial gap-filling sequences (0.00 N's per 100 kbp).

## Challenges Encountered

Several challenges were faced during the evaluation phase:

1. **BUSCO analysis failures**: Despite multiple attempts with different parameters and lineage datasets, BUSCO failed to complete successfully.

2. **Mis-assembly detection issues**: Inspector/Flagger faced memory limitations and compatibility issues with the large genome assembly.

3. **Computational resource requirements**: Processing a 1.8 Gb genome required substantial computational resources, particularly for memory-intensive tools like Merqury.

## Potential Improvements

Based on the evaluation, potential improvements to the assembly include:

1. **Hi-C Integration**: Incorporating Hi-C data for chromosome-level scaffolding
2. **Polishing**: Additional rounds of polishing with raw reads to improve base accuracy
3. **Gap Filling**: Using long reads to fill any remaining gaps
4. **Mis-assembly Correction**: Addressing any identified mis-assemblies
5. **Manual Curation**: Manually reviewing and correcting problematic regions

## Completion of Assignment Requirements

The requirements for Task 2 have been largely fulfilled:
- ✅ Task 2.1: Genome assembly using Hifiasm (completed)
- ✅ Task 2.2: Basic metrics with QUAST (completed)
- ✅ Task 2.2: QV score with Merqury (completed, QV > 40)
- ❌ Task 2.2: Gene completeness with BUSCO (unsuccessful despite multiple attempts)
- ❌ Task 2.2: Misassembly identification with Inspector (unsuccessful despite multiple attempts)

## Conclusion

The Hifiasm-based assembly of the Scincus mitranus genome has produced a high-quality reference with excellent contiguity and accuracy. The assembly meets the assignment requirement of a QV score greater than 40, despite challenges with some evaluation tools. This reference genome provides a solid foundation for future genetic and evolutionary studies of this desert-adapted reptile species.

## Repository Structure

```
lizard_assembly/
├── hifiasm_assembly.sh        # Script for genome assembly with Hifiasm
├── evaluate_assembly.sh       # Script for QUAST evaluation and report generation
├── fast_merqury.sh            # Script for optimized Merqury QV analysis
├── busco_misasm.sh            # Script for BUSCO and mis-assembly detection attempts
├── lizard_assembly_hifiasm/   # Output directory
│   ├── final_assembly.fasta   # Symbolic link to primary assembly
│   ├── lizard_hifiasm.*.gfa   # Graph files from Hifiasm
│   ├── lizard_hifiasm.*.fasta # Converted FASTA files
│   └── evaluation/            # Evaluation results
│       ├── quast/             # QUAST results
│       ├── fast_merqury/      # Merqury results
│       ├── busco/             # BUSCO attempts
│       └── misasm/            # Mis-assembly detection attempts
└── README.md                  # This file
```
