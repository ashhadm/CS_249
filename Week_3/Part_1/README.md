# Genome Assembly Algorithms Implementation

This repository contains implementations of two fundamental genome assembly algorithms: De Bruijn Graph (DBG) and Overlap-Layout-Consensus (OLC). These algorithms are used to reconstruct complete genome sequences from fragmented DNA sequence reads.

## Overview

Genome assembly is a critical step in bioinformatics pipelines. This project implements:

1. **De Bruijn Graph (DBG) Assembler**: Constructs a graph from k-mers extracted from sequencing reads and finds contigs by traversing this graph.
2. **Overlap-Layout-Consensus (OLC) Assembler**: Identifies overlaps between reads, constructs an overlap graph, finds paths in this graph, and generates consensus sequences.

Additionally, the repository includes utilities for parsing FASTQ/FASTA files, graph operations, and assembly evaluation.

## Prerequisites

- Python 3.6+
- BioPython
- NetworkX
- tqdm
- pandas
- NumPy

For external evaluation:
- [QUAST](http://quast.sourceforge.net/quast): For assembly evaluation
- [Bandage](https://rrwick.github.io/Bandage/): For assembly graph visualization

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/genome-assembly.git
cd genome-assembly

# Install dependencies
pip install biopython networkx tqdm pandas numpy
```

## Repository Structure

```
genome-assembly/
├── genome_assembly_src/
│   ├── dbg_assembler.py      # De Bruijn Graph assembler implementation
│   ├── olc_assembler.py      # Overlap-Layout-Consensus assembler implementation
│   ├── fastq_parser.py       # Utilities for reading/writing FASTQ/FASTA files
│   ├── graph_utils.py        # Utilities for graph operations (GFA output)
│   └── evaluation.py         # Functions for assembly evaluation using QUAST
├── task1_3_part1.py          # Run DBG on reads_b.fastq with k=40
├── task1_3_part2.py          # Run DBG on reads_r.fastq with k=35 and k=45, compare
├── task1_3_part3.py          # Run DBG and OLC on MERS datasets
├── task1_3_part4.py          # Run SPAdes on datasets and compare with our implementations
├── run_task1_3.sh            # Script to run all Task 1.3 components
└── cs249_hw2/                     # Place your FASTQ files here
```

## Usage

### Basic Usage

The assemblers can be run directly as Python modules:

```bash
# De Bruijn Graph (DBG) Assembler
python -m genome_assembly_src.dbg_assembler /path/to/reads.fastq --k 35 --output ./results/

# Overlap-Layout-Consensus (OLC) Assembler
python -m genome_assembly_src.olc_assembler /path/to/reads.fastq --min_overlap 30 --output ./results/
```

### Running the Task Scripts

The repository contains scripts to run various assembly experiments:

```bash
# Run all tasks in sequence
bash run_task1_3.sh

# Or run individual tasks
python task1_3_part1.py  # DBG assembly of reads_b.fastq with k=40
python task1_3_part2.py  # Compare DBG assemblies with k=35 and k=45
python task1_3_part3.py  # Compare DBG and OLC on MERS datasets
python task1_3_part4.py  # Compare with SPAdes professional assembler
```

### Output

Each assembler generates:
- FASTA files containing assembled contigs
- GFA files representing the assembly graph (for DBG)
- Evaluation reports when run through QUAST

## Algorithm Details

### De Bruijn Graph (DBG) Assembler

The DBG assembler works by:
1. Extracting k-mers from input reads
2. Building a De Bruijn graph where:
   - Nodes are (k-1)-mers
   - Edges represent k-mers connecting consecutive (k-1)-mers
3. Finding contigs by identifying non-branching paths in the graph
4. Writing the results to FASTA and GFA files

Choosing an appropriate k-mer size is critical:
- Larger k: More specific assemblies but may create more fragmented results
- Smaller k: Increases connectivity but may introduce more ambiguities

### Overlap-Layout-Consensus (OLC) Assembler

The OLC assembler works by:
1. Finding overlaps between all pairs of reads (with minimum overlap length)
2. Building an overlap graph where:
   - Nodes are reads
   - Edges represent overlaps between reads
3. Finding paths in the graph to determine the layout of reads
4. Generating consensus sequences by merging reads according to overlaps
5. Reducing redundancy through containment removal and clustering

The OLC approach is particularly effective for longer reads where overlaps can be more confidently identified.

## Task Results

### Task 1.3 Part 1: De Bruijn Graph with k=40

This task builds a DBG assembly of `reads_b.fastq` with k=40 and exports the graph in GFA format for visualization with Bandage.

Key observations:
- The graph shows how reads connect to form contigs
- Visualization helps identify bubbles, branches, and other graph structures
- Graph interpretation can guide parameter tuning for improved assemblies

### Task 1.3 Part 2: Effect of k-mer Size

This task compares DBG assemblies of `reads_r.fastq` with k=35 and k=45.

Key differences:
- k=35 assembly generally shows more connectivity between nodes
- k=45 assembly typically has fewer ambiguities but may be more fragmented
- Smaller k-values can better handle low-coverage regions
- Larger k-values provide better resolution of repeats

### Task 1.3 Part 3: MERS Dataset Assembly

This task applies both DBG and OLC algorithms to MERS virus datasets with the following variations:
- Error-free vs. error-containing reads
- Short (Illumina HiSeq) vs. long (ONT) reads

Key findings:
- Both algorithms perform better on error-free reads
- DBG typically handles short reads well but is sensitive to errors
- OLC performs better with longer reads and may be more resilient to some errors
- The choice of algorithm should depend on the data type and quality

### Task 1.3 Part 4: Comparison with Professional Tools

This task compares our custom implementations with SPAdes, a state-of-the-art genome assembler.

Key advantages of SPAdes:
1. Uses multiple k-mer sizes to build a combined de Bruijn graph
2. Incorporates error correction using BayesHammer
3. Uses sophisticated algorithms for bubble removal and repeat resolution
4. Includes advanced graph simplification techniques
5. Can utilize paired-end and mate-pair information

## Visualization

Assembly graphs can be visualized using Bandage:

```bash
# For DBG assembly graphs
Bandage load results/reads_b_dbg_k40.gfa
```

Bandage allows exploration of the graph structure, helping to identify:
- Contigs (non-branching paths)
- Bubbles (alternative paths)
- Dead-ends (potential assembly errors)
- Repetitive regions

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
