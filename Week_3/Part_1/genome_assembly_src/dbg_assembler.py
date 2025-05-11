# genome_assembly_src/dbg_assembler.py
from collections import defaultdict
from tqdm import tqdm
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genome_assembly_src.fastq_parser import read_fastq, write_fasta
from genome_assembly_src.graph_utils import write_gfa

class DBGAssembler:
    def __init__(self, k):
        self.k = k
        self.graph = defaultdict(list)
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)
        
    def build_graph(self, reads):
        """Build De Bruijn graph from reads."""
        print(f"Building De Bruijn graph with k={self.k}...")
        
        # Extract k-mers and build graph
        for _, seq, _ in tqdm(reads):
            if len(seq) < self.k:
                continue
                
            for i in range(len(seq) - self.k + 1):
                kmer = seq[i:i+self.k]
                # Check for non-ACGT characters
                if any(c not in 'ACGT' for c in kmer):
                    continue
                    
                prefix = kmer[:-1]
                suffix = kmer[1:]
                
                # Add edge
                if suffix not in self.graph[prefix]:
                    self.graph[prefix].append(suffix)
                    self.out_degree[prefix] += 1
                    self.in_degree[suffix] += 1
        
        print(f"Graph built with {len(self.graph)} nodes")
        return self.graph
    
    def find_contigs(self):
        """Find contigs by identifying paths in the graph."""
        print("Finding contigs...")
        
        contigs = []
        visited = set()
        
        # Start from nodes with imbalanced in/out degree
        start_nodes = []
        for node in self.graph:
            if self.out_degree[node] > 0 and (self.in_degree[node] != self.out_degree[node] or 
                                             self.in_degree[node] == 0):
                start_nodes.append(node)
        
        # If no suitable start nodes, use any node with outgoing edges
        if not start_nodes:
            for node in self.graph:
                if self.out_degree[node] > 0:
                    start_nodes.append(node)
                    break
        
        for start in start_nodes:
            if start in visited:
                continue
                
            # Start a new path
            path = [start]
            current = start
            
            while True:
                visited.add(current)
                
                # If current node has a unique next node and it's not a branching point
                next_nodes = self.graph[current]
                if len(next_nodes) == 1 and self.in_degree[next_nodes[0]] == 1:
                    next_node = next_nodes[0]
                    if next_node not in visited:
                        path.append(next_node[-1])  # Add just the last character
                        current = next_node
                    else:
                        break
                else:
                    break
            
            # Convert path to sequence
            contig = path[0] + ''.join(path[1:])
            if len(contig) >= self.k:
                contigs.append(contig)
        
        print(f"Found {len(contigs)} contigs")
        return contigs
    
    def assemble(self, fastq_file, output_dir="./"):
        """Assemble contigs from a FASTQ file."""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Read input file
        reads = read_fastq(fastq_file)
        print(f"Read {len(reads)} reads from {fastq_file}")
        
        # Build graph
        self.build_graph(reads)
        
        # Find contigs
        contigs = self.find_contigs()
        
        # Write contigs to FASTA
        basename = os.path.basename(fastq_file).split('.')[0]
        fasta_out = os.path.join(output_dir, f"{basename}_dbg_k{self.k}.fasta")
        write_fasta(fasta_out, contigs)
        
        # Write graph to GFA
        gfa_out = os.path.join(output_dir, f"{basename}_dbg_k{self.k}.gfa")
        write_gfa(self.graph, gfa_out)
        
        print(f"Assembly completed. Contigs written to {fasta_out}")
        print(f"Graph written to {gfa_out}")
        
        return fasta_out, gfa_out, contigs
