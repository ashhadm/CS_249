# genome_assembly_src/olc_assembler.py - Improved version
import networkx as nx
from tqdm import tqdm
import os
import sys
import numpy as np

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genome_assembly_src.fastq_parser import read_fastq, write_fasta

class OLCAssembler:
    def __init__(self, min_overlap):
        self.min_overlap = min_overlap
        self.graph = nx.DiGraph()
        
    def find_overlaps(self, reads):
        """Find all overlaps between reads."""
        print(f"Finding overlaps with minimum length {self.min_overlap}...")
        
        # Adjust read limit based on read count - process all reads for datasets with <500 reads
        read_count = len(reads)
        max_reads = min(read_count, 2000) if read_count < 500 else 500
        
        if read_count > max_reads:
            print(f"Processing {max_reads} out of {read_count} reads to save time")
        
        # Process the first max_reads
        working_reads = reads[:max_reads]
        
        # For ONT data, reads are typically longer - check if average read length > 1000
        avg_len = np.mean([len(seq) for _, seq, _ in working_reads[:100]])
        is_long_read = avg_len > 1000
        
        # Adapt min_overlap for long reads
        effective_min_overlap = self.min_overlap * 3 if is_long_read else self.min_overlap
        if is_long_read:
            print(f"Detected long reads (avg length: {avg_len:.1f}), using effective min overlap: {effective_min_overlap}")
        
        # Build overlap graph
        for i, (header_i, seq_i, _) in enumerate(tqdm(working_reads)):
            for j, (header_j, seq_j, _) in enumerate(working_reads):
                if i == j:
                    continue
                
                # Find the longest overlap
                max_overlap = 0
                # For long reads, we can use a step size to speed up search
                step = 10 if is_long_read else 1
                for k in range(effective_min_overlap, min(len(seq_i), len(seq_j)) + 1, step):
                    if seq_i[-k:] == seq_j[:k]:
                        # For step > 1, refine to find exact overlap
                        if step > 1:
                            # Search for exact overlap around k
                            start = max(k - step, effective_min_overlap)
                            end = min(k + step, min(len(seq_i), len(seq_j)))
                            for exact_k in range(start, end + 1):
                                if seq_i[-exact_k:] == seq_j[:exact_k]:
                                    max_overlap = max(max_overlap, exact_k)
                        else:
                            max_overlap = k
                        break  # Found the longest overlap starting at this position
                
                if max_overlap >= effective_min_overlap:
                    # Add edge with weight = overlap
                    self.graph.add_edge(header_i, header_j, weight=max_overlap)
        
        print(f"Found {self.graph.number_of_edges()} overlaps")
        return self.graph
    
    def find_paths(self):
        """Find paths in the graph using improved algorithm."""
        print("Finding paths in the overlap graph...")
        
        if not self.graph.nodes():
            print("Warning: Graph is empty, no paths to find.")
            return []
        
        paths = []
        visited_edges = set()  # Track visited edges instead of nodes
        
        # Find potential start nodes: prefer nodes with no incoming edges
        start_nodes = [n for n in self.graph.nodes() if self.graph.in_degree(n) == 0]
        
        # If no such nodes, use nodes with fewer in-edges than out-edges 
        if not start_nodes:
            start_nodes = [n for n in self.graph.nodes() 
                          if self.graph.in_degree(n) < self.graph.out_degree(n)]
        
        # If still no nodes, use nodes with highest out-degree
        if not start_nodes:
            nodes_with_edges = [(n, self.graph.out_degree(n)) for n in self.graph.nodes() 
                               if self.graph.out_degree(n) > 0]
            if nodes_with_edges:
                nodes_with_edges.sort(key=lambda x: x[1], reverse=True)
                start_nodes = [nodes_with_edges[0][0]]
            else:
                # If no nodes have outgoing edges, just use any node
                start_nodes = list(self.graph.nodes())[:1]
        
        # Process each start node
        for start in start_nodes:
            # Build paths greedily from this start node
            self._build_greedy_paths(start, paths, visited_edges)
        
        # Look for remaining unvisited edges and build paths from them
        all_edges = set((u, v) for u, v in self.graph.edges())
        unvisited_edges = all_edges - visited_edges
        
        while unvisited_edges:
            # Pick an unvisited edge and build a path from it
            u, v = next(iter(unvisited_edges))
            self._build_greedy_paths(u, paths, visited_edges, start_edge=(u, v))
            
            # Update unvisited edges
            unvisited_edges = all_edges - visited_edges
        
        print(f"Found {len(paths)} paths")
        return paths
    
    def _build_greedy_paths(self, start, paths, visited_edges, start_edge=None):
        """
        Build paths greedily from a starting node.
        
        Args:
            start: Starting node
            paths: List to collect paths
            visited_edges: Set to track visited edges
            start_edge: Optional specific edge to start with
        """
        # Start a new path
        current = start
        path = [current]
        
        # If we have a specific start edge, add its target to the path
        if start_edge and start_edge[0] == start:
            visited_edges.add(start_edge)
            next_node = start_edge[1]
            path.append(next_node)
            current = next_node
        
        # Extend the path greedily
        while True:
            # Get all outgoing edges
            out_edges = list(self.graph.out_edges(current))
            
            # Filter out visited edges
            unvisited_out = [(u, v) for u, v in out_edges if (u, v) not in visited_edges]
            
            if not unvisited_out:
                break  # No more unvisited outgoing edges
            
            # Choose edge with maximum overlap weight
            next_edge = max(unvisited_out, 
                           key=lambda e: self.graph[e[0]][e[1]]['weight'])
            
            visited_edges.add(next_edge)
            next_node = next_edge[1]
            
            # Avoid cycles in the path
            if next_node in path:
                break
                
            path.append(next_node)
            current = next_node
        
        # Add the path if it's long enough
        if len(path) > 1:
            paths.append(path)
    
    def generate_consensus(self, paths, read_dict):
        """Generate consensus sequences from paths."""
        print("Generating consensus sequences...")
        
        contigs = []
        
        for path in paths:
            if len(path) == 1:
                # Single read
                contigs.append(read_dict[path[0]])
                continue
                
            # Construct contig from path
            contig = read_dict[path[0]]
            
            for i in range(1, len(path)):
                # Get overlap length
                overlap = self.graph[path[i-1]][path[i]]['weight']
                
                # Make sure we don't get index errors
                if overlap >= len(read_dict[path[i]]):
                    print(f"Warning: Overlap {overlap} is >= read length {len(read_dict[path[i]])}.")
                    overlap = max(0, len(read_dict[path[i]]) - 1)
                
                # Append non-overlapping part
                contig += read_dict[path[i]][overlap:]
            
            contigs.append(contig)
        
        # Handle the case of no paths by using the longest reads
        if not contigs and read_dict:
            print("Warning: No paths found. Using longest reads as contigs.")
            # Get the 5 longest reads as contigs
            longest_reads = sorted([(h, len(s)) for h, s in read_dict.items()], 
                                 key=lambda x: x[1], reverse=True)[:5]
            for header, _ in longest_reads:
                contigs.append(read_dict[header])
        
        print(f"Generated {len(contigs)} contigs")
        return contigs
    
    def assemble(self, fastq_file, output_dir="./"):
        """Assemble contigs from a FASTQ file."""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Read input file
        reads = read_fastq(fastq_file)
        print(f"Read {len(reads)} reads from {fastq_file}")
        
        # Create read dictionary
        read_dict = {header: seq for header, seq, _ in reads}
        
        # Find overlaps
        self.find_overlaps(reads)
        
        # Find paths
        paths = self.find_paths()
        
        # Generate consensus
        contigs = self.generate_consensus(paths, read_dict)
        
        # Write contigs to FASTA
        basename = os.path.basename(fastq_file).split('.')[0]
        fasta_out = os.path.join(output_dir, f"{basename}_olc_o{self.min_overlap}.fasta")
        write_fasta(fasta_out, contigs)
        
        print(f"Assembly completed. Contigs written to {fasta_out}")
        
        return fasta_out, contigs
