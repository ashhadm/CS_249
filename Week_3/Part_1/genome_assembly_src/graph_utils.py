# genome_assembly_src/graph_utils.py
def write_gfa(graph, filename):
    """
    Write a De Bruijn graph to a GFA file for visualization.
    
    Args:
        graph: Dictionary representation of the graph
        filename: Output GFA file path
    """
    with open(filename, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments (nodes)
        for node in graph:
            f.write(f"S\t{node}\t{node}\n")
        
        # Write links (edges)
        edge_count = 0
        for node, edges in graph.items():
            for next_node in edges:
                edge_count += 1
                f.write(f"L\t{node}\t+\t{next_node}\t+\t{len(node)-1}M\n")
                
    print(f"GFA file written to {filename} with {len(graph)} nodes and {edge_count} edges")
    return filename
