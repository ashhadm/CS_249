#!/usr/bin/env python3
import os
import sys
import gzip
import io
import zipfile
import multiprocessing
from functools import partial
from time import time
from datetime import datetime
import threading
from Bio import SeqIO

# Global variables for synchronization and progress tracking
print_lock = threading.Lock()
progress_lock = threading.Lock()
processed_chromosomes = multiprocessing.Value('i', 0)
total_chromosomes = multiprocessing.Value('i', 0)
total_bases_processed = multiprocessing.Value('L', 0)
last_update_time = multiprocessing.Value('d', 0.0)

def compute_lps(pattern):
    """Compute the longest prefix suffix (LPS) array for KMP algorithm."""
    lps = [0] * len(pattern)
    length = 0
    i = 1
    
    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    
    return lps

def kmp_search(text, pattern):
    """Find exact matches of pattern in text using the KMP algorithm."""
    lps = compute_lps(pattern)
    matches = []
    
    i = j = 0
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1
        
        if j == len(pattern):
            matches.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
                
    return matches

def quick_distance(text, pattern):
    """Optimized edit distance that stops if more than 1 mismatch found."""
    if abs(len(text) - len(pattern)) > 1:
        return (2, None)
        
    edits = 0
    p_len = len(pattern)
    t_len = len(text)
    
    if p_len == t_len:  # Check for substitution
        for i in range(p_len):
            if text[i] != pattern[i]:
                edits += 1
                if edits > 1:
                    return (2, None)
        return (edits, 'S')
        
    elif p_len == t_len + 1:  # Deletion
        i = j = 0
        while i < t_len and j < p_len:
            if text[i] != pattern[j]:
                edits += 1
                j += 1
                if edits > 1:
                    return (2, None)
            else:
                i += 1
                j += 1
        return (1, 'D')
        
    else:  # Insertion
        i = j = 0
        while i < t_len and j < p_len:
            if text[i] != pattern[j]:
                edits += 1
                i += 1
                if edits > 1:
                    return (2, None)
            else:
                i += 1
                j += 1
        return (1, 'I')

def find_mismatched_matches(text, pattern, exact_matches):
    """Find approximate matches with optimized distance calculation."""
    matches = []
    matched_indices = set()
    
    # Skip exact match regions
    for match in exact_matches:
        for i in range(match, match + len(pattern)):
            matched_indices.add(i)
    
    i = 0
    while i < len(text) - len(pattern) + 1:
        # Skip if in matched region
        if any(idx in matched_indices for idx in range(i, i + len(pattern))):
            i += 1
            continue
            
        # Try each window size
        for length in range(len(pattern) - 1, len(pattern) + 2):
            if i + length <= len(text):
                substring = text[i:i + length]
                dist, mtype = quick_distance(substring, pattern)
                
                if dist == 1:
                    matches.append((i, mtype))
                    break
        i += 1
    
    return matches

def read_fasta(file_path):
    """Read and yield sequences from FASTA file, supporting plain FASTA, BGZIP, and ZIP formats."""
    print(f"Attempting to read: {file_path}")
    
    if file_path.endswith(".zip"):
        with zipfile.ZipFile(file_path, 'r') as zipf:
            print(f"Files in zip: {zipf.namelist()}")
            for filename in zipf.namelist():
                if 'GCF' in filename and filename.endswith('.fna'):
                    print(f"Processing FASTA file: {filename}")
                    with zipf.open(filename) as f:
                        f = io.TextIOWrapper(f)
                        for record in SeqIO.parse(f, "fasta"):
                            seq_str = str(record.seq).upper()
                            print(f"\nSequence verification for {record.id}:")
                            print(f"First 50 bases: {seq_str[:50]}")
                            print(f"Sequence length: {len(seq_str)}")
                            yield record.id, seq_str
    
    elif file_path.endswith(".bgz") or file_path.endswith(".gz"):
        with gzip.open(file_path, 'rt') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_str = str(record.seq).upper()
                print(f"\nSequence verification for {record.id}:")
                print(f"First 50 bases: {seq_str[:50]}")
                print(f"Sequence length: {len(seq_str)}")
                yield record.id, seq_str

    else:  # Handle plain FASTA files
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_str = str(record.seq).upper()
                print(f"\nSequence verification for {record.id}:")
                print(f"First 50 bases: {seq_str[:50]}")
                print(f"Sequence length: {len(seq_str)}")
                yield record.id, seq_str


def process_exact_matches(args):
    """Process a single chromosome for exact matches only."""
    chrom_id, sequence, pattern = args
    start_time = time()
    
    with print_lock:
        print(f"\nProcessing exact matches for {chrom_id}:")
        print(f"Sequence length: {len(sequence):,} bases")
    
    exact_matches = kmp_search(sequence, pattern)
    
    processing_time = time() - start_time
    bases_per_second = len(sequence) / processing_time
    
    current_time = time()
    with progress_lock:
        processed_chromosomes.value += 1
        total_bases_processed.value += len(sequence)
        completion = (processed_chromosomes.value / total_chromosomes.value) * 100
        total_gb = total_bases_processed.value / 1_000_000_000
        
        if current_time - last_update_time.value >= 10:
            print(f"\nExact Match Progress ({datetime.now()}):")
            print(f"Processed {processed_chromosomes.value}/{total_chromosomes.value} chromosomes ({completion:.2f}%)")
            print(f"Total bases processed: {total_gb:.2f} Gb")
            print(f"Current speed: {bases_per_second:,.0f} bases/second")
            if exact_matches:
                print(f"Found {len(exact_matches)} exact matches in {chrom_id}")
            last_update_time.value = current_time
    
    return chrom_id, exact_matches

def process_mismatched_matches(args):
    """Process a single chromosome for mismatched matches using exact match results."""
    chrom_id, sequence, pattern, exact_matches = args
    start_time = time()
    
    with print_lock:
        print(f"\nProcessing mismatched matches for {chrom_id}:")
        print(f"Sequence length: {len(sequence):,} bases")
    
    mismatched_matches = find_mismatched_matches(sequence, pattern, exact_matches)
    
    processing_time = time() - start_time
    bases_per_second = len(sequence) / processing_time
    
    current_time = time()
    with progress_lock:
        processed_chromosomes.value += 1
        total_bases_processed.value += len(sequence)
        completion = (processed_chromosomes.value / total_chromosomes.value) * 100
        total_gb = total_bases_processed.value / 1_000_000_000
        
        if current_time - last_update_time.value >= 10:
            print(f"\nMismatch Progress ({datetime.now()}):")
            print(f"Processed {processed_chromosomes.value}/{total_chromosomes.value} chromosomes ({completion:.2f}%)")
            print(f"Total bases processed: {total_gb:.2f} Gb")
            print(f"Current speed: {bases_per_second:,.0f} bases/second")
            if mismatched_matches:
                print(f"Found {len(mismatched_matches)} mismatched matches in {chrom_id}")
            last_update_time.value = current_time
    
    return chrom_id, mismatched_matches

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 script.py genome.zip pattern.fna num_cores")
        sys.exit(1)
    
    overall_start_time = time()
    genome_file = sys.argv[1]
    pattern_file = sys.argv[2]
    num_cores = int(sys.argv[3])
    
    # Read pattern
    with open(pattern_file, "r") as f:
        pattern = ''.join(line.strip() for line in f if not line.startswith('>'))
    print("\nPattern verification:")
    print(f"Pattern length: {len(pattern)}")
    print(f"Pattern: {pattern}")
    
    # Read all chromosomes
    chromosomes = []
    for chrom_id, sequence in read_fasta(genome_file):
        chromosomes.append((chrom_id, sequence, pattern))
    
    total_chromosomes.value = len(chromosomes)
    print(f"\nFound {len(chromosomes)} chromosomes")
    
    # Phase 1: Process exact matches
    print("\nPhase 1: Finding exact matches...")
    exact_start_time = time()
    processed_chromosomes.value = 0
    total_bases_processed.value = 0
    last_update_time.value = time()
    
    exact_results = []
    with multiprocessing.Pool(num_cores) as pool:
        for result in pool.imap_unordered(process_exact_matches, chromosomes, chunksize=1):
            exact_results.append(result)
    
    # Organize exact match results
    exact_matches_by_chrom = {chrom_id: matches for chrom_id, matches in exact_results}
    total_exact = sum(len(matches) for matches in exact_matches_by_chrom.values())
    
    exact_time = time() - exact_start_time
    print(f"\nPhase 1 Complete:")
    print(f"Total exact matches found: {total_exact}")
    print(f"Time taken: {exact_time:.2f} seconds")
    
    # Phase 2: Process mismatched matches
    print("\nPhase 2: Finding mismatched matches...")
    mismatch_start_time = time()
    processed_chromosomes.value = 0
    total_bases_processed.value = 0
    last_update_time.value = time()
    
    # Prepare data for mismatched matching
    mismatch_inputs = [
        (chrom_id, sequence, pattern, exact_matches_by_chrom[chrom_id])
        for chrom_id, sequence, pattern in chromosomes
    ]
    
    mismatched_results = []
    with multiprocessing.Pool(num_cores) as pool:
        for result in pool.imap_unordered(process_mismatched_matches, mismatch_inputs, chunksize=1):
            mismatched_results.append(result)
    
    # Organize mismatched results
    mismatched_matches_by_chrom = {chrom_id: matches for chrom_id, matches in mismatched_results}
    total_mismatched = sum(len(matches) for matches in mismatched_matches_by_chrom.values())
    
    mismatch_time = time() - mismatch_start_time
    total_time = time() - overall_start_time
    total_gb = total_bases_processed.value / 1_000_000_000
    
    print(f"\nPhase 2 Complete:")
    print(f"Total mismatched matches found: {total_mismatched}")
    print(f"Time taken: {mismatch_time:.2f} seconds")
    
    print("\nOverall Results:")
    print(f"Total exact matches: {total_exact}")
    print(f"Total mismatched matches: {total_mismatched}")
    print(f"Total runtime: {total_time:.2f} seconds")
    
    # Write detailed results
    with open("all_matches.txt", "w") as f:
        f.write(f"Pattern used: {pattern}\n\n")
        f.write(f"Phase 1 (Exact matches) time: {exact_time:.2f} seconds\n")
        f.write(f"Phase 2 (Mismatched matches) time: {mismatch_time:.2f} seconds\n")
        f.write(f"Total runtime: {total_time:.2f} seconds\n\n")
        
        f.write("Exact matches by chromosome:\n")
        for chrom_id, matches in exact_matches_by_chrom.items():
            if matches:
                f.write(f"{chrom_id}: {len(matches)} matches at positions {matches}\n")
        
        f.write("\nMismatched matches by chromosome:\n")
        for chrom_id, matches in mismatched_matches_by_chrom.items():
            if matches:
                f.write(f"{chrom_id}: {len(matches)} matches with details:\n")
                for pos, mtype in matches:
                    f.write(f"  Position {pos}: {mtype} (S=substitution, I=insertion, D=deletion)\n")

if __name__ == "__main__":
    main()
