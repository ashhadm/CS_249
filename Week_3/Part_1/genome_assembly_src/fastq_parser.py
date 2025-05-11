# genome_assembly_src/fastq_parser.py
from Bio import SeqIO
import os

def read_fastq(filename):
    """Read a FASTQ file and return a list of (header, sequence, quality) tuples."""
    records = []
    for record in SeqIO.parse(filename, "fastq"):
        records.append((record.id, str(record.seq), record.letter_annotations["phred_quality"]))
    return records

def read_fasta(filename):
    """Read a FASTA file and return a list of (header, sequence) tuples."""
    records = []
    for record in SeqIO.parse(filename, "fasta"):
        records.append((record.id, str(record.seq)))
    return records

def write_fasta(filename, contigs):
    """Write contigs to a FASTA file."""
    with open(filename, 'w') as f:
        for i, seq in enumerate(contigs):
            f.write(f">contig_{i+1}\n{seq}\n")
    return filename
