#This file was originally created to fraction metagenomes; however, it can fraction most multiline FASTA files

import argparse
#import os
from Bio import SeqIO

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("metagenome_name", help="Name of the input metagenome file (in FASTA format).")
    parser.add_argument("maxsize", type=int, help="Maximum number of sequences per output file.")
    return parser.parse_args()

def writer(filename, seq_list):
    """Write a list of sequences to a FASTA file."""
    with open(filename, "w+") as f:
        SeqIO.write(seq_list, f, "fasta")

def split_metagenome(filehandle, metagenome_name, maxsize):
    """Split the metagenome into smaller FASTA files with a maximum of maxsize sequences."""
    record = SeqIO.parse(filehandle, "fasta") 
    file_counter = 0  # Counter for naming output files
    seq_list = []  
    sequence_counter = 0  # Counter to track the number of sequences in the current file

    for sequence in record:
        if sequence_counter >= maxsize:  
            file_counter += 1 
            output_filename = f"{metagenome_name}fraction{file_counter:03d}.fasta" 
            writer(output_filename, seq_list)  
            seq_list = [sequence]  
            sequence_counter = 1  # Reset the sequence counter
        else:
            seq_list.append(sequence)  # Add the sequence to the current list
            sequence_counter += 1  # Increment the sequence counter

    # Write any remaining sequences to a final output file
    final_filename = f"{metagenome_name}fractionx.fasta"
    writer(final_filename, seq_list)

def main():
    args = parse_args()

    with open(args.metagenome_name, "r") as filehandle:
        split_metagenome(filehandle, args.metagenome_name, args.maxsize)

if __name__ == "__main__":
    main()
