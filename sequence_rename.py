#Some bioinformatics tools have challenges dealing with complex headers in fasta filese
#This program allows you to rename genomes/scaffolds/contigs with simple naming schemes, while retaining the original names in a table for later

import argparse
import Bio
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence")
parser.add_argument("--table_outfile")
parser.add_argument("--sequence_outfile")
args = parser.parse_args()

sequence_outfile_handle = open(args.sequence_outfile, 'a')
table_outfile_handle = open(args.table_outfile, 'a')

#need to update with fstring formatting
counter = 1
for record in SeqIO.parse(args.sequence, 'fasta'):
    record_original_header = str(record.id)
    record_new_header = "sequence_%i" %(counter)
    table_outfile_handle.write(record_original_header + "\t" + record_new_header + "\n")
    sequence_outfile_handle.write(">" + record_new_header + "\n" + str(record.seq) +"\n")
    counter += 1
