#This program counts 3-mers within input scaffolds/contigs to set up classifying scaffolds (in my use case, as bacteriophage or bacterium-derived) based on 3-mer composition. 
#Ambigious nucleotides will cause the script to terminate, with the exception of "N" as it is present in all scaffolds
#pandas can easily replace sqlite3 if inputs are suitably small

#import itertools
import argparse
from Bio import SeqIO
import sqlite3
#import pandas as pd
#from collections import Counter


#nucleotides = ['A', 'T', 'G', 'C']
# Generate all possible 3-mers using itertools.product
#mers_3 = itertools.product(nucleotides, repeat=3)
#create a dictionary where each 3-mer from the mers_3 iterable is a key and gets a default value of 0
#mers_3_dict = {''.join(m): 0 for m in mers_3}

parser = argparse.ArgumentParser()
#parser.add_argument("-k", "--kmer_length", type=int, default=3, help="Length of the k-mer") #in an alternative version, k-kmer length is provided by the user and the sequence_dict is created accordingly, but here the k-mer length has to be 3 as the sequence_dict is preconfigured. I will have to update this (the above itertools code gets it started).
parser.add_argument("-s", "--sequence", required=True)
parser.add_argument("-f", "--sqlite_filename", required=True)
args = parser.parse_args()

#connect to sqlite database file and initialize cursor
conn = sqlite3.connect(args.sqlite_filename + '.sqlite')
cur = conn.cursor()
cur.execute('DROP TABLE IF EXISTS Counts')
cur.execute('''CREATE TABLE Counts (scaffold_name TEXT, length INTEGER, scaffold_type INTEGER)''')

#define kmer_length as sequence_dict is preconfigured
kmer_length = 3 

for record in SeqIO.parse(args.sequence, 'fasta'):
        sequence_proper = str(record.seq) #get sequence (i.e., nucleotides without header)
        sequence_length = len(sequence_proper)
        scaffold_name = str(record.id) #get scaffold/sequence header/name
        #create dictionary of all possible 3-mers, excluding ambigious nucleotide calls.
        sequence_dict={"ATT":0,"TTG":0,"TGA":0,"GAA":0,"AAA":0,"AAT":0,"ATA":0,"TAA":0,"AAC":0,"ACA":0,"CAG":0,"AGA":0,"TTC":0,"TCT":0,"CTT":0,"GAC":0,"ACC":0,"CCT":0,"TCC":0,"CCC":0,"CCG":0,"CGC":0,"GCG":0,"CGA":0,"ACT":0,"CTG":0,"TGT":0,"GTA":0,"CAC":0,"CAT":0,"ATG":0,"TGG":0,"GGC":0,"GCT":0,"ACG":0,"GCC":0,"CCA":0,"GTG":0,"GGA":0,"GAG":0,"AGC":0,"GCA":0,"TTT":0,"TCA":0,"AGT":0,"GTC":0,"TTA":0,"TAT":0,"ATC":0,"CGT":0,"GGT":0,"TGC":0,"CGG":0,"CTA":0,"TAC":0,"CTC":0,"GAT":0,"TCG":0,"AGG":0,"GTT":0,"AAG":0,"GGG":0,"CAA":0,"TAG":0}
        cur.execute('''INSERT INTO Counts (scaffold_name, length, scaffold_type)
                    VALUES (?,?,?)''', (scaffold_name, sequence_length, 1))
        for i in range(len(sequence_proper)-kmer_length+1):
            sequence_segment = sequence_proper[i:i+kmer_length]
            #this will terminate the script as soon as a non ATGCN character is observed.
            assert all(char in "ATGCN" for char in sequence_segment), "Invalid character in sequence!" 
            #print(list(Counter(sequence_segment))) #this allowed me to examine the composition of each kmer
            #this accounts for Ns in these scaffolds, as no 3-mers in my sequence_dict contain "N" and I don't want such k-mers to be appended to it
            if "N" not in sequence_segment: 
                 sequence_dict[sequence_segment] = sequence_dict.get(sequence_segment, 0) + 1 
        for ky in sequence_dict:
            try:
                cur.execute("ALTER TABLE Counts ADD COLUMN '%s' INTEGER" % ky)
            except:
                pass
            cur.execute('''UPDATE Counts SET '%s' = '%i' WHERE scaffold_name = ?''' % (ky, sequence_dict[ky]), (scaffold_name,))


#commit changes to sqlite database file and close connection to it
conn.commit()
conn.close()
