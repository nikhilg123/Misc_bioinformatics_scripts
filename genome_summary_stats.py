#this script is to perform genome QC on all isolate genomes obtained from a collaborator
#output (currently prints to stdout) includes the name of the genome, the number of sequences that make up the genome, the %GC, N50, the shortest and longest sequence, and the total amount of nucleotides/bp in the genome
#unlike v2, this implements parallelization using concurrent.futures library for a 2-3 fold speed up

from Bio import SeqIO
#from Bio.SeqUtils import gc_fraction #I don't need gc_fraction now as I compute GC manually
import glob
import concurrent.futures
#collections and Counter import not necessary if GC is determined via biopython function "gc_fraction" and similar
from collections import Counter
import time
import os
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("file_path", help="absolute file path to directory containing genome files (in FASTA format).")
args = parser.parse_args()

#docstring and comments could be improved for clarity
def assess_genome_quality(file_path):
	"""quantifies and prints genome stats including genome name, number of sequences, %GC,
	N50, the shortest and longest sequence in the genome, and total bp/nucleotides
	
	Args:
	file_path: an absolute path to the directory containing genomes of interest
	"""
	genome_name = (os.path.basename(file_path)).rstrip(".fasta")
    
	number_of_sequences = 0 
	sequences = []
	for sequence in SeqIO.parse(file_path, 'fasta'):
		number_of_sequences += 1
		sequences.append(str(sequence.seq)) #only the sequence, i.e., no header, is being appended
	sequence_lengths = [len(sequence) for sequence in sequences] #list comprehension to get a list of all sequence lengths
	sequence_lengths.sort(reverse=True) #sorting occurs in place, which allows the shortest and longest sequences to be easily retrieved, and aids the calculation of N50
	total_bp = sum(sequence_lengths)
	nucleotide_counts = Counter("".join(sequences)) #all sequences within genome are concatenated and characters are counted 
	GC = ((nucleotide_counts.get('C') + nucleotide_counts.get('G'))/total_bp)*100 #calculating GC content using output from Counter
	try:
		#this ensures that genomes with ambigious nucleotides (aside from "N") are removed. Other genomes will still run as we are spawning and running multiple processes with concurrent.futures
		assert total_bp == nucleotide_counts.get('C') + nucleotide_counts.get('G') + nucleotide_counts.get('T') + nucleotide_counts.get('A') + nucleotide_counts.get('N'), f"Assertion Error: {genome_name}.fasta contains invalid characters and was ignored"
	except AssertionError as e:
			time.sleep(5)       
			print (e)
			raise
    
#calculating N50 
	counted_bp = 0
	for lengths in sequence_lengths: 
		counted_bp += lengths
		if counted_bp >= (total_bp*0.5):
			N50 = lengths
			break
	#currently output prints to stdout
	#need to revisit print formatting, perhaps implement as dataframe
	#need to revisit formatting and write to an outfile so stdout doesn't need to be redirected to file within shell		
	print(f"genome name:{genome_name}\t number of sequences:{number_of_sequences:,}\t GC:{GC:.2f}%\t N50:{N50:,}\t longest sequence:{sequence_lengths[0]:,}\t shortest sequence:{sequence_lengths[-1]:,}\t total bp:{total_bp:,}")

def main():
	fasta_files = [files for files in glob.glob(f"{args.file_path}/*.fasta") if "GD" in files] #GD is an internal keyword. Modify or remove as needed.
	#print(fasta_files)								 
	with concurrent.futures.ProcessPoolExecutor() as executor: 
		executor.map(assess_genome_quality, fasta_files) #calls the function "assess_genome_quality" for each element in "fasta_files"


# Ensure the script is executed only when run directly, not when imported
if __name__ == '__main__':
    main()