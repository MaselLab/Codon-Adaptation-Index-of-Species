import pandas as pd  
from pathlib import Path
import csv, sys, gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.Seq import Seq, MutableSeq
from BCBio import GFF

'''

Author: Andrew Wheeler
Date Updated: 10/12/2022

Script to identify local genomic regions and calculate the corresponding GC content for each. Reads in a zipped fasta file of a genome assembly and a zipped gff annotation file for the same assembly.
Identifies flanking sequences around each CDS feature in the genome. Additional CDS regions overlapping these flanking sequences are grouped into the same local region. Calculates the GC content of each local window.
Outputs a csv file containing a UID for each local window, the GC content of each window, the in-frame coding sequence within each window, and an integer species UID given by the user.

Use: python calculate_local_gc.py [zipped_annotation_file.gff.gz] [zipped_genome_assembly_sequence.fna.gz] [species UID]

'''

#Breaks the assembly into scaffolds, sorts all CDS features by which scaffold they are found in, and lists the start and stop coordinates of each CDS feature
#Outputs a dictionary with each scaffold and its associated CDS features, indluding the id, location, and reading frame information.
def parse_input(gff_file,assembly_file):
	scaffold_dict = {}
	gff_handle = gzip.open(gff_file, "rt")

	with gzip.open(assembly_file,"rt") as handle:
		for seq_record in SeqIO.parse(handle, "fasta"):
			scaffold = seq_record.id
			seq = (seq_record.seq)
			scaffold_dict[scaffold] = [seq,[[],[],[]]]

			gff_handle.seek(0)
			limit_info = dict(gff_id=[scaffold],gff_type=["CDS"])

			for rec in GFF.parse(gff_handle,limit_info=limit_info,target_lines=1):
				for feature in rec.features:
					scaffold_dict[scaffold][1][0].append(feature.id)
					scaffold_dict[scaffold][1][1].append(feature.location)
					scaffold_dict[scaffold][1][2].append(int(feature.qualifiers["phase"][0]))
	gff_handle.close()
	handle.close()
	return scaffold_dict

#Masks the coding regions within a sequence
def mask_sequences(sequence,features):
	seq = sequence
	newseq = MutableSeq(str(seq))
	if not features[0]:
		pass
	else:
		for i in range(len(features[0])):
			newseq = GeneMasker(newseq,features[1][i].start,features[1][i].end)
	return newseq

#Masking subfunction
def GeneMasker(sequence,gene_start,gene_end):
	gene_index_list = range(int(gene_start), int(gene_end))
	for index in gene_index_list:
		sequence[index] = "-"
	return sequence

#Identifies the start and stop coordinates of local sequence 'windows' surrounding CDS features and calculates local GC content
#The size of the flanking non-coding regions can be specified, by default it is 1500 base pairs on either side of the features within the window 
#Output a dictionary of each local window in a sequence, its GC content, and the in-frame coding sequence within the window
def create_gene_dict(masked_sequence,sequence,species_UID,features):
	gene_dict = {}
	window_size = 1500
	index = 0
	UID = 0

	#divide the scaffold sequence into local windows around each CDS feature (denoted by masked sequence regions)
	#if local windows overlap, restart the count for the flanking sequence
	while index < len(masked_sequence):
		letter = masked_sequence[index]
		if letter == "-":
			windowstart = index - window_size
			if windowstart < 0:
				windowstart = 0
			while letter == "-" and index < len(masked_sequence):
				letter = masked_sequence[index]
				index += 1
			count = 0
			while count < window_size and index < len(masked_sequence):
				letter = masked_sequence[index]
				if letter == "-":
					count = 0
				else:
					count += 1
				index += 1
			windowend = index

			#calculate intergenic GC content within the window
			window_intergenic_seq = str(masked_sequence[windowstart:windowend]).replace("-","")
			gc = GC(Seq(window_intergenic_seq))

			# identify genes that are within this window
			window_features = gene_assigner(windowstart,windowend,features)
			if not window_features[0]:
				print("no genes in window")

			#merge together features within the same window, trim off codons that are split between windows
			else:
				merging_dict = {}
				for i in range(len(window_features[0])):
					gene_id = window_features[0][i]
					start = window_features[1][i].start
					end = window_features[1][i].end
					strand = window_features[1][i].strand
					phase = window_features[2][i]

					if strand == -1:
						gene_seq = sequence[start:end].reverse_complement()
						if gene_id in merging_dict:
							merging_dict[gene_id] = [gene_seq + merging_dict[gene_id][0],phase]
						else:
							merging_dict[gene_id] = [gene_seq,phase]
					else:
						gene_seq = sequence[start:end]
						if gene_id in merging_dict:
							merging_dict[gene_id][0] = merging_dict[gene_id][0] + gene_seq
						else:
							merging_dict[gene_id] = [gene_seq,phase]

				window_coding_seq = Seq('')
				for gene in merging_dict:
					merging_dict[gene][0] = merging_dict[gene][0][merging_dict[gene][1]:]
					overhang = len(merging_dict[gene][0]) % 3
					if overhang != 0:
						merging_dict[gene][0] = merging_dict[gene][0][:len(merging_dict[gene][0])-overhang]
					window_coding_seq = window_coding_seq + merging_dict[gene][0]

				#create dictionary of the coding sequence within each window, and its local GC content
				gene_dict[UID] = [species_UID,gc,window_coding_seq]
				UID += 1
		index += 1
	return gene_dict

# Assigns features from a genomic feature table to the sequence windows they are found within
def gene_assigner(windowstart,windowend,features):
	window_features = [[],[],[]]
	if not features[1]:
		print("no genes in window")
	else:
		for i in range(len(features[1])):
			if features[1][i].start-1 > windowstart and features[1][i].start-1 < windowend:
				window_features[0].append(features[0][i])
				window_features[1].append(features[1][i])
				window_features[2].append(features[2][i])
	return window_features

#read in the input files, generate local windows for each scaffold in the assembly, write the overall output to a csv file
def main(argv):
	scaffold_dict = parse_input(sys.argv[1],sys.argv[2])
	species_UID = sys.argv[3]
	intergenic_gc_table = pd.DataFrame(columns = ['Species_UID','GC','Sequence'])

	for scaffold in scaffold_dict:
		sequence = scaffold_dict[scaffold][0]
		features = scaffold_dict[scaffold][1]
		masked_sequence = mask_sequences(sequence,features)
		gene_dict = create_gene_dict(masked_sequence,sequence,species_UID,features)

		df = pd.DataFrame.from_dict(data=gene_dict, orient='index',columns = ['Species_UID','GC','Sequence'])	
		intergenic_gc_table = pd.concat([intergenic_gc_table,df],ignore_index = True,sort=False)

	filename = Path(species_UID)
	filename = str(filename.with_suffix('')) + "_local_gc.csv"
	intergenic_gc_table.to_csv(filename)

if __name__ == "__main__":
	main(sys.argv[1:])