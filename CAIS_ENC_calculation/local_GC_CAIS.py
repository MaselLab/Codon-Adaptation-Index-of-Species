"""

@author: Andrew
"""

import os, sys, json, csv, datetime,math


# Verbose prints out the progress of the script for the user when set to True
Verbose = False

maxInt = sys.maxsize

while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

Total_AA_freqTable = {'F':0.03620088531983065,
                'L':0.09712159116035492,
                'S':0.08444485585103126,
                'Y':0.027075491234830648,
                '*':0.0015910750475627488,
                'C':0.022000259378660177,
                'W':0.011688204745587587,
                'P':0.06074804607959609,
                'H':0.026040014458008767,
                'Q':0.04793455554691821,
                'R':0.05583288289673258,
                'I':0.04446309540942312,
                'M':0.021993612231572347,
                'T':0.05393328314898733,
                'N':0.03749414829125922,
                'K':0.058767417665229714,
                'V':0.06066089459555564,
                'A':0.06762234528582944,
                'D':0.049382150317773404,
                'E':0.07124589639253541,
                'G':0.06375929494272072}

def calc_CAIS(window):
	RawCount = {'F':{'TTT':0,'TTC':0},
	            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
	            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
	            'Y':{'TAT':0,'TAC':0},
	            '*':{'TAA':0,'TAG':0,'TGA':0},
	            'C':{'TGT':0,'TGC':0},
	            'W':{'TGG':0},
	            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
	            'H':{'CAT':0,'CAC':0},
	            'Q':{'CAA':0,'CAG':0},
	            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
	            'I':{'ATT':0,'ATC':0,'ATA':0},
	            'M':{'ATG':0},
	            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
	            'N':{'AAT':0,'AAC':0},
	            'K':{'AAA':0,'AAG':0},
	            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
	            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
	            'D':{'GAT':0,'GAC':0},
	            'E':{'GAA':0,'GAG':0},
	            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

	Sum = {'F':0,
	       'L':0,
	       'S':0,
	       'Y':0,
	       '*':0,
	       'C':0,
	       'W':0,
	       'P':0,
	       'H':0,
	       'Q':0,
	       'R':0,
	       'I':0,
	       'M':0,
	       'T':0,
	       'N':0,
	       'K':0,
	       'V':0,
	       'A':0,
	       'D':0,
	       'E':0,
	       'G':0}

	ProbTable = {'F':{'TTT':0,'TTC':0},
	            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
	            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
	            'Y':{'TAT':0,'TAC':0},
	            '*':{'TAA':0,'TAG':0,'TGA':0},
	            'C':{'TGT':0,'TGC':0},
	            'W':{'TGG':0},
	            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
	            'H':{'CAT':0,'CAC':0},
	            'Q':{'CAA':0,'CAG':0},
	            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
	            'I':{'ATT':0,'ATC':0,'ATA':0},
	            'M':{'ATG':0},
	            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
	            'N':{'AAT':0,'AAC':0},
	            'K':{'AAA':0,'AAG':0},
	            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
	            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
	            'D':{'GAT':0,'GAC':0},
	            'E':{'GAA':0,'GAG':0},
	            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

	Expected = {'F':{'TTT':0,'TTC':0},
	            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
	            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
	            'Y':{'TAT':0,'TAC':0},
	            '*':{'TAA':0,'TAG':0,'TGA':0},
	            'C':{'TGT':0,'TGC':0},
	            'W':{'TGG':0},
	            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
	            'H':{'CAT':0,'CAC':0},
	            'Q':{'CAA':0,'CAG':0},
	            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
	            'I':{'ATT':0,'ATC':0,'ATA':0},
	            'M':{'ATG':0},
	            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
	            'N':{'AAT':0,'AAC':0},
	            'K':{'AAA':0,'AAG':0},
	            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
	            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
	            'D':{'GAT':0,'GAC':0},
	            'E':{'GAA':0,'GAG':0},
	            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

	Observed = {'F':{'TTT':0,'TTC':0},
	            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
	            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
	            'Y':{'TAT':0,'TAC':0},
	            '*':{'TAA':0,'TAG':0,'TGA':0},
	            'C':{'TGT':0,'TGC':0},
	            'W':{'TGG':0},
	            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
	            'H':{'CAT':0,'CAC':0},
	            'Q':{'CAA':0,'CAG':0},
	            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
	            'I':{'ATT':0,'ATC':0,'ATA':0},
	            'M':{'ATG':0},
	            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
	            'N':{'AAT':0,'AAC':0},
	            'K':{'AAA':0,'AAG':0},
	            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
	            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
	            'D':{'GAT':0,'GAC':0},
	            'E':{'GAA':0,'GAG':0},
	            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}


	#define parameters for the current window
	if window != []:
		window_UID = window[0]
		Species_UID = int(window[1])
		gene_clump_nucleotide_sequence = window[3]
		sequenceLength = len(gene_clump_nucleotide_sequence)
		GC_local_prob = float(window[2])/100
		if GC_local_prob == 0:
			GC_local_prob = 1/100
		notGC_local_prob = 1- GC_local_prob

	# We partition each coding sequence into a list codons and count the codons

	codonList = [gene_clump_nucleotide_sequence[n:n+3] for n in range(0,sequenceLength,3)]

	for AA in RawCount:
		for Codon in RawCount[AA]:
			RawCount[AA][Codon] += codonList.count(Codon)

	for AA in RawCount:
		Sum[AA] = sum(RawCount[AA].values())
		for Codon in RawCount[AA]:
			if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
				Prob = (notGC_local_prob*notGC_local_prob*notGC_local_prob)*0.125
			elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
				Prob = (GC_local_prob*GC_local_prob*GC_local_prob)*0.125
			elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC','TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
				Prob = (GC_local_prob*notGC_local_prob*notGC_local_prob)*0.125
			elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT','AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
				Prob = (GC_local_prob*GC_local_prob*notGC_local_prob)*0.125
			else:
				print("storm the castle")
			ProbTable[AA][Codon] = Prob

	for AA in ProbTable:
		for Codon in ProbTable[AA]:
			Expected[AA][Codon] = ProbTable[AA][Codon]/sum(ProbTable[AA].values())

	#calculate observed codon frequencies
	for AA in RawCount:
		for Codon in RawCount[AA]:
			if Sum[AA] == 0:
				Observed[AA][Codon] = 0
			else:
				Observed[AA][Codon] = RawCount[AA][Codon]/Sum[AA]

	#CAIS
	CAIS=0
	for AA in RawCount:
		for Codon in RawCount[AA]:
			if Observed[AA][Codon] == 0:
				CAIS+=0
			else:
				CAIS += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])

	#weighted_CAIS
	weighted_CAIS = 0
	for AA in RawCount:
		if Sum[AA] == 0:
			aa_weight = 1
		else:
			aa_weight = Total_AA_freqTable[AA]/sum(Observed[AA].values())
		for Codon in RawCount[AA]:
			if Observed[AA][Codon] == 0:
				CAIS+=0
			else:
				weighted_CAIS += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])*aa_weight

	return [CAIS,weighted_CAIS,len(codonList)]


#Geometric mean of local CAIS values
def main(argv):
	log_CAIS_sum = 0
	log_weighted_CAIS_sum = 0
	length_sum = 0
	with open(sys.argv[1], newline='') as csvfile:
		reader = csv.reader(csvfile, delimiter = ",")
		next(reader)
		for row in reader:
			winstats = calc_CAIS(row)
			current_CAIS = winstats[0]
			current_CAIS_weighted = winstats[1]
			current_length = winstats[2]
			if current_CAIS == 0:
				log_CAIS_sum += 0
				log_weighted_CAIS_sum += 0
				length_sum += 0
			else:
				log_CAIS_sum += math.log(current_CAIS)
				log_weighted_CAIS_sum += math.log(current_CAIS_weighted)
				length_sum += current_length

	CAIS = math.exp(log_CAIS_sum/length_sum)
	CAIS_weighted = math.exp(log_weighted_CAIS_sum/length_sum)
	Species = sys.argv[1].split("_")[0]
	print(Species + "," + str(CAIS) + "," + str(CAIS_weighted))


if __name__ == "__main__":
	main(sys.argv[1:])