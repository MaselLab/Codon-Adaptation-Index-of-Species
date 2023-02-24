import csv, sys, math
from pathlib import Path
import pandas as pd

'''

Author: Andrew Wheeler
Date Updated: 10/12/2022

script to calculate RSCUS values for each local region of a genome
reads output file from "local_GC_windows.py"
Outputs a csv file with the 64 RSCUS values and codon counts for each local window

Use: python RSCUs_intergenic.py [output_from_local_GC_windows.csv]

'''

maxInt = sys.maxsize

while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def calculate_RSCUS(window):

  #initiate empty tables
    
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

  obs_freq_table = {'F':{'TTT':0,'TTC':0},
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

  exp_freq_table = {'F':{'TTT':0,'TTC':0},
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

  local_RSCUS_table = {'F':{'TTT':0,'TTC':0},
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

  # calculate the expected probabilities of each codon given the local GC content 

  for AA in RawCount:
    for Codon in RawCount[AA]:

      if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
          Prob = (notGC_local_prob*notGC_local_prob*notGC_local_prob)
      elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
          Prob = (GC_local_prob*GC_local_prob*GC_local_prob)
      elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
          Prob = (GC_local_prob*notGC_local_prob*notGC_local_prob)
      elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
          Prob = (GC_local_prob*notGC_local_prob*notGC_local_prob)
      elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
          Prob = (GC_local_prob*GC_local_prob*notGC_local_prob)
      elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
          Prob = (GC_local_prob*GC_local_prob*notGC_local_prob)
      else:
          print("storm the castle")
      ProbTable[AA][Codon] = Prob
                        
  # calculate expected probability that a given codon is encoded by a given amino acid

  for AA in RawCount:
    for Codon in RawCount[AA]:
      exp_freq_table[AA][Codon] = ProbTable[AA][Codon]/sum(ProbTable[AA].values())

  # calculate the observed frequency that a given codon is encoded by a given amino acid

  for AA in RawCount:
    for Codon in RawCount[AA]:
      if sum(RawCount[AA].values()) == 0:
        obs_freq_table[AA][Codon] = 1
        exp_freq_table[AA][Codon] = 1
      else:
        obs_freq_table[AA][Codon] = RawCount[AA][Codon]/sum(RawCount[AA].values())

  # calculate the RSCUS score for each codon
                  
  for AA in RawCount:
      for Codon in RawCount[AA]:
          local_RSCUS_table[AA][Codon] = obs_freq_table[AA][Codon]/exp_freq_table[AA][Codon]

  RSCUS_entry = [Species_UID,sequenceLength,GC_local_prob,str(local_RSCUS_table),str(RawCount)]

  return RSCUS_entry

def main(argv):
  RSCUS_df = pd.DataFrame(columns = ['Species_UID','Sequence_length','GC','RSCUS_table','Raw_counts'])
  with open(sys.argv[1], newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter = ",")
    next(reader)
    for row in reader:
      RSCUS_entry = calculate_RSCUS(row)
      RSCUS_df.loc[len(RSCUS_df.index)] = RSCUS_entry

  filename = sys.argv[1].split("_")[0] + "_RSCUS.csv"
  RSCUS_df.to_csv(filename)

if __name__ == "__main__":
  main(sys.argv[1:])
