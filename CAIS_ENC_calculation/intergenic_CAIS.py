import csv, sys, json, math
import pandas as pd

'''

Author: Andrew Wheeler
Date Updated: 10/12/2022

Script to calculate local GC CAIS values for a species
Takes the grand geometric mean of RSCUS values for each codon across all local GC windows
Uses exponential of an arithmetic mean for computational efficiency
Uses the output file from intergenic_RSCUS.py
Prints the final CAIS value

Use: python CAIS_intergenic.py [output_from_intergenic_RSCUS.csv]

'''
maxInt = sys.maxsize

while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def main(argv):
    total_codons = 0
    log_RSCUS_sum = 0
    with open(sys.argv[1], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter = ",")
        next(reader)

        #for each local GC window

        for row in reader:
            window_log_RSCUS_sum = 0
            window_codons = 0
            local_RSCUS = json.loads(row[4].replace("'", "\""))
            RawCount = json.loads(row[5].replace("'", "\""))
            for AA in RawCount:

                #for each codon

                for codon in RawCount[AA]:
                    if RawCount[AA][codon] != 0:
                        total_codons += RawCount[AA][codon]

                        #multiply the log(RSCUS) of this codon by the number of times this codon appears in the current window
                        #add to the total sum

                        log_RSCUS_sum += math.log(local_RSCUS[AA][codon])*RawCount[AA][codon]

        #divide the total sum by the total number of codons in the genome

        CAIS = math.exp(log_RSCUS_sum/total_codons)
        
    Species = sys.argv[1].split("_")[0]
    print(Species + " : " + str(CAIS))
    return CAIS

if __name__ == "__main__":
  main(sys.argv[1:])