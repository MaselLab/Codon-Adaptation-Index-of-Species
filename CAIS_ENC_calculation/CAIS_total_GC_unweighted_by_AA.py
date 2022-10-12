# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 19:15:16 2020

@author: Catarina
"""

import os, sys, json, csv, mysql.connector,datetime,math

"""
The purpose of this file is to determine the codon actually adaptation indices (CAAI) for individual species.CAAI is Codon Adaptation Index (1987 version) contolled for Amino Acid composition and potential mutation bias within codon usuage.
This is a metric that describes codon usage patterns. 
CAAI is calculated in the following way:
    
    I)
  
        1) Calculate the probability that a given nucleotide is a G or a C.
            -done easily, by counting the number of Gs and Cs in each codon, then counting the total number of each codon in a species' proteome
            -p(G or C) = total number of G or C in proteome / (3* number of codons in proteome)
        
        2) Calculate probability of a codon given p(G or C)), i.e.
            -p(GCA|p(G or C))) = (p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2)
    
        3) Calculate probability of codon given occurence of its amino acid, i.e.
            -p(GCA|A)=p(GCA| GCA or GCT or GCC or GCG) = ((p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2))/((p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2) + (p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2) +
                                                           (p(G or C)/2)*((p(G or C)/2)*(p(G or C)/2) + (p(G or C)/2)*((p(G or C)/2)*(p(G or C)/2))
            -this is the null expectation
    
        4) Calculate the total number of each codon
        
        5) Calculate the frequency of a codon given the occurence of its amino acid
            -# of times GCA occurs / # of times GCA and synonymous codons occur i.e.
            -f(GCA|A)
    
        6) Calculate the relative adaptiveness of each codon i , w_i
            - w_i = f(i|AA)/p(i|AA)
        
    II) Calculate the Amino Acid controlled codon frequency
        
        7)                                                   
        
"""


####################################################################################
#                                 User Information                                 #
####################################################################################

# SQL Connection Information
Database = ''
User = ''
Host = ''
Password =''

# The tables where the coding sequences are stored
# ** These are stored as a list so that we can iterate through them **
CodingTables = [ 'NCBIGenomes_Coding_Complete']
#,'Backup_EnsemblGenomes_Coding_Complete_05062019'
# Verbose prints out the progress of the script for the user when set to True
Verbose = False


####################################################################################
#                             Program Executes Below                               #
####################################################################################

# A connection to the SQL server is established
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

#########For Total Dataset metrics
                    
# We'll pull species from all coding data tables one at a time
for CodingTable in CodingTables:

    # We find the min/max value of the species UIDs in the table and then
    # iterate through all integers in that range so we only deal with one
    # species at a time-Sara W

    SelectMaxSpeciesUID = "SELECT MAX(SpeciesUID) FROM %s"%CodingTable
    mycursor.execute(SelectMaxSpeciesUID)
    MaxSpeciesUID = mycursor.fetchone()[0]

    SelectMinSpeciesUID = "SELECT MIN(SpeciesUID) FROM %s"%CodingTable
    mycursor.execute(SelectMinSpeciesUID)
    MinSpeciesUID = mycursor.fetchone()[0]
    
    #The values in this table are calculated in another script (WHICH TAKES & HOURS TO RUN)
    #USED FOR THE RE_WEIGHTING OF CAI

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
    
    Species_Total_GC_content = {4:0.4105900,
                                5:0.3984647,
                                461:0.4053087,
                                463:0.4118122,
                                407:0.4134434,
                                408:0.4168525,
                                467:0.3942816,
                                409:0.4250931,
                                410:0.4162005,
                                411:0.4182618,
                                468:0.4231701,
                                469:0.4184607,
                                9:0.4066691,
                                470:0.4194617,
                                471:0.4119122,
                                412:0.4136108,
                                413:0.4134932,
                                414:0.4111889,
                                472:0.4097088,
                                12:0.4072864,
                                415:0.4007741,
                                14:0.4175329,
                                15:0.4116991,
                                416:0.4147482,
                                16:0.4138114,
                                474:0.4182248,
                                475:0.4319899,
                                18:0.4120919,
                                477:0.4327206,
                                417:0.4070456,
                                478:0.4432461,
                                20:0.4136309,
                                480:0.4096289,
                                418:0.4163353,
                                481:0.4085882,
                                482:0.4107121,
                                22:0.4200345,
                                483:0.4340620,
                                484:0.4182700,
                                485:0.3830205,
                                486:0.3702761,
                                23:0.3641617,
                                25:0.4276346,
                                420:0.4111542,
                                421:0.4332301,
                                28:0.4151455,
                                423:0.4094320,
                                424:0.4209174,
                                488:0.4266520,
                                489:0.4210899,
                                491:0.4161310,
                                29:0.4188767,
                                31:0.4096652,
                                492:0.4094208,
                                493:0.4103237,
                                496:0.4533417,
                                497:0.4167646,
                                35:0.4168004,
                                498:0.4061100,
                                499:0.4148764,
                                36:0.4148158,
                                425:0.4172848,
                                38:0.4098491,
                                502:0.4003387,
                                40:0.4260456,
                                503:0.4102434,
                                41:0.4090246,
                                506:0.4111120,
                                42:0.3938096,
                                427:0.4170514,
                                43:0.4106600,
                                45:0.4120511,
                                46:0.4125298,
                                47:0.4110758,
                                429:0.4067810,
                                430:0.4274666,
                                49:0.4303131,
                                50:0.4112322,
                                51:0.4291126,
                                431:0.4170559,
                                53:0.4256563,
                                54:0.4215828,
                                55:0.4267324,
                                57:0.4213925,
                                433:0.4220165,
                                434:0.4239676,
                                58:0.4231795,
                                59:0.4188740,
                                60:0.4095077,
                                61:0.4273248,
                                437:0.4200084,
                                438:0.4176694,
                                439:0.4226807,
                                62:0.4563402,
                                63:0.4419925,
                                65:0.4151206,
                                67:0.4125521,
                                68:0.4099843,
                                69:0.4178094,
                                442:0.4215060,
                                72:0.4338484,
                                73:0.4293511,
                                445:0.4186339,
                                446:0.4091846,
                                75:0.4225226,
                                447:0.3966267,
                                76:0.4317794,
                                77:0.4210380,
                                450:0.4009270,
                                80:0.4154242,
                                452:0.4396255,
                                82:0.4186202,
                                83:0.4129064,
                                84:0.4547498,
                                453:0.4055075,
                                455:0.4194594,
                                86:0.3973327,
                                87:0.3914836}

############################ACTUAL CAI CALCULATING CODE#####################################
    
    for i in range(406,507):
        if i in Species_Total_GC_content.keys():
            if Verbose == True:
                print('Extracting SpeciesUID: %s'%i)
                sys.stdout.flush()
        
            totalCodonCount = 0
            totalnucleotidecount=0
            totalGCcount=0

        # We'll keep dictionaries of all the codons that we count, their RSCU values,
        # and their relative adaptedness values, sorted by which amino acid they
        # correspond to. This will make calculating the CAI relatively easy and painless- Sara W

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
        
            Codon_freqTable_raw = {'F':{'TTT':0,'TTC':0},
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

            Codon_Summed_freqTable = {'F':0,
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
            
            Sum_ProbTable= {'F':0,
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
        
            relative_Codon_freqTable = {'F':{'TTT':0,'TTC':0},
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
        
            Codon_freqTable_weighted = {'F':{'TTT':0,'TTC':0},
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

            RelativeAdaptednessTable = {'F':{'TTT':0,'TTC':0},
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
            
            RelativeAdaptednessTable_weighted = {'F':{'TTT':0,'TTC':0},
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

######################################PUTTING STUFF INTO THE DICTIONARIES

        # We start by attempting to extract a set of coding sequences sharing the specified species UID
            SelectionStatement = "SELECT CodingSequence FROM %s WHERE SpeciesUID = %s"%(CodingTable,i)
            mycursor.execute(SelectionStatement)
            speciesSpecificResults = mycursor.fetchall()

        # If the species UID doesn't exist in the table, we move on, otherwise we iterate through all
        # extracted coding sequences to count the occurences of each codon

            if speciesSpecificResults != []:
                
                for CodingSequence in speciesSpecificResults:
                    CodingSequence = CodingSequence[0]
                    sequenceLength = len(CodingSequence)

                # We partition each coding sequence into a list codons and then compare them to our dictionary
                    codonList = [CodingSequence[n:n+3] for n in range(0,sequenceLength,3)]

                # There are some "codons" that we will be ignoring. Specifically, anything that has an x or some other
                # irregular character in it
                    for AA in RawCount:
                        for Codon in RawCount[AA]:
                            CodonCount = codonList.count(Codon)
                            RawCount[AA][Codon] += CodonCount
                            totalCodonCount += CodonCount
                            totalnucleotidecount += CodonCount*3
                            
                            if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                                totalGCcount +=0
                            elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                                totalGCcount += CodonCount*3
                            elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
                                totalGCcount += CodonCount
                            elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                                totalGCcount += CodonCount
                            elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
                                totalGCcount +=CodonCount*2
                            elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                                totalGCcount +=CodonCount*2
                            
                GC_genic_prob = totalGCcount/totalnucleotidecount
                notGC_genic_prob = 1 - GC_genic_prob
                
                GC_total_prob = float(Species_Total_GC_content[i])
                notGC_total_prob = 1- GC_genic_prob
            #print(GCprob)

            # Once the raw counts are found for each codon, we then start calculating RSCU_i values (see the beginning of
            # this script for the specifics).
                Summed_codon_instances = 0 
 
                for AA in RawCount:
                    Summed_codon_instances += sum(RawCount[AA].values())
                    Sum[AA] = sum(RawCount[AA].values())
                    for Codon in RawCount[AA]:

                    # In some weird cases an amino acid never shows up in a particular genome... if this is the case, we set
                    # the RSCU_i value to one for all codons corresponding to that amino acid. We do this because later we'll
                    # take the geometric mean of all relative adaptedness values and multiplying by one will not affect this value
                    # In essence, setting them to one "silences" these.
                        if Sum[AA] != 0:
                            if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                                Prob = (notGC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
                            elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                                Prob = (GC_total_prob*GC_total_prob*GC_total_prob)*0.125
                            elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
                                Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
                            elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                                Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
                            elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
                                Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)*0.125
                            elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                                Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)*0.125
                            else:
                                print("storm the castle")
                            ProbTable[AA][Codon] = Prob
                            
                            
            if Verbose == True: 
                print("-------SumTable-----------------")    
                print(Sum)      
                print("-------ProbTable-----------------")    
                print(ProbTable)   
                print("-------SUMProbTable-----------------")    
                print(Sum_ProbTable)
                sys.stdout.flush()
                            
            ####################UPDATE TABLE WITH RELATIVE PROBABILITIES###########################################
                            
            for AA in ProbTable:
                Sum_ProbTable[AA] = sum(ProbTable[AA].values())
                for Codon in ProbTable[AA]:
                    if ProbTable[AA] != 0:
                        Prob = ProbTable[AA][Codon]/float(Sum_ProbTable[AA])
                    else:
                        print("storm the castle")
                    ProbTable[AA][Codon] = Prob
            #print(Sum)
            #################################calculate codon frequencies############################
            
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    Codon_freq_raw = RawCount[AA][Codon]/Summed_codon_instances
                    Codon_freqTable_raw[AA][Codon] = Codon_freq_raw
                    
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    relative_Codon_freq = RawCount[AA][Codon]/Sum[AA]
                    relative_Codon_freqTable[AA][Codon] = relative_Codon_freq
                    
            # Sum of codon frequencies per AA
            for AA in Codon_freqTable_raw:
                Codon_freq_sum = 0
                for Codon in Codon_freqTable_raw[AA]:
                    Codon_freq_sum += Codon_freqTable_raw[AA][Codon]
                Codon_Summed_freqTable[AA]= Codon_freq_sum
                
            
            for AA in Codon_freqTable_raw:
                codon_count = len(RawCount[AA])       
                for Codon in Codon_freqTable_raw[AA]:
                    if codon_count*Codon_freqTable_raw[AA][Codon]!=0:
                        codon_weight = Total_AA_freqTable[AA]/(Codon_Summed_freqTable[AA])
                        Codon_freqTable_weighted[AA][Codon] = Codon_freqTable_raw[AA][Codon]*codon_weight
                    else:
                        codon_weight = 0
                        Codon_freqTable_weighted[AA][Codon] = Codon_freqTable_raw[AA][Codon]*codon_weight
                        
                        
             ######## Once all Prob_i values are found and relative codon frequencies, we can calculate the relative adaptedness values and save those to our dictionary
            for AA in ProbTable:
                for Codon in ProbTable[AA]:
                    RelativeAdaptedness = relative_Codon_freqTable[AA][Codon]/ProbTable[AA][Codon]
                    RelativeAdaptednessTable[AA][Codon] = RelativeAdaptedness
    

            for AA in ProbTable:
                MaxProb = max(ProbTable[AA].values())
                for Codon in ProbTable[AA]:
                    RelativeAdaptedness = relative_Codon_freqTable[AA][Codon]/ProbTable[AA][Codon]
                    NormalizedRelativeAdaptedness = RelativeAdaptedness/MaxProb
                    RelativeAdaptednessTable_weighted[AA][Codon] = NormalizedRelativeAdaptedness
            
            if Verbose == True: 
                print("-------Codon_freqTable_raw-----------------")
                print(Codon_freqTable_raw)
                print("---------relative_Codon_freqTable---------------")    
                print(relative_Codon_freqTable)
                print("--------ProbTable----------------")    
                print(ProbTable)   
                print("-----------RelativeAdaptednessTable-------------")    
                print(RelativeAdaptednessTable)           
                sys.stdout.flush()

            #CAI
            LogOfCAAI = 0
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    # Here, we're just finding the number of occurrences of each relative adaptedness value, taking the log of w_i, and
                    # adding it to the total Log of CAI
                    for k in range(0,RawCount[AA][Codon]):
                        LogOfCAAI += math.log(RelativeAdaptednessTable[AA][Codon])      
            # We divide by the total number of codons in all coding sequences
            LogOfCAAI = (1/totalCodonCount)*LogOfCAAI
            # and we invert the log to get the CAI
            CAAI = math.exp(LogOfCAAI)            
            

            #Weighted wi CAI
            Weighted_wi_Log_ofCAAI = 0
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    #for k in range(0,RawCount[AA][Codon]):
                    if Verbose == True:
                        print(Codon)
                        sys.stdout.flush()

                    if RelativeAdaptednessTable[AA][Codon] != 0:
                        Weighted_wi_Log_ofCAAI += Codon_freqTable_raw[AA][Codon]*math.log(RelativeAdaptednessTable_weighted[AA][Codon])

                    else:
                        if Verbose == True:
                            print("%s,%s"%(RawCount[AA][Codon],RelativeAdaptednessTable_weighted[AA][Codon]))
                            sys.stdout.flush()
                        else:
                            pass

                    if Verbose == True:
                        print(Weighted_wi_Log_ofCAAI)
                        sys.stdout.flush()
  
            # and we invert the log to get the CAI
            Weighted_wi_CAAI = math.exp(Weighted_wi_Log_ofCAAI)
            if Verbose == True:
                    print(Weighted_wi_Log_ofCAAI)
                    sys.stdout.flush()




           #Weighted codon frequency CAI
            Weighted_fi_Log_ofCAAI = 0
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    #for k in range(0,RawCount[AA][Codon]):
                    if Verbose == True:
                        print(Codon)
                        sys.stdout.flush()

                    if RelativeAdaptednessTable[AA][Codon] != 0:
                        Weighted_fi_Log_ofCAAI += Codon_freqTable_weighted[AA][Codon]*math.log(RelativeAdaptednessTable[AA][Codon])

                    else:
                        if Verbose == True:
                            print("%s,%s"%(RawCount[AA][Codon],RelativeAdaptednessTable[AA][Codon]))
                            sys.stdout.flush()
                        else:
                            pass

                    if Verbose == True:
                        print(Weighted_fi_Log_ofCAAI)
                        sys.stdout.flush()
  
            # and we invert the log to get the CAI
            Weighted_fi_CAAI = math.exp(Weighted_fi_Log_ofCAAI)
            if Verbose == True:
                    print(Weighted_fi_Log_ofCAAI)
                    sys.stdout.flush()

            #print("%s,%s,%s,%s"%(i,CAAI,Weighted_wi_CAAI,Weighted_fi_CAAI))
            print("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s"%(i,Weighted_fi_CAAI,Codon_Summed_freqTable['F'],Codon_Summed_freqTable['L'],Codon_Summed_freqTable['S'],Codon_Summed_freqTable['Y'],Codon_Summed_freqTable['*'],Codon_Summed_freqTable['C'],Codon_Summed_freqTable['W'],Codon_Summed_freqTable['P'],Codon_Summed_freqTable['H'],Codon_Summed_freqTable['Q'],Codon_Summed_freqTable['R'],Codon_Summed_freqTable['I'],Codon_Summed_freqTable['M'],Codon_Summed_freqTable['T'],Codon_Summed_freqTable['N'],Codon_Summed_freqTable['K'],Codon_Summed_freqTable['V'],Codon_Summed_freqTable['A'],Codon_Summed_freqTable['D'],Codon_Summed_freqTable['E'],Codon_Summed_freqTable['G']))

cnx.close()

                        
                        
             
