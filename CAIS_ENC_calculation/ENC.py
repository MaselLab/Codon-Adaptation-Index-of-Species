# -*- coding: utf-8 -*-
"""
Created on Fri May 29 14:57:40 2020

@author: Catarina
"""

import os, sys, json, csv, mysql.connector, datetime,math

"""
The purpose of this file is to determine the effective number of codons (ENC) for individual species.This is a metric that describes codon usage patterns. This script will calculate several versions of ENC
the most fundamental of which is Sweall Wright's ENC from wright et. al, 1990.


i) Wright's N_{c}: 2 + 9/ avg(F_{2}) + 2/ avg(F_{3}) + 5/ avg(F_{4}) + 3/ avg(F_{6}), where avg(F_{i}) is the average F value for amino acids with i synonymous codons and 
    F_{aa}: (N_{aa}*sum from 1 to k of (p_{i}^2 )-1)/(N_{aa} -1), where p_{i} is the frequency of the ith synonymous codon relative to occurences of the amino acid, k is the number of synonymous codons for aa, and N_{aa} is the
    observed number of codons for aa
    
ii) Novembre's N_{c}: 2 + 9/ avg(F_{2}) + 1/ avg(F_{3}) + 5/ avg(F_{4}) + 3/ avg(F_{6}), where avg(F_{i}) is the average F value for amino acids with i synonymous codons and  
    F_{aa}: (X_{aa}^2 + N_{aa} - k)/(k*(N_{aa} -1)), where X_{aa}^2: sum from 1 to k of((N_{aa}(p_{i} - e_{i})^2)/e_{i}), with e_{i} is the expected usage of the ith codon from GC content (in this case)
    
    
Below are the steps in the order written in the script: 
    
    I) This is used for finding e_{i}, N_{aa}, p_{i}:
    
        1) Calculate the total number of each codon
                                                           
        2) Calculate the total number of occurences of each aa
        
        3) Calculate probability of a codon given p(G or C)), i.e.
            -p(GCA|p(G or C))) = (p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2)
    
        4) Calculate probability of codon given occurence of its amino acid, i.e.
            -p(GCA|A)=p(GCA| GCA or GCT or GCC or GCG) = ((p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2))/((p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2) + (p(G or C)/2)*((p(G or C)/2)*((1-p(G or C))/2) +
                                                           (p(G or C)/2)*((p(G or C)/2)*(p(G or C)/2) + (p(G or C)/2)*((p(G or C)/2)*(p(G or C)/2))
            -this is the null expectation (e_{i} or expected usage)
        
        5) Calculate the frequency of a codon given
            a)the occurence of its amino acid
            -# of times GCA occurs / # of times GCA and synonymous codons occur i.e.
            -f(GCA|A): p_{i}
            b) all other codons
    
    II) This is used for finding X_{aa}^2, F_{aa}:
        
                                              
        
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
Verbose = True


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

############################ACTUAL ENC CALCULATING CODE
    
    for i in range(410,411):
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
            
            
            Chi_squared_AA = {'F':0,
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
            
            F_Wright_AA = {'F':0,
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
            
            F_Novembre_AA = {'F':0,
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
                notGC_total_prob = 1- GC_total_prob
            #print(GCprob)

            # EXPECTED OCCURENCE OF GIVEN CODONS (Analagous to RSCU values per codon)
                Summed_codon_instances = 0 
 
                for AA in RawCount:
                    Summed_codon_instances += sum(RawCount[AA].values())
                    Sum[AA] = sum(RawCount[AA].values())
                    for Codon in RawCount[AA]:

                    # In some weird cases an amino acid never shows up in a particular genome... if this is the case, we set
                    # the RSCU_i value to one for all codons corresponding to that amino acid. We do this because later we'll
                    # take the geometric mean of all relative adaptedness values and multiplying by one will not affect this value
                    # In essence, setting them to one "silences" these.
                        if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                            Prob = (notGC_total_prob*notGC_total_prob*notGC_total_prob)/8
                        elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                            Prob = (GC_total_prob*GC_total_prob*GC_total_prob)/8
                        elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
                            Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)/8
                        elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                            Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)/8
                        elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
                            Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)/8
                        elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                            Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)/8
                        else:
                            print(Codon)

                        ProbTable[AA][Codon] = Prob

            #################################calculate codon frequencies
            
                for AA in RawCount:
                    for Codon in RawCount[AA]:
                        Codon_freq_raw = RawCount[AA][Codon]/totalCodonCount
                        Codon_freqTable_raw[AA][Codon] = Codon_freq_raw
                    
                for AA in RawCount:
                    for Codon in RawCount[AA]:
                        if Sum[AA] != 0:
                            relative_Codon_freq = RawCount[AA][Codon]/Sum[AA]
                            relative_Codon_freqTable[AA][Codon] = relative_Codon_freq
                        else:
                            relative_Codon_freq = 0
                            relative_Codon_freqTable[AA][Codon] = relative_Codon_freq   
                     
            
                if Verbose == True:    
                    print("table of total instances of amino acid (from codon count)")
                    print(Sum)
                    print("------------------------")
                    print("Probability of G/C for any nucleotide")
                    print(GC_total_prob)
                    print("------------------------")
                    print("expectation of relative codon frequencies from GC content")
                    print(ProbTable)   
                    print("------------------------")
                    print("Table of relative codon frequencies")
                    print(relative_Codon_freqTable)   
                    print("------------------------")
                        
                        
             ######## Once all Prob_i values  and relative codon frequencies are found, we can calculate the X^2 and F values
             
             #X^2
                for AA in ProbTable:
                    chi_squared = 0
                    for Codon in ProbTable[AA]:
                        if ProbTable[AA] != 0: #Novembre paper is unclear as to whether p_{i} is from relative or total frequencies- values only make sense from total
                            chi_squared += Sum[AA]*(Codon_freqTable_raw[AA][Codon] - ProbTable[AA][Codon])*(Codon_freqTable_raw[AA][Codon] - ProbTable[AA][Codon])/(ProbTable[AA][Codon])
                            #chi_squared += Sum[AA]*(relative_Codon_freqTable[AA][Codon] - ProbTable[AA][Codon])*(relative_Codon_freqTable[AA][Codon] - ProbTable[AA][Codon])/(ProbTable[AA][Codon])
                        else:
                            chi_squared = 0
                    Chi_squared_AA[AA] = chi_squared

            #F_Wright DOUBLE CHECK WITH JENNY 
                for AA in relative_Codon_freqTable:
                    sub_sum = 0
                    if Sum[AA] != 1:
                        for Codon in relative_Codon_freqTable[AA]:
                            sub_sum += Sum[AA]*relative_Codon_freqTable[AA][Codon]*relative_Codon_freqTable[AA][Codon]
                        F_wright =(sub_sum -1)/(Sum[AA] -1)
                        F_Wright_AA[AA]=F_wright
                    else:
                       F_Wright_AA[AA] = 0 
                    
            #F_Novembre
                for AA in Chi_squared_AA:
                    k = len(Codon_freqTable_raw[AA])
                    if Sum[AA] != 1:
                        F_nov = (Chi_squared_AA[AA] + Sum[AA] - k)/(k*(Sum[AA] -1))
                        F_Novembre_AA[AA]=F_nov
                    else:
                        F_Novembre_AA[AA] = 0
                    
            
                if Verbose == True: 
                    print("Chi_squared Table")
                    print(Chi_squared_AA)
                    print("------------------------")
                    print("Novembre F values")
                    print(F_Novembre_AA )
                    print("------------------------")
                    print("Wright F values")
                    print(F_Wright_AA)
                    print("------------------------")
                    for AA in relative_Codon_freqTable:
                        print("codon degeneracy of %s = %s"%(AA,len(Codon_freqTable_raw[AA])))
                    sys.stdout.flush()
        

            #Nc_wright
                Nc_wright = 0
                F2w=0
                F3w=0
                F4w=0
                F6w=0
                for AA in F_Wright_AA :
                    if len(Codon_freqTable_raw[AA])==2:
                        F2w += F_Wright_AA[AA]/9
                    elif len(Codon_freqTable_raw[AA])==3:
                        F3w += F_Wright_AA[AA]/2
                    elif len(Codon_freqTable_raw[AA])==4:
                        F4w += F_Wright_AA[AA]/5
                    elif len(Codon_freqTable_raw[AA])==6:
                        F6w += F_Wright_AA[AA]/3
                    else:
                        pass   
                if Verbose == True:    
                    print("Wright F value averages")
                    print("------------------------") 
                    print("F2= %s"%F2w)
                    print("------------------------")    
                    print("F3= %s"%F3w)
                    print("------------------------")    
                    print("F4= %s"%F4w)     
                    print("------------------------")    
                    print("F6= %s"%F6w)  
                    sys.stdout.flush() 
                    
                Nc_wright = 2+ 9/F2w + 2/F3w + 5/F4w + 3/F6w
                
                
            #Nc_Novembre
                Nc_nov = 0
                F2n=0
                F3n=0
                F4n=0
                F6n=0
                for AA in F_Wright_AA :
                    if len(Codon_freqTable_raw[AA])== 2:
                        F2n += F_Novembre_AA[AA]/9
                    elif len(Codon_freqTable_raw[AA])==3:
                        F3n += F_Novembre_AA[AA]/2
                    elif len(Codon_freqTable_raw[AA])==4:
                        F4n += F_Novembre_AA[AA]/5
                    elif len(Codon_freqTable_raw[AA])==6:
                        F6n += F_Novembre_AA[AA]/3
                    else:
                        pass   
                if Verbose == True: 
                    print("Novembre F value averages")
                    print("------------------------") 
                    print("F2= %s"%F2n)
                    print("------------------------")    
                    print("F3= %s"%F3n)
                    print("------------------------")    
                    print("F4= %s"%F4n)     
                    print("------------------------")    
                    print("F6= %s"%F6n)  
                    sys.stdout.flush()             
                Nc_nov = 2+ 9/F2n + 2/F3n + 5/F4n + 3/F6n               
                    
            


                print("%s,%s,%s"%(i,Nc_wright,Nc_nov))
        
            
cnx.close()
