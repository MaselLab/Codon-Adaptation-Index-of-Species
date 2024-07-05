"""

@author: Andrew
"""

import os, sys, json, csv, mysql.connector,datetime,math

"""
The purpose of this file is to determine the codon adaptation index of species (CAIS) for individual species.
This is a metric that describes codon bias patterns. 
CAIS is calculated in the following way:
    
    I) Calculate the probability of a codon given the genomic GC content and amino acid frequencies (null expectation)
    
    II) Count the number of each codon in a sequence
        
    II) Calculate CAIS metric for the species
         
"""


####################################################################################
#                                 User Information                                 #
####################################################################################

# SQL Connection Information
Database =
User = 
Host =
Password =

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
    
    #The amino acid frequency values in this table are calculated in another script (WHICH TAKES & HOURS TO RUN)
    #USED FOR THE RE_WEIGHTING OF CAIS

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
    
    #genome wide GC content for each species
    Species_Total_GC_content = {3:0.53,
                                4:0.4929,
                                5:0.4842,
                                6:0.515,
                                8:0.4332,
                                9:0.5129,
                                10:0.5258,
                                11:0.5323,
                                12:0.5123,
                                13:0.5258,
                                14:0.5269,
                                15:0.5167,
                                16:0.5166,
                                17:0.5325,
                                18:0.5205,
                                19:0.4266,
                                20:0.5116,
                                22:0.5076,
                                23:0.4983,
                                24:0.5257,
                                25:0.5254,
                                26:0.5405,
                                27:0.5052,
                                28:0.5161,
                                29:0.5274,
                                30:0.5206,
                                31:0.5187,
                                32:0.5579,
                                33:0.5054,
                                34:0.561,
                                35:0.5152,
                                36:0.5233,
                                38:0.5222,
                                39:0.5152,
                                40:0.5273,
                                41:0.4667,
                                42:0.5302,
                                43:0.5184,
                                44:0.5139,
                                45:0.5146,
                                46:0.5168,
                                47:0.5107,
                                48:0.4871,
                                49:0.5254,
                                50:0.5283,
                                51:0.5253,
                                52:0.4839,
                                53:0.5222,
                                54:0.5169,
                                55:0.5244,
                                56:0.5225,
                                57:0.5314,
                                58:0.5329,
                                59:0.5165,
                                60:0.5098,
                                61:0.5317,
                                62:0.5385,
                                63:0.5386,
                                64:0.528,
                                65:0.5162,
                                66:0.5264,
                                67:0.5132,
                                68:0.5167,
                                69:0.5273,
                                70:0.5052,
                                71:0.5155,
                                72:0.4834,
                                73:0.5293,
                                74:0.5956,
                                75:0.5301,
                                76:0.5181,
                                77:0.5175,
                                78:0.5153,
                                79:0.3961,
                                80:0.5146,
                                81:0.4804,
                                82:0.5269,
                                83:0.5082,
                                84:0.5402,
                                85:0.5596,
                                86:0.4629,
                                87:0.5346,
                                135:0.5295,
                                136:0.4552,
                                137:0.4256,
                                138:0.4389,
                                139:0.4288,
                                140:0.521,
                                141:0.4601,
                                142:0.4595,
                                143:0.4622,
                                145:0.6979,
                                146:0.5327,
                                147:0.4351,
                                148:0.5657,
                                149:0.4338,
                                150:0.3978,
                                151:0.4339,
                                152:0.4342,
                                153:0.4342,
                                154:0.5095,
                                155:0.4236,
                                156:0.4304,
                                157:0.4075,
                                158:0.4962,
                                159:0.5226,
                                160:0.5221,
                                161:0.6103,
                                163:0.5286,
                                164:0.5276,
                                165:0.5305,
                                166:0.5864,
                                167:0.4375,
                                168:0.4913,
                                169:0.435,
                                170:0.4431,
                                171:0.5434,
                                173:0.4185,
                                174:0.5316,
                                175:0.4349,
                                176:0.3968,
                                177:0.5354,
                                178:0.5252,
                                179:0.4471,
                                180:0.4893,
                                181:0.3909,
                                182:0.4992,
                                183:0.4057,
                                184:0.5417,
                                185:0.5652,
                                186:0.3905,
                                187:0.4302,
                                188:0.4746,
                                189:0.4779,
                                190:0.3929,
                                191:0.4383,
                                192:0.446,
                                194:0.4751,
                                195:0.4293,
                                196:0.4896,
                                197:0.4565,
                                198:0.5539,
                                199:0.4664,
                                200:0.473,
                                201:0.4608,
                                202:0.5482,
                                203:0.5471,
                                204:0.5196,
                                206:0.5289,
                                208:0.5633,
                                292:0.6274,
                                294:0.4401,
                                295:0.5145,
                                296:0.402,
                                297:0.3911,
                                298:0.4032,
                                299:0.5075,
                                300:0.4546,
                                301:0.4463,
                                302:0.454,
                                303:0.4938,
                                304:0.5377,
                                305:0.5296,
                                306:0.4516,
                                307:0.4737,
                                308:0.451,
                                309:0.425,
                                310:0.4388,
                                311:0.4605,
                                312:0.4382,
                                313:0.5346,
                                314:0.5653,
                                315:0.5548,
                                316:0.5135,
                                317:0.5594,
                                318:0.5233,
                                319:0.5576,
                                320:0.5572,
                                321:0.561,
                                322:0.5367,
                                323:0.5582,
                                324:0.5455,
                                325:0.545,
                                326:0.5434,
                                327:0.5598,
                                328:0.4521,
                                329:0.5014,
                                330:0.4238,
                                331:0.4346,
                                332:0.5233,
                                333:0.4324,
                                334:0.341,
                                335:0.4103,
                                336:0.4747,
                                337:0.4665,
                                338:0.4637,
                                339:0.4458,
                                340:0.3965,
                                341:0.4639,
                                342:0.4805,
                                343:0.3883,
                                344:0.4503,
                                345:0.5014,
                                346:0.477,
                                347:0.4751,
                                348:0.4727,
                                349:0.3888,
                                350:0.4573,
                                351:0.4502,
                                352:0.5341,
                                353:0.4676,
                                354:0.4181,
                                355:0.4441,
                                356:0.4442,
                                357:0.448,
                                358:0.447,
                                359:0.4968,
                                360:0.4946,
                                361:0.482,
                                362:0.458,
                                363:0.4898,
                                364:0.4453,
                                365:0.4425,
                                366:0.4216,
                                367:0.4391,
                                368:0.4106,
                                369:0.4336,
                                370:0.4317,
                                371:0.4311,
                                372:0.475,
                                373:0.4581,
                                374:0.474,
                                375:0.452,
                                376:0.4489,
                                377:0.433,
                                378:0.4329,
                                379:0.4273,
                                380:0.4503,
                                381:0.4204,
                                382:0.4452,
                                383:0.4185,
                                384:0.4525,
                                385:0.6659,
                                386:0.4467,
                                387:0.4613,
                                388:0.4433,
                                389:0.4234,
                                390:0.4238,
                                391:0.4243,
                                392:0.5241,
                                393:0.5881,
                                394:0.4799,
                                395:0.4348,
                                396:0.4439,
                                397:0.4461,
                                398:0.4553,
                                399:0.4607,
                                400:0.4253,
                                401:0.451,
                                402:0.4199,
                                403:0.4286,
                                404:0.4785,
                                405:0.4344,
                                407:0.5222,
                                408:0.5155,
                                409:0.5274,
                                410:0.527,
                                411:0.5301,
                                412:0.5244,
                                413:0.5234,
                                414:0.5188,
                                415:0.5171,
                                416:0.5221,
                                417:0.5114,
                                418:0.5255,
                                419:0.5401,
                                420:0.5231,
                                421:0.5372,
                                422:0.517,
                                423:0.5133,
                                424:0.5299,
                                425:0.5326,
                                426:0.526,
                                427:0.5329,
                                428:0.5304,
                                429:0.5156,
                                430:0.5339,
                                431:0.5194,
                                433:0.5165,
                                434:0.5203,
                                435:0.5338,
                                436:0.5406,
                                437:0.5275,
                                438:0.5271,
                                439:0.5332,
                                440:0.5136,
                                442:0.5363,
                                444:0.4918,
                                445:0.5222,
                                446:0.5121,
                                447:0.5145,
                                448:0.5166,
                                449:0.5199,
                                450:0.5236,
                                452:0.5476,
                                453:0.4996,
                                454:0.504,
                                455:0.515,
                                456:0.5281,
                                457:0.491,
                                458:0.5371,
                                459:0.4917,
                                460:0.4857,
                                461:0.4789,
                                462:0.4759,
                                463:0.482,
                                464:0.488,
                                465:0.4915,
                                466:0.4773,
                                467:0.5201,
                                468:0.4869,
                                469:0.5038,
                                470:0.4872,
                                471:0.4923,
                                472:0.4766,
                                473:0.4917,
                                474:0.4875,
                                475:0.4871,
                                476:0.478,
                                477:0.4885,
                                478:0.5561,
                                479:0.481,
                                480:0.4905,
                                481:0.4922,
                                482:0.4944,
                                483:0.4804,
                                484:0.4914,
                                485:0.5212,
                                486:0.5023,
                                487:0.4876,
                                488:0.5466,
                                489:0.4883,
                                490:0.4845,
                                491:0.4863,
                                492:0.4754,
                                493:0.5517,
                                494:0.4709,
                                495:0.4771,
                                496:0.5013,
                                497:0.5079,
                                498:0.4638,
                                499:0.4892,
                                500:0.5224,
                                501:0.5513,
                                502:0.5156,
                                503:0.5348,
                                504:0.5356,
                                505:0.5301,
                                506:0.5017,
                                507:0.4844,
                                508:0.5064,
                                509:0.4957,
                                510:0.5204,
                                511:0.4885,
                                512:0.4893,
                                513:0.4906,
                                514:0.4732,
                                515:0.5212,
                                516:0.4807,
                                517:0.4976,
                                518:0.529,
                                519:0.5392,
                                520:0.4966,
                                521:0.5446,
                                522:0.4856,
                                523:0.5047,
                                524:0.4769,
                                525:0.4754,
                                526:0.4709,
                                527:0.506,
                                528:0.5411,
                                529:0.4831,
                                530:0.4743,
                                531:0.5135,
                                532:0.4786,
                                533:0.5235,
                                534:0.5218,
                                535:0.479,
                                536:0.5526,
                                537:0.5469,
                                538:0.5045,
                                539:0.5354,
                                540:0.5091,
                                541:0.5074,
                                542:0.5091,
                                543:0.5436,
                                544:0.4721,
                                545:0.5047,
                                546:0.4824,
                                547:0.4643,
                                548:0.4654,
                                549:0.4576,
                                550:0.4975}

############################ACTUAL CAIS CALCULATING CODE#####################################
    
    for i in range(0,550):
        if i in Species_Total_GC_content.keys():
            if Verbose == True:
                print('Extracting SpeciesUID: %s'%i)
                sys.stdout.flush()
        
            GC_total_prob = Species_Total_GC_content[i]
            notGC_total_prob = 1-GC_total_prob


        # We'll keep dictionaries of all the codons that we count, and their expected probabilities,
        # sorted by which amino acid they correspond to. This will make calculating the CAIS relatively easy and painless- Sara W

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
            Sum_ProbTable = {'F':0,
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


######################################

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

            # Calculate expected distribution given GC content and amino acid frequency
                Summed_codon_instances = 0 
 
                for AA in RawCount:
                    Sum[AA] = sum(RawCount[AA].values())
                    for Codon in RawCount[AA]:
                        if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                            Prob = (notGC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
                        elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                            Prob = (GC_total_prob*GC_total_prob*GC_total_prob)*0.125
                        elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC','TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                            Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
                        elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT','AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                            Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)*0.125
                        else:
                            print("storm the castle")
                        ProbTable[AA][Codon] = Prob

                for AA in ProbTable:
                    Sum_ProbTable[AA]= sum(ProbTable[AA].values())
                    for Codon in ProbTable[AA]:
                        Expected[AA][Codon] = ProbTable[AA][Codon]/Sum_ProbTable[AA]
                            
            #calculate observed codon frequencies
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    Observed[AA][Codon] = RawCount[AA][Codon]/Sum[AA]

            #CAIS
            CAIS = 0
            for AA in RawCount:
                aa_weight = Total_AA_freqTable[AA]/sum(Observed[AA].values())
                print(aa_weight)
                for Codon in RawCount[AA]:
                    CAIS += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])*aa_weight

            #CAIS_unweighted
            CAIS_unweighted = 0
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    CAIS_unweighted += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])

            if Verbose == True: 
                print("-------SumTable-----------------")    
                print(Sum)      
                print("-------ProbTable-----------------")    
                print(ProbTable)
                print("-------Sum_ProbTable--------------")
                print(Sum_ProbTable)
                print("-------Expected-----------------")
                print(Expected)
                print("-------Observed-----------------")
                print(Observed)
                sys.stdout.flush()        
            
            #Print results for each species
            print("%s,%s,%s"%(i,CAIS,CAIS_unweighted))
            
cnx.close()

                        
                        
             
