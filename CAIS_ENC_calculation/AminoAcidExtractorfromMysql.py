# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:44:49 2020

@author: Catarina
"""

import mysql.connector
from itertools import chain


def LoadConnection():
	Database='Database'
	user= 'username'
	host= 'DNS'
	password = password'
	return Database, user,host,password

#################################Write the heading#######################
with open("AA_for_linearmodeling_species_differences.txt", "a") as f1:
            f1.write( 'SpeciesUID,PfamUID,A,R,N,D,C,E,Q,G,H,O,I,L,K,M,F,P,U,S,T,W,Y,V' + '\n')

f1.close()	

##################connect to MYSQL##################################
cnx = mysql.connector.connect(user = 'user',password= 'passwrd' ,host='DNS' , database= 'PFAMphylostratigraphy')

mycursor = cnx.cursor(buffered = True)


####################################limiting set of extracted species to keys in the dictionary below#######################
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

#####################BEGIN THE EXTRACTION#################################

for i in range(4,507):
    if i in Species_Total_GC_content.keys():
        print("extracting species %s"%i)
        #extractor
#*=entire row
        extractionStatement= 'SELECT SpeciesUID,PfamUID,PercentAminoAcidComposition FROM Genomes_Multicellular_DomainMetrics WHERE SpeciesUID =%s'%i

        mycursor.execute(extractionStatement)
        
        #########################################WRITE THE EXTRACTED STUFF TO THE FILE###################

        for entry in mycursor:
	
            string = ','.join(map(str,entry))
            print(string)
            
            with open("AA_for_linearmodeling_species_differences.txt", "a") as f1:
                f1.write(string + '\n')

f1.close()
