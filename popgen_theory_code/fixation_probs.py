# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:45:15 2020

@author: Catarina
"""

import matplotlib.pyplot as plt
import numpy as np

def probfix(popsize,s):
    
    #probability of fixationof weakly deleterious mutations
    probfix = (1-np.exp(-0.5*s))*((1-np.exp(-popsize*s))**-1)
    
    return probfix

def delfix(popsize,s):
    
    #probability of fixation of weakly deleterious mutations
    delfix = (1-np.exp(0.5*s))*((1-np.exp(popsize*s))**-1)
    
    return delfix

def ratio(popsize,s):
    
    #probability of fixationof weakly deleterious mutations
    probfix = (1-np.exp(-0.5*s))*((1-np.exp(-popsize*s))**-1)
    
    #probability of fixation of weakly deleterious mutations
    delfix = (1-np.exp(0.5*s))*((1-np.exp(popsize*s))**-1)
    
    #ratio of fixation of deleterious/nondeleterious mutations
    R = delfix*(probfix**-1)
  
    return R

def flat(popsize):
    f = 1
    return f

#selection coefficients
s = 1
s1 = 0.01
s2=0.012
s3=0.013
s4=0.014
s5 = 0.015
s6=0.016
s7=0.017
s8=0.018
s9 = 0.05

#population size ranges
popsize1=np.arange(10, 200, 0.1)
popsize=np.arange(0, 600, 0.1)
_popsize=np.arange(0, 6000, 1)
popsize2=np.arange(2100, 2500, 0.1)

#selection coefficient ranges
select=np.arange(-0.0025, 0.05, 0.0005)

#plots
plt.plot( select, probfix(20000, select), 'r', label = 'N = 20000, beneficial')
plt.plot(select, delfix(20000, select),'r--',label = 'N = 20000, deleterious')
plt.plot( select, probfix(50, select), 'c', label = 'N = 50, beneficial')
plt.plot(select, delfix(50, select),'c--',label = 'N = 50, deleterious')
plt.plot( select, probfix(200, select), 'b', label = 'N = 200, beneficial')
plt.plot(select, delfix(200, select),'b--',label = 'N = 200, deleterious')
plt.plot( select, probfix(20, select), 'k',label = 'N = 20, beneficial')
plt.plot(select, delfix(20, select),'k--', label = 'N = 20, deleterious')
plt.legend(loc="upper right")
plt.title('Probality of Fixation depending on Selection Coefficient')
plt.xlabel('Selection Coefficient s')
plt.ylabel('Probability of Fixation')
plt.show()

plt.plot(popsize1, probfix(popsize1, s1), 'r',label = 's = 0.01, beneficial')
plt.plot(popsize1, delfix(popsize1, s1), 'r--',label = 's = 0.01, deleterious')
plt.plot(popsize1, probfix(popsize1, s8), 'b',label = 's = 0.018, beneficial')
plt.plot(popsize1, delfix(popsize1, s8), 'b--',label = 's = 0.018, deleterious')
plt.plot(popsize1, probfix(popsize1, s9), 'k',label = 's = 0.05, beneficial')
plt.plot(popsize1, delfix(popsize1, s9), 'k--',label = 's = 0.05, deleterious')
plt.legend(loc="upper right")
plt.title('Probality of Fixation depending on Population size')
plt.xlabel('Population Size N')
plt.ylabel('Probability of Fixation')
plt.show()


plt.plot(popsize, ratio(popsize, s1), 'r--',label = 's = 0.01')
plt.plot(popsize, ratio(popsize, s2), 'k--',label = 's = 0.012')
plt.plot(popsize, ratio(popsize, s3), 'y--',label = 's = 0.013')
plt.plot(popsize, ratio(popsize, s4), 'g--',label = 's = 0.014')
plt.plot(popsize, ratio(popsize, s5), 'c--',label = 's = 0.015')
plt.plot(popsize, ratio(popsize, s6), 'b--',label = 's = 0.016')
plt.plot(popsize, ratio(popsize, s7), 'm--',label = 's = 0.017')
plt.plot(popsize, ratio(popsize, s8), 'r-',label = 's = 0.018')
plt.legend(loc="upper right")
plt.title('Ratio of Fixation as function of Population size')
plt.xlabel('Population Size N')
plt.ylabel('Ratio of fixation \n of deleterious vs beneficial mutations')
plt.show()


plt.plot(popsize*s1, ratio(popsize, s1), 'r--',label = 's = 0.01')
plt.plot(popsize*s2, ratio(popsize, s2), 'k--',label = 's = 0.012')
plt.plot(popsize*s3, ratio(popsize, s3), 'y--',label = 's = 0.013')
plt.plot(popsize*s4, ratio(popsize, s4), 'g--',label = 's = 0.014')
plt.plot(popsize*s5, ratio(popsize, s5), 'c--',label = 's = 0.015')
plt.plot(popsize*s6, ratio(popsize, s6), 'b--',label = 's = 0.016')
plt.plot(popsize*s7, ratio(popsize, s7), 'm--',label = 's = 0.017')
plt.plot(popsize*s8, ratio(popsize, s8), 'r-',label = 's = 0.018')
plt.legend(loc="upper right")
plt.title('Ratio of Fixation as function of Population size')
plt.xlabel('sN')
plt.ylabel('Ratio of fixation \n of deleterious vs beneficial mutations')
plt.show()

plt.plot(popsize*0.001, ratio(_popsize, 0.001), 'r--',label = 's = 0.001')
plt.plot(popsize*0.0001, ratio(_popsize, 0.0001), 'y--',label = 's = 0.0001')
plt.plot(popsize*0.00001, ratio(_popsize, 0.00001), 'g--',label = 's = 0.00001')
plt.plot(popsize*0.000001, ratio(_popsize, 0.000001), 'c--',label = 's = 0.000001')
plt.legend(loc="upper right")
plt.title('Ratio of Fixation as function of Population size times s')
plt.xlabel('sN')
plt.ylabel('Ratio of fixation \n of deleterious vs beneficial mutations')
plt.show()

plt.plot(popsize2, ratio(popsize2, 0.01), 'r--',label = 's = 0.01')
plt.plot(popsize2, ratio(popsize2, 0.00998), 'k--',label = 's = 0.00998')
plt.plot(popsize2, ratio(popsize2, 0.00996), 'y--',label = 's = 0.00996')
plt.plot(popsize2, ratio(popsize2, 0.00994), 'g--',label = 's = 0.00994')
plt.plot(popsize2, ratio(popsize2, 0.00992), 'm--',label = 's = 0.00992')
#plt.hlines(y=1,xmin=27000,xmax=28000, color='c', linestyle='-',label =' - and + fixation equal probability ')
plt.legend(loc="upper right")
plt.title('Ratio of Fixation as function of Population size \n Zoomed in')
plt.xlabel('Population Size N')
plt.ylabel('Ratio of fixation \n of deleterious vs beneficial mutations (1e-9)') 
plt.show()

plt.plot(popsize2*0.01, ratio(popsize2, 0.01), 'r--',label = 's = 0.01')
plt.plot(popsize2*0.00998, ratio(popsize2, 0.00998), 'k--',label = 's = 0.00998')
plt.plot(popsize2*0.00996, ratio(popsize2, 0.00996), 'y--',label = 's = 0.00996')
plt.plot(popsize2*0.00994, ratio(popsize2, 0.00994), 'g--',label = 's = 0.00994')
plt.plot(popsize2*0.00992, ratio(popsize2, 0.00992), 'm--',label = 's = 0.00992')
#plt.hlines(y=1,xmin=27000,xmax=28000, color='c', linestyle='-',label =' - and + fixation equal probability ')
plt.legend(loc="upper right")
plt.title('Ratio of Fixation as function of Population size \n Zoomed in')
plt.xlabel('sN')
plt.ylabel('Ratio of fixation \n of deleterious vs beneficial mutations (1e-9)') 
plt.show()

#plt.plot(select, ratio(2000, select), 'r--',label = 'N = 2000')
#plt.plot(select, ratio(5000,select), 'k--',label = 'N = 5000')
#plt.plot(select, ratio(10000, select), 'y--',label = 'N = 10000')
#plt.plot(select, ratio(20000, select), 'g--',label = 'N = 20000')
#plt.legend(loc="upper left")
#plt.title('Ratio of Fixation as function of selection coefficents, Zoomed in')
#plt.xlabel('Selection Coefficient s')
#plt.ylabel('Ratio of fixation \n of nondeleterious vs deleterious mutations') 
#plt.show()
