#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 10:08:30 2019

@author: lu
"""
import math
import numpy as np
from scipy import optimize





deltaEtaFile="/home/lu/Code/selection.csv"
deltaMFile="/home/lu/Code/mutation.csv"
codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

inverseTable=  {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
                'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
              'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
              'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
              'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
              'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
              'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}
synonymonDict={}
for key in inverseTable:
    valueList=inverseTable[key]
    for value in valueList:
        synonymonDict[value]=valueList

deltaEtaFile=open(deltaEtaFile)
lines=deltaEtaFile.readlines()
etaDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        etaDict[splitList[1]]=splitList[2]  

deltaMFile=open(deltaMFile)
lines=deltaMFile.readlines()
mDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        mDict[splitList[1]]=splitList[2]

phiDict=dict()
def calculatePhi(codonList,windowSize=9):
    global mDict
    global etaDict   
    global phiDict
    
    
    phiList=[]

    import sys

    rangeList=[]
    for i in range(1,101):
        rangeList.append(i/10.0)
    windCnt=0
    for i in range(len(codonList)): #This one loops through windows
        windCnt+=1
        print "move to window %d"%windCnt
        if i+windowSize>=len(codonList):
            return phiList
        windowCodonList=codonList[i:i+windowSize+1]
        windowProb=1.0
        maxPhi=0.0
        maxProb=-sys.maxsize -1
        for phi in rangeList:
            for codon in windowCodonList:
                if codon not in mDict or codon not in etaDict:
                    windowProb*=1
                else:
                    #calculate top part of division
                    deltaM_i=float(mDict[codon])
                    deltaEta_i=float(etaDict[codon])
                    divisor=np.exp(-1.0*deltaM_i-(deltaEta_i*phi))
#                    print divisor
                    #calculate bot part of division
                    dividant=0.0
                    synoList=synonymonDict[codon]
                    for syno in synoList:
                        if syno not in mDict or syno not in etaDict:
                            break
                        else:
                            dividant+=np.exp((-1*float(mDict[syno]))-(float(etaDict[syno]))*phi)
                    if dividant==0:
                        windowProb*=1
                    else:
                        windowProb*=divisor/dividant
                windowProb=float(windowProb)
            if windowProb>maxProb:
                maxProb=windowProb

                maxPhi=phi
                print "prob",
                print maxProb
                print "phi",
                print maxPhi



            
                        
                        
                        
        
    
    


def loadSequence(sequence):
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    actualCodonList=[]
    started=False
    for codon in codonList:
        if codon in stopCodonList:
            break
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList




def test():
    sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
    codonList=loadSequence(sequence)
    testResults=calculatePhi(codonList)

    
test()


    