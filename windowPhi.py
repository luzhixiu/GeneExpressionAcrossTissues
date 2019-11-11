#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
from statistics import mean
from scipy.stats.mstats import gmean

windowSize=11 #change window size here

deltaEtaFile="/home/lu/Code/selection.csv"
deltaMFile="/home/lu/Code/mutation.csv"


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





phiDict=dict()
def method4(codonList):
    global mDict
    global etaDict   
    global phiDict
    phiList=[]

    for codon in codonList:
        if codon in phiDict:
             phiList.append(phiDict[codon])
        else:
            maxprob=0.0
            selectedPhi=0.0
            rangeList=[]
            for i in range(1,10001):
                rangeList.append(i/100.0)
            for phi in rangeList:
                if codon in mDict:    
                    deltaM=float(mDict[codon])
                else:
                    deltaM=1.0
                if codon in etaDict:
                    deltaEta=float(etaDict[codon])
                else:
                    deltaEta=1.0
                global synonymonDict
                synoList=synonymonDict[codon]
                divisor=np.exp(-1.00*deltaM-(deltaEta*phi))
                dividant=0
                for syno in synoList:
                    if syno in mDict:
                        deltaM=float(mDict[syno])
                    else:
                        deltaM=1.0
                    if syno in etaDict:
                        deltaEta=float(etaDict[syno])
                    else:
                        deltaEta=1.0
                    tmp=np.exp((-1.00*deltaM)-(deltaEta*phi))
    #                print((-1.00*deltaM)-(deltaEta*phi))
    #                print (tmp)
                    dividant+=tmp
                if not dividant==0:
                    prob=divisor/dividant
                else:
                    prob=0
                
                if prob>maxprob:
                    maxprob=prob
                    selectedPhi=phi
            if not selectedPhi==0:
                phiDict[codon]=selectedPhi
                phiList.append(selectedPhi)
            else:
                phiList.append(0.001)
                print ("found 0")
    #        print "Prob: %s"%str(maxprob)
    #        print "Phi: %s"%str(selectedPhi)   
    return phiList

def calPhiForGene(sequence):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)
    avg=gmean(phiList)
    return avg

sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
codonList=loadSequence(sequence)
testResults=method4(codonList)
print(testResults)
print (calPhiForGene(sequence))
