##!/usr/bin/env python2
## -*- coding: utf-8 -*-
#"""
#Created on Fri Oct 25 11:16:41 2019
#
#@author: lu
#"""
#
import findSequenceById
import windowPhi as methods
import numpy as np
import scipy.stats as ss
import math
import random as rd

rd.seed(0)
tissueSamples=25

methods.deltaEtaFile="/home/lu/Desktop/data_modified/Crei_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Crei_Mutation.csv"
sequenceDict=findSequenceById.findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")




f=open("cEl.csv","r+")
lines=f.readlines()

#21 tissues or lifestages

tissueMatrix=[]
phiList=[]

for i in range(tissueSamples):
    tissueMatrix.append([])

for line in lines:
    randomN=rd.randrange(1,101)
    if randomN>=10:
        continue
    line=line.rstrip()
    splitList=line.split(",")
    geneName=splitList[0]
    if geneName in sequenceDict:
        phi=methods.calPhiForGene(sequenceDict[geneName])
        if math.isnan(phi):
            phiList.append(0)  
        else:
            phiList.append(phi)
    else:
        continue
    for i in range(1,len(splitList)):
        if len(splitList[i])==0:
            tissueMatrix[i-1].append(0)
        else:
            tissueMatrix[i-1].append(float(splitList[i]))

def corelate(list1,list2):
    return ss.pearsonr(list1,list2)[0]


for tissueList in tissueMatrix:
    print ss.pearsonr(tissueList,phiList)

larva_L1=tissueMatrix[7]
adult=tissueMatrix[12]
embryos=tissueMatrix[16]
corelationList=[]

weights=range(-10,10)
for i in weights:
    for j in weights:
        for p in weights:
            newList1=[k*i for k in larva_L1]
            newList2=[k*j for k in adult]
            newList3=[k*p for k in embryos]
            sumList=[(x + y)/2.0 for x, y in zip(newList1, newList2)]
            corelation=ss.pearsonr(phiList,sumList)[0]
            if not corelation==0:
                corelationList.append(corelation)


print corelate(larva_L1,phiList)
print corelate(adult,phiList)
print max(corelationList)
    











