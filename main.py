#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 13:07:28 2019

@author: lu
"""
import findSequenceById
import calculateTAI
from Bio.SeqUtils import GC
import validator
import matplotlib.pyplot as plt
import calculateCAI
import math
import random
from sklearn.preprocessing import minmax_scale

lowExpCut="off"

def removeMarked(list1,markedList)
    myList=[]
    for i in range(len(list1)):
        if i not in markedList:
            myList.append(list1[i])
    return myList

def meanGC(lst): #takes in a list of sequence
    gcList=[]
    for seq in lst:
        gcList.append(GC(seq))
    return ss.mstats.gmean(gcList)

#def run(sequenceFileName,expressionFile,speciesName,geneIdType="locus_tag"): #example: "atha.fasta", "atha_mRNAExpressionData.csv", "atha"
#    geneDict=findSequenceById.findSequenceByID(sequenceFileName,idType=geneIdType)
#    f=open(expressionFile)
##    speciesName="atha"
#    expressionDict=dict() #key is the sequence and the expression value is the result
#    lines=f.readlines()              
#    header=lines[0]
#    lines=lines[1:]
#    for i in range(len(lines)):
#        line=lines[i].rstrip()
#        splitList=line.split(",")
#        geneId=splitList[0].rstrip()
#
#        expression=splitList[1]
#        if expression=="NA":
#            expression==0.01
#        else:
#            expression=float(expression)
#        if geneId in geneDict:
#            sequence=geneDict[geneId]
#            expressionDict[(geneId,sequence)]=expression
#    print expressionDict
#        try:
#            phi=methods.calPhiForGene(sequence)
#        except:
#            markedIndexList.append(i)
#            phi=0
#        phiList.append(float(phi))
#        sequenceList.append(sequence)
#        
#    tAIList=calculateTAI.calculate(sequenceList,speciesName)
#    
#    phiList=removeMarked(phiList,markedIndexList)
#    expressionList=removeMarked(expressionList,markedIndexList)
#    tAIList=removeMarked(tAIList,markedIndexList)
#
#    metricList=[]
#    metricList.append(phiList)
#    metricList.append(tAIList)
#    metricList.append(expressionList)
#    
#    print ("phi vs expression")
#
#    typeSet=set()
#    for phi in phiList:
#        typeSet.add(str(type(phi)))
#    
#    print("tAI vs expression")
##    print len(tAIList)
##    print len(expressionList)
##    print (phiList)    #print phiList
#    from scipy import stats
#    print "="
#    print stats.pearsonr(list(tAIList),(phiList))
#    print "="
#    validator.validate(phiList,expressionList,corelationFunction="spearman",logScale="yes")
#    validator.validate(tAIList,expressionList,corelationFunction="spearman",logScale="yes")
#    validator.validate(tAIList,phiList,corelationFunction="spearman",logScale="yes")


def setUp(sequenceFileName,expressionFile,speciesName,geneIdType="locus_tag"):
    geneDict=findSequenceById.findSequenceByID(sequenceFileName,idType=geneIdType)
    f=open(expressionFile)
#    speciesName="atha"
    expressionDict=dict() #key is the sequence and the expression value is the result
    lines=f.readlines()              
    header=lines[0]
    lines=lines[1:]
    for i in range(len(lines)):
        line=lines[i].rstrip()
        splitList=line.split(",")
        geneId=splitList[0].rstrip()
        expression=splitList[1]
        if 'NA' in expression or len(expression)==0 or float(expression)==0:
            continue
            expression=0.01
        else:
            expression=float(expression)
#species="yeast"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Scer_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Scer_Mutation.csv"
##run("/home/lu/Desktop/sequences/s288c.fasta", "/home/lu/Desktop/Expression/Yeast/yeast_mrna_expression.csv", "yeast",geneIdType="locus_tag")
#expressionDict=setUp("/home/lu/Desktop/sequ
        if geneId in geneDict:
            sequence=geneDict[geneId]
            expressionDict[(geneId,sequence)]=expression
    return expressionDict
    
          
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


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def replaceStrageValues(myList):
    newList=[]
    for i in myList:
        if not is_number(i):
            newList.append(0.001)
        elif math.isnan(i):
            newList.append(0.001)
        else:
            newList.append(i)
    return newList
        
            



import methods

#run("atha.fasta", "atha_mRNAExpressionData.csv", "atha")
        



#run("/home/lu/Desktop/sequences/dmel.fasta", "/home/lu/Desktop/Expression/Dmelangaster/dmel_mrna_expression.csv", "dmel",geneIdType="gene")


#
#species="yeast"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Scer_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Scer_Mutation.csv"
##run("/home/lu/Desktop/sequences/s288c.fasta", "/home/lu/Desktop/Expression/Yeast/yeast_mrna_expression.csv", "yeast",geneIdType="locus_tag")
#expressionDict=setUp("/home/lu/Desktop/sequences/s288c.fasta", "/home/lu/Desktop/Expression/Yeast/yeast_mrna_expression.csv", "yeast",geneIdType="locus_tag")
#
#species="ecoli"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Ecoli_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Ecoli_Mutation.csv"
#expressionDict=setUp("/home/lu/Desktop/sequences/ecoli.fasta", "/home/lu/Desktop/Expression/EColi/Mrna_Expression.csv", "ecoli",geneIdType="gene")
#
#
##
species="fruitfly"
methods.deltaEtaFile="/home/lu/Desktop/data_modified/Drer_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Drer_Mutation.csv"
expressionDict=setUp("/home/lu/Desktop/sequences/dmel.fasta", "/home/lu/Desktop/Expression/Dmelangaster/dmel_mrna_expression.csv", "",geneIdType="gene")


#species="frog"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Xtro_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Xtro_Mutation.csv"
#expressionDict=setUp("/home/lu/Desktop/sequences/xtro.fasta", "/home/lu/Desktop/Expression/x.tropicalis/xtro_mrna_expression.csv", "frog",geneIdType="gene")

#species="crei"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Crei_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Crei_Mutation.csv"
#expressionDict=setUp("/home/lu/Desktop/sequences/c_elegan.fasta", "/home/lu/Desktop/Expression/cElegans/mRNA_Crei.csv", "crei",geneIdType="gene")

#species="human"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Crei_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Crei_Mutation.csv"
#expressionDict=setUp("/home/lu/Desktop/data_modified/human.fasta", "/home/lu/Desktop/Expression/cElegans/mRNA_Crei.csv", "crei",geneIdType="gene")
#
#




#/home/lu/Desktop/sequences/c_elegan.fasta



#species="atha"
#methods.deltaEtaFile="/home/lu/Desktop/data_modified/Atha_Selection.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Atha_Mutation.csv"
#expressionDict=setUp("atha.fasta", "/home/lu/Desktop/Expression/Arabidopsis/atha_mRNAExpressionData.csv", "atha")
#



codonOrder=['ATA','ATC', 'ATT', 'ATG',
'ACA', 'ACC', 'ACG', 'ACT',
'AAC', 'AAT', 'AAA', 'AAG',
'AGC', 'AGT', 'AGA', 'AGG',
'CTA', 'CTC', 'CTG', 'CTT',
'CCA', 'CCC', 'CCG', 'CCT',
'CAC', 'CAT', 'CAA', 'CAG',
'CGA', 'CGC', 'CGG', 'CGT',
'GTA', 'GTC', 'GTG', 'GTT',
'GCA', 'GCC', 'GCG', 'GCT',
'GAC', 'GAT', 'GAA', 'GAG',
'GGA', 'GGC', 'GGG', 'GGT',
'TCA', 'TCC', 'TCG', 'TCT',
'TTC', 'TTT', 'TTA', 'TTG',
'TAC', 'TAT', 'TAA', 'TAG',
'TGC', 'TGT','TGG', 'TGA'
]



referenceList=[]
expressionCut=[]
expressionList=[]
for key in expressionDict:

    sequence=key[1]
    expression=expressionDict[key]
    expressionList.append(float(expression))
    
expressionList=sorted(expressionList)
expressionCut= expressionList[int(0.95*len(expressionList))]
low_expressionCut=expressionList[int(0.90*len(expressionList))]

for key in expressionDict:
    sequence=key[1]
    expression=expressionDict[key]
    if expression>=expressionCut:
        if len(sequence)%3==0:
            referenceList.append(sequence)






expressionList=[]
sequenceList=[]
featureMatrix=[]
testList=[]
phiList=[]
tAIList=[]
CAIList=[]
lengthList=[]
deltaMList=[]
gcList=[]
calculateCAI.setUpReference(referenceList)

phiList=[]
for key in expressionDict:
    sequence=key[1]
    sequenceList.append(sequence)
    expression=expressionDict[key]
    if len(sequence)%3!=0:
        continue
    if "on" in lowExpCut:
        if expression<low_expressionCut:
            continue
    
    if "N" in sequence or "M" in sequence or "K" in sequence:
        continue
    codonList=loadSequence(sequence)
    codonCntDict=dict()
    for codon in codonList:
        if codon in codonCntDict:
            codonCntDict[codon]+=1
        else:
            codonCntDict[codon]=1
    featureList=[]
    
#    for codon in codonOrder:
#        featureList.append(random.randint(0,100))

    for codon in codonOrder:
        if codon in codonCntDict:
            featureList.append(float(codonCntDict[codon]))
        else:
            featureList.append(0)
    gcList.append(GC(sequence))
#    featureList.append(GC(sequence))
    lengthList.append(len(sequence))
#    featureList.append(len(sequence))
    phi=methods.calPhiForGene(sequence)
#    featureList.append(phi)
#    featureList.append()
#    featureList.append(methods.calMForGene(sequence))
#    if phi<0.1:
#        continue


    phiList.append(phi)
    tAIList.append(calculateTAI.calculateOneGene(sequence,species))
#    CAIList.append(calculateCAI.calculateCAIForAGene(sequence))
    featureList=minmax_scale(featureList) # Here does the scale
    featureMatrix.append(featureList)
    expressionList.append(expression)
    

print (len(expressionList))

expressionCuts=[]
sortedExpressionList=sorted(expressionList)
#print sortedExpressionList
binNumbers=4 #This is the number of labels
for i in range(1,binNumbers):
    cutIndex=int((float(i)/binNumbers)*len(sortedExpressionList))
    expressionCuts.append(sortedExpressionList[cutIndex])

expressionLabel=[]
expressionCuts.append(9223372036854775807)#maxInt in the end
for expression in expressionList:
    for i in range(len(expressionCuts)):
        expressionCut=expressionCuts[i]
        if expression<expressionCut:
            expressionLabel.append(i) # here imagine we have [0.3, 0.8, 1.8], since its sorted, we go from left to right
            #the first one we see that's larger' index is the label
            break


import scipy.stats as ss

tAIRank=ss.rankdata(tAIList)
phiRank=ss.rankdata(phiList)
expressionRank=ss.rankdata(expressionList)

lowerSequence=[] #lower is good for phi
higherSequence=[]

print (len(phiRank))

diffList=[]
for i in range(len(phiRank)):
    diff=abs(tAIRank[i]-expressionRank[i])-abs(phiRank[i]-expressionRank[i])
    diffList.append(diff)
    
diffList=sorted(diffList)
topSampleList=[1000,2000,3000,4000,5000,6000,7000,8000,9000,100000]

for topSamples in topSampleList:
    diffCut=(diffList[topSamples])    #choose 100 genes with the most prediction difference
    max100SequenceList=[]
    max100GC_List=[]
    for i in range(len(phiRank)):
        diff=abs(phiRank[i]-expressionRank[i])-abs(tAIRank[i]-expressionRank[i])
        if diff<=diffCut:
            max100SequenceList.append(sequenceList[i])
            max100GC_List.append(GC(sequenceList[i]))
    print (len(max100SequenceList))
    print("GC for the top selected %d"%topSamples)
    print(meanGC(max100SequenceList))
    




#
##print expressionList[:10]
##print expressionLabel[:10]
#import pandas as pd
#import confusionMatrixPlot
#from sklearn.pipeline import Pipeline
#from sklearn.feature_selection import SelectFromModel
#from sklearn.svm import SVC
#from sklearn.svm import LinearSVC
#import numpy as np
#from sklearn.model_selection import train_test_split 
#from sklearn.metrics import confusion_matrix 
#from sklearn.metrics import accuracy_score
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.neighbors import KNeighborsClassifier 
#from sklearn.tree import DecisionTreeClassifier 
#from sklearn.naive_bayes import GaussianNB 
#from sklearn import svm
#
#import random
#X=featureMatrix
#print (len(featureMatrix))
#print (len(featureMatrix[0]))
#print("-------")
# 
##random.shuffle(X)
#
#y=expressionLabel
##random.shuffle(y)
#
#trainCut=int(len(X)*(3.0/4))
#print ("trainCut: %d"%trainCut)
#testCut=len(X)-trainCut
#print ("testCut: %d"%testCut)
#X_train=X[:trainCut]
#y_train=y[:trainCut]
#X_test=X[testCut:]
#y_test=y[testCut:]
#
#
#
#classifier= KNeighborsClassifier(n_neighbors=binNumbers)
##classifier= GaussianNB() 
##
##classifier= svm.SVC(gamma='scale', decision_function_shape='ovo')
##
#knn=classifier.fit(X_train,y_train)
##
#knn_predictions = knn.predict(X_test)  
#cm = confusion_matrix(y_test, knn_predictions) 
#accuracy=accuracy_score(knn_predictions,y_test)
#print (accuracy)
#
##print (knn_predictions[:20])
##print (y_test[:20])
#confusionMatrixPlot.plot_confusion_matrix(y_test, knn_predictions,np.array(["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"]))
#print (len(gcList))
##
##medianGC=0.420642 #median gc for dmel
##gcList=np.asanyarray(gcList)-medianGC
##gcList=minmax_scale(gcList)
##scaledTAI=minmax_scale(tAIList)
##scaledPhi=minmax_scale(phiList)
##print(expressionList[:100])
##expressionList=minmax_scale(expressionList)
##deltaList= abs(scaledTAI-expressionList) - abs(scaledPhi-expressionList)
##
#if(not len(expressionList)-len(tAIList)==0):
#    print ("wrong assertion------------------------")
#cnt=0.0
#
#
#
#for i in range(len(expressionList)):
#    if abs(scaledPhi[i]-expressionList[i])<=abs(scaledTAI[i]-expressionList[i]):
#        cnt+=1
#print(cnt)
#print (len(expressionList))
#print (cnt/len(expressionList))
#validator.validate(gcList,deltaList,yLabel="Delta of performance difference between tAi and Phi",xLabel="GC bias relative to the entire genome")
##validator.validate(CAIList,gcList)
#validator.validate(tAIList,gcList,xLabel="Delta of performance difference between tAi and Phi",yLabel="GC bias relative to the entire genome")
#sortedTai=np.asarray([scaledTAI for _,scaledTAI in sorted(zip(expressionList,scaledTAI))])
#sortedPhi=np.asarray([scaledPhi for _,scaledPhi in sorted(zip(expressionList,scaledPhi))])
#sortedExpressionList=np.asarray(sorted(expressionList))
#xRank=range(len(sortedTai))
#plt.figure()
#
#plt.plot(xRank,expressionList)
#plt.figure()
#plt.plot(xRank,sortedPhi-sortedExpressionList,color="green")
#plt.figure()
#plt.plot(xRank,sortedTai-sortedExpressionList,color="blue")



#confusionMatrixPlot.plot_confusion_matrix(y_test, tAI_discrete,np.array(["E1","E2","E3","E4","E5","E6"]))

import lu_io

#phiList=replaceStrageValues(phiList)
#lu_io.writeListToCsv(phiList,"test_phi.csv")
#validator.validate(phiList,expressionList,corelationFunction="pearson",logScale="yes",xLabel="MLE-Phi (log)",yLabel="Expression (log)" )
#validator.validate(tAIList,expressionList,corelationFunction="pearson",logScale="yes",xLabel="tAI (log)",yLabel="Expression (log)")
#validator.validate(CAIList,expressionList,corelationFunction="pearson",logScale="yes",xLabel="CAI (log)",yLabel="Expression (log)")
#validator.validate(lengthList,expressionList)
#phiDiscrete=pd.cut((phiList), bins=binNumbers, labels=np.arange(binNumbers), right=False)
#expressionDisrete=pd.cut((expressionList), bins=binNumbers, labels=np.arange(binNumbers), right=False)




#tAIList=calculateTAI.calculate(sequenceList,"yeast")
#validator.validate(tAIList,expressionList)

#methods.deltaEtaFile="/home/lu/Desktop/data_modified/.csv"
#methods.deltaMFile="/home/lu/Desktop/data_modified/Xtro_Mutation.csv"
#run("/home/lu/Desktop/sequences/xtro.fasta", "/home/lu/Desktop/Expression/x.tropicalis/xtro_mrna_expression.csv", "frog",geneIdType="gene")

