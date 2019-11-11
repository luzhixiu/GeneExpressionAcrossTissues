#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:22:45 2019

@author: lu
"""
import findSequenceById
import calculateTAI
import validator


class Record:
    geneId=""
    sequence=""
    expression=""  

    

def run(sequenceFileName,expressionFile,speciesName,geneIdType="locus_tag"): #example: "atha.fasta", "atha_mRNAExpressionData.csv", "atha"
    geneDict=findSequenceById.findSequenceByID(sequenceFileName,idType=geneIdType)
    f=open(expressionFile)
#    speciesName="atha"
    
    lines=f.readlines()
    expressionList=[]
    markedIndexList=[]
    sequenceList=[]
    recordList=[] 
    for i in range(len(lines)):
        line=lines[i].rstrip()
        splitList=line.split(",")
        geneId=splitList[0].rstrip()
        expression=splitList[1]
        try:
            expressionList.append(float(expression))
        except:
            expressionList.append(0)
            markedIndexList.append(i)
        sequence=""
        if geneId in geneDict:
            sequence=geneDict[geneId]
        else:
            markedIndexList.append(i)
        sequenceList.append(sequence)
    for i in range(len(expressionList)):
        if i not in markedIndexList:
            expression=expressionList[i]
            sequence=sequenceList[i]
            myRecord=Record()
            myRecord.sequence=sequence
            myRecord.expression=expression
            recordList.append(myRecord)

    return recordList
    
import methods

methods.deltaEtaFile="/home/lu/Desktop/data_modified/Atha_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Atha_Mutation.csv"
athaRecords=run("atha.fasta", "atha_mRNAExpressionData.csv", "atha")

print athaRecords[0].sequence
        
methods.deltaEtaFile="/home/lu/Desktop/data_modified/Drer_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Drer_Mutation.csv"
run("/home/lu/Desktop/sequences/dmel.fasta", "/home/lu/Desktop/Expression/Dmelangaster/dmel_mrna_expression.csv", "dmel",geneIdType="gene")

methods.deltaEtaFile="/home/lu/Desktop/data_modified/Scer_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Scer_Mutation.csv"
run("/home/lu/Desktop/sequences/s288c.fasta", "/home/lu/Desktop/Expression/Yeast/yeast_mrna_expression.csv", "yeast",geneIdType="locus_tag")

methods.deltaEtaFile="/home/lu/Desktop/data_modified/Xtro_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Xtro_Mutation.csv"
run("/home/lu/Desktop/sequences/xtro.fasta", "/home/lu/Desktop/Expression/x.tropicalis/xtro_mrna_expression.csv", "frog",geneIdType="gene")

methods.deltaEtaFile="/home/lu/Desktop/data_modified/.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Xtro_Mutation.csv"
run("/home/lu/Desktop/sequences/xtro.fasta", "/home/lu/Desktop/Expression/x.tropicalis/xtro_mrna_expression.csv", "frog",geneIdType="gene")