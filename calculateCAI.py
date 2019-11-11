#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:02:15 2019

@author: lu
"""


from CAI import CAI


referenceList=[]


def setUpReference(inputList):
    global referenceList
    referenceList=inputList


def calculateCAIForAGene(sequence):
    global referenceList
    if len(referenceList)==0:
        print ("Call setUpReference function to set up the reference list for calculation CAI")
        exit()
    else: 
        return CAI(sequence,reference=referenceList)
    

def test():
    sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
    setUpReference([sequence,sequence])
    print (calculateCAIForAGene("ATGAAA"))
    
#test()