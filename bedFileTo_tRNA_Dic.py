#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 13:08:55 2019

@author: lu
"""
#this code converts the tRNA bed file into tRNA dictionary availiable for use in codes.
import sys

f= "sacCer3-tRNAs.bed"
try:
    f=sys.argv[1]
except:
    print "no bed file given"
    
    
    

lines=open(f,"r+").readlines()
antiCodonDict=dict() # records the tRNA count by anticodon
for line in lines:
    splitList=line.split("\t")
#    print splitList
    infoHunk=splitList[3]#this is the hunk contaning copies of tRNA information, ex: "tRNA-Pro-TGG-1-1"
#    print infoHunk
    hunkSplit=infoHunk.split("-")

    antiCodon=hunkSplit[2]
    if "NNN" in antiCodon:
        continue
    if antiCodon in antiCodonDict:
        antiCodonDict[antiCodon]+=1
    else:
        antiCodonDict[antiCodon]=1

    
from Bio.Seq import Seq
codontRNADict=dict()
for antiCodon in antiCodonDict:
    cnt=antiCodonDict[antiCodon]
    seq= Seq(antiCodon)  
    codon=str(seq.reverse_complement())
    codontRNADict[codon]=cnt
    
print codontRNADict