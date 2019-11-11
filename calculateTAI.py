#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 08:44:28 2019

@author: lu
"""

from tAI import tAI
import pandas as pd
#import methods

import scipy.stats as stats

#sequence="ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAAAATCACTTTGAATCAGCTGTGGGATATATCTGGTAAATATTTTGATTTGTCTGATAAAAAAGTTAAACAGTTCGTGCTTTCATGCGTGATATTGAAAAAGGACATTGAGGTGTATTGTGATGGTGCTATAACAACTAAAAATGTGACTGATATTATAGGCGACGCTAATCATTCATACTCGGTTGGGATTACTGAGGACAGCCTATGGACATTATTAACGGGATACACAAAAAAGGAGTCAACTATTGGAAATTCTGCATTTGAACTACTTCTCGAAGTTGCCAAATCAGGAGAAAAAGGGATCAATACTATGGATTTGGCGCAGGTAACTGGGCAAGATCCTAGAAGTGTGACTGGACGTATCAAGAAAATAAACCACCTGTTAACAAGTTCACAACTGATTTATAAGGGACACGTCGTGAAGCAATTGAAGCTAAAAAAATTCAGCCATGACGGGGTGGATAGTAATCCCTATATTAATATTAGGGATCATTTAGCAACAATAGTTGAGGTGGTAAAACGATCAAAAAATGGTATTCGCCAGATAATTGATTTAAAGCGTGAATTGAAATTTGACAAAGAGAAAAGACTTTCTAAAGCTTTTATTGCAGCTATTGCATGGTTAGATGAAAAGGAGTACTTAAAGAAAGTGCTTGTAGTATCACCCAAGAATCCTGCCATTAAAATCAGATGTGTAAAATACGTGAAAGATATTCCAGACTCTAAAGGCTCGCCTTCATTTGAGTATGATAGCAATAGCGCGGATGAAGATTCTGTATCAGATAGCAAGGCAGCTTTCGAAGATGAAGACTTAGTCGAAGGTTTAGATAATTTCAATGCGACTGATTTATTACAAAATCAAGGCCTTGTTATGGAAGAGAAAGAGGATGCTGTAAAGAATGAAGTTCTTCTTAATCGATTTTATCCACTTCAAAATCAGACTTATGACATTGCAGATAAGTCTGGCCTTAAAGGAATTTCAACTATGGATGTTGTAAATCGAATTACCGGAAAAGAATTTCAGCGAGCTTTTACCAAATCAAGCGAATATTATTTAGAAAGTGTGGATAAGCAAAAAGAAAATACAGGGGGGTATAGGCTTTTTCGCATATACGATTTTGAGGGAAAGAAGAAGTTTTTTAGGCTGTTCACAGCTCAGAACTTTCAAAAGTTAACAAATGCGGAAGACGAAATATCCGTTCCAAAAGGGTTTGATGAGCTAGGCAAATCTCGTACCGATTTGAAAACTCTCAACGAGGATAATTTCGTCGCACTCAACAACACTGTTAGATTTACAACGGACAGCGATGGACAGGATATATTCTTCTGGCACGGTGAATTAAAAATTCCCCCAAACTCAAAAAAAACTCCGAATAAAAACAAACGGAAGAGGCAGGTTAAAAACAGTACTAATGCTTCTGTTGCAGGAAACATTTCGAATCCCAAAAGGATTAAGCTAGAGCAGCATGTCAGCACTGCACAGGAGCCGAAATCTGCTGAAGATAGTCCAAGTTCAAACGGAGGCACTGTTGTCAAAGGCAAGGTGGTTAACTTCGGCGGCTTTTCTGCCCGCTCTTTGCGTTCACTACAGAGACAGAGAGCCATTTTGAAAGTTATGAATACGATTGGTGGGGTAGCATACCTGAGAGAACAATTTTACGAAAGCGTTTCTAAATATATGGGCTCCACAACGACATTAGATAAAAAGACTGTCCGTGGTGATGTTGATTTGATGGTAGAAAGCGAAAAATTAGGAGCCAGAACAGAGCCTGTATCAGGAAGAAAAATTATTTTTTTGCCCACTGTTGGAGAGGACGCTATCCAAAGGTACATCCTGAAAGAAAAAGATAGTAAAAAAGCAACCTTTACTGATGTTATACATGATACGGAAATATACTTCTTTGACCAAACGGAAAAAAATAGGTTTCACAGAGGAAAGAAATCAGTTGAAAGAATTCGTAAGTTTCAGAACCGCCAAAAGAATGCTAAGATCAAAGCTTCAGATGACGCTATCTCTAAGAAGAGTACGTCGGTCAACGTATCAGATGGAAAGATCAAAAGGAGAGACAAAAAAGTGTCTGCTGGTAGGACAACGGTGGTCGTGGAAAATACTAAAGAAGACAAAACTGTCTATCATGCAGGCACTAAAGATGGTGTTCAGGCTTTAATCAGAGCTGTTGTAGTTACTAAAAGTATTAAAAATGAAATAATGTGGGACAAAATAACAAAATTATTTCCTAATAATTCTTTAGATAACCTAAAAAAGAAATGGACGGCACGGCGAGTAAGAATGGGTCATAGTGGTTGGAGGGCATATGTCGATAAGTGGAAAAAAATGCTCGTTCTAGCCATTAAAAGTGAAAAGATTTCACTGAGGGATGTTGAAGAACTAGATCTTATCAAATTGCTTGATATTTGGACCTCTTTTGATGAAAAGGAAATAAAAAGGCCGCTCTTTCTTTATAAGAACTACGAAGAGAATAGAAAAAAATTTACTCTGGTACGTGATGACACACTTACACATTCTGGCAACGATCTGGCCATGTCTTCTATGATTCAAAGAGAGATCTCTTCTTTAAAAAAAACTTACACTAGAAAGATTTCCGCTTCTACTAAGGACTTATCGAAGAGTCAAAGCGACGATTATATTCGCACAGTGATCCGGTCCATATTAATAGAAAGTCCTTCGACCACTAGAAATGAAATAGAGGCGTTGAAGAACGTTGGAAACGAATCAATAGATAACGTCATCATGGATATGGCTAAGGAAAAGCAAATTTATCTCCATGGCTCAAAACTTGAATGTACTGATACTTTACCAGACATTTTGGAAAATAGAGGAAATTATAAAGATTTTGGTGTAGCTTTTCAGTATAGATGTAAGGTTAATGAATTATTGGAGGCCGGAAACGCTATTGTTATCAATCAAGAGCCGTCCGATATATCCTCTTGGGTTTTAATTGATTTGATTTCGGGAGAGCTATTGAATATGGATGTAATTCCAATGGTGAGAAATGTTCGACCTTTAACGTATACTTCAAGGAGATTTGAAATACGAACATTAACTCCCCCTCTGATTATATATGCCAATTCTCAGACAAAATTGAATACAGCAAGGAAGTCTGCTGTCAAAGTTCCACTGGGCAAACCATTTTCTCGTTTATGGGTGAATGGATCTGGTTCCATTAGGCCAAACATATGGAAGCAGGTAGTTACTATGGTCGTTAACGAAATAATATTTCATCCAGGGATAACATTGAGTAGATTGCAATCTAGGTGTCGTGAAGTACTTTCGCTTCATGAAATATCAGAAATATGCAAATGGCTCCTAGAAAGACAAGTATTAATAACTACTGATTTTGATGGCTATTGGGTCAATCATAATTGGTATTCTATATATGAATCTACATAA"
#sequenceList=[]
#sequenceList.append(sequence)
#sequenceList.append(sequence)




def calculate(sequenceList,species):
    
    dmel_trna={'CTT': 4, 'ATG': 12, 'AAG': 13, 'AAA': 6, 'AAC': 10, 'ATA': 2, 'AGG': 3, 'CCT': 7, 'AGC': 6, 'ACA': 6, 'AGA': 3, 'ATT': 10, 'CTG': 8, 'CTA': 2, 'ACT': 9, 'CAC': 6, 'ACG': 3, 'CCG': 5, 'CAG': 8, 'CAA': 4, 'CGA': 10, 'CCA': 5, 'TCT': 8, 'TGC': 7, 'TGA': 1, 'GGA': 6, 'TGG': 8, 'GGC': 14, 'TAC': 10, 'GAG': 15, 'TCG': 4, 'TTA': 4, 'TTG': 4, 'CGT': 10, 'GAA': 6, 'TCA': 2, 'GCA': 2, 'GTA': 2, 'GCG': 3, 'GTG': 7, 'TTC': 8, 'GTT': 6, 'GCT': 12, 'GAC': 14}
    atha= {'CTT': 12, 'ATG': 25, 'AAG': 18, 'AAA': 14, 'AAC': 16, 'ATA': 5, 'AGG': 8, 'CCT': 16, 'CTC': 1, 'AGC': 13, 'ACA': 8, 'AGA': 9, 'ATT': 17, 'CTG': 3, 'CTA': 10, 'ACT': 10, 'CAC': 10, 'ACG': 5, 'CCG': 5, 'CAG': 9, 'CAA': 8, 'CGA': 6, 'CCA': 45, 'TCT': 37, 'CGG': 4, 'TGC': 16, 'GGG': 5, 'GGA': 13, 'TGG': 14, 'GGC': 23, 'TAC': 76, 'GAG': 13, 'TCG': 4, 'TTA': 6, 'TTG': 10, 'CGT': 9, 'GAA': 14, 'TAA': 1, 'GCA': 10, 'GTA': 7, 'GCG': 7, 'GTG': 8, 'TTC': 16, 'GTT': 15, 'GCT': 16, 'GAC': 27, 'TCC': 1, 'TCA': 9}
    yeast={'ATG': 10, 'AAG': 14, 'AAA': 7, 'AAC': 10, 'ATA': 2, 'AGG': 1, 'CCT': 2, 'CTC': 1, 'AGC': 2, 'ACA': 4, 'AGA': 11, 'ATT': 13, 'CTA': 3, 'ACT': 11, 'CAC': 7, 'ACG': 1, 'CAA': 9, 'CAG': 1, 'CCA': 10, 'TCT': 11, 'CGG': 1, 'TGC': 4, 'GGG': 2, 'GGA': 3, 'TGG': 6, 'GGC': 16, 'TAC': 8, 'GAG': 2, 'TCG': 1, 'TTA': 7, 'TTG': 10, 'CGT': 6, 'GAA': 14, 'TCA': 3, 'GCA': 5, 'GTA': 2, 'GTG': 2, 'TTC': 10, 'GTT': 14, 'GCT': 11, 'GAC': 16}
    frog={'CTT': 40, 'ATG': 205, 'AAG': 107, 'AAA': 110, 'AAC': 91, 'ATA': 43, 'AGG': 38, 'CCT': 63, 'AGC': 91, 'ACA': 90, 'AGA': 143, 'CAT': 1, 'AAT': 1, 'ATT': 76, 'CTG': 60, 'CTA': 25, 'ACT': 158, 'CAC': 50, 'ACG': 31, 'CCG': 29, 'CAG': 38, 'CAA': 26, 'GGT': 1, 'TGT': 3, 'CGA': 44, 'CCA': 48, 'TCT': 41, 'CGG': 6, 'TTT': 1, 'TGC': 88, 'GGG': 154, 'TGA': 4, 'GGA': 71, 'TGG': 41, 'GGC': 110, 'TAC': 137, 'GAG': 68, 'TCG': 29, 'TTA': 50, 'TTG': 30, 'CGT': 40, 'GAA': 70, 'TCA': 29, 'GCA': 73, 'GTA': 32, 'GCG': 31, 'GTG': 49, 'TTC': 62, 'GTT': 72, 'GCT': 61, 'GAC': 80, 'TAA': 1, 'CGC': 2}
    human={'CTT': 13, 'ATG': 20, 'AAG': 21, 'AAA': 19, 'ATC': 3, 'AAC': 35, 'ATA': 5, 'AGG': 8, 'CCT': 11, 'AGC': 8, 'ACA': 6, 'AGA': 6, 'AAT': 2, 'ATT': 16, 'CTG': 9, 'CTA': 4, 'ACT': 10, 'CAC': 10, 'ACG': 6, 'CCG': 4, 'AGT': 1, 'CAG': 22, 'CAA': 9, 'CCC': 1, 'TAT': 1, 'GCG': 5, 'TGT': 1, 'CGA': 6, 'CCA': 9, 'TCT': 11, 'CGG': 4, 'TGC': 31, 'GGG': 10, 'TGA': 3, 'GGA': 11, 'TGG': 7, 'GGC': 15, 'TAC': 15, 'GAG': 12, 'TCG': 4, 'TTA': 5, 'TTG': 7, 'CGT': 7, 'GAA': 16, 'TCA': 4, 'GCA': 10, 'GTA': 5, 'TAA': 1, 'GTG': 18, 'TTC': 17, 'GTT': 11, 'GCT': 30, 'GAC': 19}
    fruitfly={'CTT': 19, 'ATG': 19, 'AAG': 31, 'AAA': 14, 'AAC': 20, 'ATA': 7, 'AGG': 4, 'CCT': 7, 'AGC': 9, 'ACA': 11, 'AGA': 8, 'ATT': 21, 'CTG': 5, 'CTA': 3, 'ACT': 17, 'CAC': 18, 'ACG': 7, 'CCG': 4, 'CAG': 7, 'CAA': 20, 'CGA': 8, 'CCA': 31, 'TCT': 15, 'CGG': 1, 'TGC': 13, 'GGG': 3, 'TGA': 1, 'GGA': 35, 'TGG': 12, 'GGC': 15, 'TAC': 19, 'GAG': 24, 'TCG': 6, 'TTA': 4, 'TTG': 7, 'CGT': 19, 'GAA': 17, 'TCA': 9, 'GCA': 9, 'GTA': 5, 'GCG': 4, 'GTG': 6, 'TTC': 15, 'GTT': 19, 'GCT': 23, 'GAC': 27}
    ecoli={'ACC': 2, 'ATG': 8, 'ACA': 1, 'AAA': 6, 'ATC': 3, 'AAC': 4, 'AGG': 1, 'AGC': 1, 'AGA': 1, 'CTG': 4, 'CTA': 1, 'CTC': 1, 'CAC': 1, 'ACG': 2, 'CCG': 1, 'CCA': 1, 'CAA': 2, 'CCC': 1, 'CAG': 2, 'CGG': 1, 'TGC': 1, 'GGG': 1, 'TGA': 1, 'GGA': 1, 'TGG': 1, 'GGC': 4, 'TAC': 3, 'TTC': 2, 'TCG': 1, 'TTA': 1, 'TTG': 1, 'TCC': 2, 'GAA': 4, 'TCA': 1, 'GCA': 3, 'GTA': 5, 'GCC': 2, 'GTC': 2, 'GAC': 3, 'CGT': 4}
    crei={'CTT': 19, 'ATG': 19, 'AAG': 31, 'AAA': 14, 'AAC': 20, 'ATA': 7, 'AGG': 4, 'CCT': 7, 'AGC': 9, 'ACA': 11, 'AGA': 8, 'ATT': 21, 'CTG': 5, 'CTA': 3, 'ACT': 17, 'CAC': 18, 'ACG': 7, 'CCG': 4, 'CAG': 7, 'CAA': 20, 'CGA': 8, 'CCA': 31, 'TCT': 15, 'CGG': 1, 'TGC': 13, 'GGG': 3, 'TGA': 1, 'GGA': 35, 'TGG': 12, 'GGC': 15, 'TAC': 19, 'GAG': 24, 'TCG': 6, 'TTA': 4, 'TTG': 7, 'CGT': 19, 'GAA': 17, 'TCA': 9, 'GCA': 9, 'GTA': 5, 'GCG': 4, 'GTG': 6, 'TTC': 15, 'GTT': 19, 'GCT': 23, 'GAC': 27}
    atha_abundances = pd.Series(atha)
    atha_tai_abundance = tAI(atha_abundances)
    
    dmel_abundances = pd.Series(dmel_trna)
    dmel_tai_abundance = tAI(dmel_abundances)
    
    yeast_abundances = pd.Series(yeast)
    yeast_tai_abundance = tAI(yeast_abundances)
    
    frog_abundances = pd.Series(frog)
    frog_tai_abundance = tAI(frog_abundances)
    
    human_abundances = pd.Series(human)
    human_tai_abundance = tAI(human_abundances)
    
    fruitfly_abundances = pd.Series(fruitfly)
    fruitfly_tai_abundance = tAI(fruitfly_abundances)
    
    ecoli_abundances = pd.Series(ecoli)
    ecoli_tai_abundance = tAI(ecoli_abundances)
    
    crei_abundance=pd.Series(crei)
    crei_tai_abundance=tAI(crei_abundance)
    
    tAIList=[]
    for sequence in sequenceList:
        if "atha" in species:
            tAIList.append(atha_tai_abundance.calc(sequence))
        elif"dmel" in species:
            tAIList.append(dmel_tai_abundance.calc(sequence))
        elif"yeast" in species:
            tAIList.append(yeast_tai_abundance.calc(sequence))
        elif"frog" in species:
            tAIList.append(frog_tai_abundance.calc(sequence))
        elif"human" in species:
            tAIList.append(human_tai_abundance.calc(sequence))                        
        elif"fruitfly" in species:
            tAIList.append(fruitfly_tai_abundance.calc(sequence))
        elif "ecoli" in species:
            tAIList.append(ecoli_tai_abundance.calc(sequence))
        elif "crei" in species:
            tAIList.append(crei_tai_abundance.calc(sequence))
        else: 
            print ("Species %s not supported"%(species))
            exit(0)
    return tAIList
dmel_trna={'CTT': 4, 'ATG': 12, 'AAG': 13, 'AAA': 6, 'AAC': 10, 'ATA': 2, 'AGG': 3, 'CCT': 7, 'AGC': 6, 'ACA': 6, 'AGA': 3, 'ATT': 10, 'CTG': 8, 'CTA': 2, 'ACT': 9, 'CAC': 6, 'ACG': 3, 'CCG': 5, 'CAG': 8, 'CAA': 4, 'CGA': 10, 'CCA': 5, 'TCT': 8, 'TGC': 7, 'TGA': 1, 'GGA': 6, 'TGG': 8, 'GGC': 14, 'TAC': 10, 'GAG': 15, 'TCG': 4, 'TTA': 4, 'TTG': 4, 'CGT': 10, 'GAA': 6, 'TCA': 2, 'GCA': 2, 'GTA': 2, 'GCG': 3, 'GTG': 7, 'TTC': 8, 'GTT': 6, 'GCT': 12, 'GAC': 14}
atha= {'CTT': 12, 'ATG': 25, 'AAG': 18, 'AAA': 14, 'AAC': 16, 'ATA': 5, 'AGG': 8, 'CCT': 16, 'CTC': 1, 'AGC': 13, 'ACA': 8, 'AGA': 9, 'ATT': 17, 'CTG': 3, 'CTA': 10, 'ACT': 10, 'CAC': 10, 'ACG': 5, 'CCG': 5, 'CAG': 9, 'CAA': 8, 'CGA': 6, 'CCA': 45, 'TCT': 37, 'CGG': 4, 'TGC': 16, 'GGG': 5, 'GGA': 13, 'TGG': 14, 'GGC': 23, 'TAC': 76, 'GAG': 13, 'TCG': 4, 'TTA': 6, 'TTG': 10, 'CGT': 9, 'GAA': 14, 'TAA': 1, 'GCA': 10, 'GTA': 7, 'GCG': 7, 'GTG': 8, 'TTC': 16, 'GTT': 15, 'GCT': 16, 'GAC': 27, 'TCC': 1, 'TCA': 9}
yeast={'ATG': 10, 'AAG': 14, 'AAA': 7, 'AAC': 10, 'ATA': 2, 'AGG': 1, 'CCT': 2, 'CTC': 1, 'AGC': 2, 'ACA': 4, 'AGA': 11, 'ATT': 13, 'CTA': 3, 'ACT': 11, 'CAC': 7, 'ACG': 1, 'CAA': 9, 'CAG': 1, 'CCA': 10, 'TCT': 11, 'CGG': 1, 'TGC': 4, 'GGG': 2, 'GGA': 3, 'TGG': 6, 'GGC': 16, 'TAC': 8, 'GAG': 2, 'TCG': 1, 'TTA': 7, 'TTG': 10, 'CGT': 6, 'GAA': 14, 'TCA': 3, 'GCA': 5, 'GTA': 2, 'GTG': 2, 'TTC': 10, 'GTT': 14, 'GCT': 11, 'GAC': 16}
frog={'CTT': 40, 'ATG': 205, 'AAG': 107, 'AAA': 110, 'AAC': 91, 'ATA': 43, 'AGG': 38, 'CCT': 63, 'AGC': 91, 'ACA': 90, 'AGA': 143, 'CAT': 1, 'AAT': 1, 'ATT': 76, 'CTG': 60, 'CTA': 25, 'ACT': 158, 'CAC': 50, 'ACG': 31, 'CCG': 29, 'CAG': 38, 'CAA': 26, 'GGT': 1, 'TGT': 3, 'CGA': 44, 'CCA': 48, 'TCT': 41, 'CGG': 6, 'TTT': 1, 'TGC': 88, 'GGG': 154, 'TGA': 4, 'GGA': 71, 'TGG': 41, 'GGC': 110, 'TAC': 137, 'GAG': 68, 'TCG': 29, 'TTA': 50, 'TTG': 30, 'CGT': 40, 'GAA': 70, 'TCA': 29, 'GCA': 73, 'GTA': 32, 'GCG': 31, 'GTG': 49, 'TTC': 62, 'GTT': 72, 'GCT': 61, 'GAC': 80, 'TAA': 1, 'CGC': 2}
human={'CTT': 13, 'ATG': 20, 'AAG': 21, 'AAA': 19, 'ATC': 3, 'AAC': 35, 'ATA': 5, 'AGG': 8, 'CCT': 11, 'AGC': 8, 'ACA': 6, 'AGA': 6, 'AAT': 2, 'ATT': 16, 'CTG': 9, 'CTA': 4, 'ACT': 10, 'CAC': 10, 'ACG': 6, 'CCG': 4, 'AGT': 1, 'CAG': 22, 'CAA': 9, 'CCC': 1, 'TAT': 1, 'GCG': 5, 'TGT': 1, 'CGA': 6, 'CCA': 9, 'TCT': 11, 'CGG': 4, 'TGC': 31, 'GGG': 10, 'TGA': 3, 'GGA': 11, 'TGG': 7, 'GGC': 15, 'TAC': 15, 'GAG': 12, 'TCG': 4, 'TTA': 5, 'TTG': 7, 'CGT': 7, 'GAA': 16, 'TCA': 4, 'GCA': 10, 'GTA': 5, 'TAA': 1, 'GTG': 18, 'TTC': 17, 'GTT': 11, 'GCT': 30, 'GAC': 19}
fruitfly={'CTT': 19, 'ATG': 19, 'AAG': 31, 'AAA': 14, 'AAC': 20, 'ATA': 7, 'AGG': 4, 'CCT': 7, 'AGC': 9, 'ACA': 11, 'AGA': 8, 'ATT': 21, 'CTG': 5, 'CTA': 3, 'ACT': 17, 'CAC': 18, 'ACG': 7, 'CCG': 4, 'CAG': 7, 'CAA': 20, 'CGA': 8, 'CCA': 31, 'TCT': 15, 'CGG': 1, 'TGC': 13, 'GGG': 3, 'TGA': 1, 'GGA': 35, 'TGG': 12, 'GGC': 15, 'TAC': 19, 'GAG': 24, 'TCG': 6, 'TTA': 4, 'TTG': 7, 'CGT': 19, 'GAA': 17, 'TCA': 9, 'GCA': 9, 'GTA': 5, 'GCG': 4, 'GTG': 6, 'TTC': 15, 'GTT': 19, 'GCT': 23, 'GAC': 27}
ecoli={'ACC': 2, 'ATG': 8, 'ACA': 1, 'AAA': 6, 'ATC': 3, 'AAC': 4, 'AGG': 1, 'AGC': 1, 'AGA': 1, 'CTG': 4, 'CTA': 1, 'CTC': 1, 'CAC': 1, 'ACG': 2, 'CCG': 1, 'CCA': 1, 'CAA': 2, 'CCC': 1, 'CAG': 2, 'CGG': 1, 'TGC': 1, 'GGG': 1, 'TGA': 1, 'GGA': 1, 'TGG': 1, 'GGC': 4, 'TAC': 3, 'TTC': 2, 'TCG': 1, 'TTA': 1, 'TTG': 1, 'TCC': 2, 'GAA': 4, 'TCA': 1, 'GCA': 3, 'GTA': 5, 'GCC': 2, 'GTC': 2, 'GAC': 3, 'CGT': 4}
crei={'CTT': 19, 'ATG': 19, 'AAG': 31, 'AAA': 14, 'AAC': 20, 'ATA': 7, 'AGG': 4, 'CCT': 7, 'AGC': 9, 'ACA': 11, 'AGA': 8, 'ATT': 21, 'CTG': 5, 'CTA': 3, 'ACT': 17, 'CAC': 18, 'ACG': 7, 'CCG': 4, 'CAG': 7, 'CAA': 20, 'CGA': 8, 'CCA': 31, 'TCT': 15, 'CGG': 1, 'TGC': 13, 'GGG': 3, 'TGA': 1, 'GGA': 35, 'TGG': 12, 'GGC': 15, 'TAC': 19, 'GAG': 24, 'TCG': 6, 'TTA': 4, 'TTG': 7, 'CGT': 19, 'GAA': 17, 'TCA': 9, 'GCA': 9, 'GTA': 5, 'GCG': 4, 'GTG': 6, 'TTC': 15, 'GTT': 19, 'GCT': 23, 'GAC': 27}

atha_abundances = pd.Series(atha)
atha_tai_abundance = tAI(atha_abundances)

dmel_abundances = pd.Series(dmel_trna)
dmel_tai_abundance = tAI(dmel_abundances)

yeast_abundances = pd.Series(yeast)
yeast_tai_abundance = tAI(yeast_abundances)

frog_abundances = pd.Series(frog)
frog_tai_abundance = tAI(frog_abundances)

human_abundances = pd.Series(human)
human_tai_abundance = tAI(human_abundances)

fruitfly_abundances = pd.Series(fruitfly)
fruitfly_tai_abundance = tAI(fruitfly_abundances)

ecoli_abundances = pd.Series(ecoli)
ecoli_tai_abundance = tAI(ecoli_abundances)

crei_abundance=pd.Series(crei)
crei_tai_abundance=tAI(crei_abundance)

def calculateOneGene(sequence,species):

    
    global atha_tai_abundance,dmel_tai_abundance,yeast_tai_abundance,frog_tai_abundance,human_tai_abundance,fruitfly_tai_abundance,ecoli_tai_abundance
    

    if "atha" in species:
        return atha_tai_abundance.calc(sequence)
    elif"dmel" in species:
        return dmel_tai_abundance.calc(sequence)
    elif"yeast" in species:
        return yeast_tai_abundance.calc(sequence)
    elif"frog" in species:
        return frog_tai_abundance.calc(sequence)
    elif"human" in species:
        return human_tai_abundance.calc(sequence)                        
    elif"fruitfly" in species:
        return fruitfly_tai_abundance.calc(sequence)
    elif "ecoli" in species:
        return ecoli_tai_abundance.calc(sequence)
    elif "crei" in species:
        return (crei_tai_abundance.calc(sequence))
    else: 
        print ("Species %s not supported")
        exit(0)

            
            
