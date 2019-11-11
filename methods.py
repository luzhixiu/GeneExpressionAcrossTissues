import math
import numpy as np
from scipy import optimize
from scipy.optimize import minimize_scalar
from sklearn.preprocessing import minmax_scale
from statistics import mean





deltaEtaFile="/home/lu/Code/selection.csv"
deltaMFile="/home/lu/Code/mutation.csv"


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





##read in a fasta file, return the condon list from the first gene
#def readSequence(fname):
#    f=open(fname)
#    lines=f.readlines()
#    sequence= lines[1]
#    sequence="ATGGGTGTTGAACAAATCTTAAAGAGAAAGACCGGTGTCATCGTTGGTGAAGATGTCCACAACTTATTCACTTACGCTAAGGAACACAAGTTCGCTATTCCAGCTATTAACGTCACCTCTTCTTCTACTGCCGTCGCTGCTTTAGAAGCTGCTAGAGACAGCAAGTCCCCAATCATTTTGCAAACCTCTAACGGTGGTGCTGCTTACTTCGCTGGTAAGGGTATCTCTAACGAAGGTCAAAATGCTTCCATCAAGGGTGCTATTGCCGCTGCCCACTACATCAGATCCATTGCTCCAGCTTACGGTATCCCAGTTGTCTTACACTCTGACCACTGTGCCAAGAAGTTGTTGCCATGGTTCGATGGTATGTTGGAAGCTGATGAAGCTTACTTCAAGGAACACGGTGAACCATTATTCTCCTCCCACATGTTGGATTTGTCTGAAGAAACCGATGAAGAAAACATCTCTACTTGTGTCAAGTACTTCAAGAGAATGGCCGCTATGGACCAATGGTTAGAAATGGAAATCGGTATTACCGGTGGTGAAGAAGATGGTGTTAACAACGAAAACGCTGACAAGGAAGACTTGTACACCAAGCCAGAACAAGTTTACAACGTCTACAAGGCTTTGCACCCAATCTCTCCAAACTTCTCCATTGCTGCTGCTTTCGGTAACTGTCACGGTTTGTACGCTGGTGACATCGCTTTGAGACCAGAAATCTTGGCTGAACACCAAAAGTACACCAGAGAACAAGTTGGTTGCAAGGAAGAAAAGCCATTGTTCTTGGTCTTCCACGGTGGTTCCGGTTCTACTGTCCAAGAATTCCACACTGGTATTGACAACGGTGTTGTCAAGGTCAACTTGGACACTGACTGTCAATACGCTTACTTGACTGGTATCAGAGACTACGTCTTGAACAAGAAGGACTACATAATGTCCCCAGTCGGTAACCCAGAAGGTCCAGAAAAGCCAAACAAGAAGTTCTTCGACCCAAGAGTCTGGGTTAGAGAAGGTGAAAAGACCATGGGTGCTAAGATCACCAAGTCTTTGGAAACTTTCCGTACCACTAACACTTTATAA"
#    startCodon="ATG"
#    stopCodonList=["TAG","TAA","TGA"]
#    codonList=[]
#    i=0
#    while(i<len(sequence)):
#        codon=sequence[i:i+3]
#        if len(codon)==3:
#            codonList.append(codon)
#        i+=3
#    actualCodonList=[]
#    started=False
#    for codon in codonList:
#        if codon in stopCodonList:
#            break
#        if started:
#            actualCodonList.append(codon)
#        if codon==startCodon:
#            started=True
#    codonList=actualCodonList
#   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
#    return codonList

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

#this method removes sequence that cant be handled by mDict, etaDict and "TGG",
def parseSequence(sequence):
    i=0
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    parsedSequence=""
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            if (codon in mDict and codon in etaDict) or (codon in startCodon) or (codon in stopCodonList) and ("TGG" not in codon):
                parsedSequence+=codon
        i+=3
    return parsedSequence
    

def cutSequence(seq):
    sequence=seq
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
    started=FalsewindowSize
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


def roundList(lst):
    roundedList=[]
    decimalPlaces=6
    for i in lst:
        roundedList.append(round(i,decimalPlaces))
    return roundedList
    



# use method 1 from micheals method, simply use rate/pausing time, return the data points for plotting 
# for option, zero returns the rate and 1 returns the pausing time result.
def method1(codonList,option):
    popData=open("Pop.et.al_codon.specific.translation.rates.csv")
    lines=popData.readlines()
    popRateDict=dict()
    for line in lines[1:]:
        splitList=line.split(",")
        if len(splitList[0])>0 and len(splitList[1])>0:
            popRateDict[splitList[0]]=splitList[1]
    
    rateList=[]    
    pausingTimeList=[]    
    for codon in codonList:
            rateList.append(popRateDict[codon])
            pausingTimeList.append(1/float(popRateDict[codon]))
    if "0" in option:
        return roundList(rateList)
    if "1" in option:
        return roundList(pausingTimeList)
        
# use method 2 from micheals method, considers the difference between pausing time, return the data points for plotting 
def method2(codonList):
    global etaDict
    etaList=[]
    for codon in codonList:
        if codon in etaDict:
            eta=etaDict[codon]
        else:
            eta=0
        etaList.append(float(eta))
    return roundList(etaList)
    
    
def method3(codonList):
    global mDict
    global etaDict
    etaList=[]
    for codon in codonList:
        if codon in etaDict:
            eta=etaDict[codon]
        else:
            eta=0
        etaList.append(eta)

    mList=[]
    for codon in codonList:
        if codon in mDict:
            m=mDict[codon]
        else:
            m=0
        mList.append(m)
    
    global synonymonDict
    expectedEtaList=[]
    for codon in codonList:
        if codon in mDict:
            mValue=float(mDict[codon])*(-1.0)
#           take exp of the result            
            mValue= np.exp(mValue)
            codonSynoList=synonymonDict[codon]
            mValueSum=0
            etaValueSumJ=0# added this weighting factor accroding to the write up
            for syno in codonSynoList:
                if syno in mDict:
                    mValueSum+=np.exp(float(mDict[syno])*(-1.0))
                if syno in etaDict:
                    etaValueSumJ+=float(etaDict[syno])
            expectedEta=mValue/mValueSum*etaValueSumJ
            expectedEtaList.append(float(expectedEta))
            
        else:       
            expectedEtaList.append(1)
    scoreList=[]            
    for i in range(0,len(codonList)):
        score=float(etaList[i])-float(expectedEtaList[i])
        scoreList.append(score)
    return roundList(scoreList)


    
    


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
    
0
    
#============================================================================== MinMax
#The function for calculating %MinMax
def calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize):    
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
        freqDict = freqDict
        aaFreqDict = aaFreqDict
        windowSize = windowSize
        mapDict = mapDict
        codonSeq = sequence
        minMaxValues = [] #list to be returned of the %MinMax values
        
        for i in range(int(windowSize/2)): #%MinMax is undefined for the first and last (windowSize/2) condons
            minMaxValues.append(0)
        
        #Using the specified sliding window size (windowSize/2 - 1 on either side of the central codon), min/max is calculated
        for i in range(len(codonSeq)-windowSize+1):
            window = codonSeq[i:i+windowSize] #list of the codons in the current window
            Actual = 0.0     #average of the actual codon frequencies
            Max = 0.0        #average of the min codon frequencies
            Min = 0.0        #average of the max codon frequencies
            Avg = 0.0        #average of the averages of all the frequencies associated with each amino acid
            
            #Sum the frequencies
            for codon in window:
                frequencies = aaFreqDict[mapDict[codon]] #list of all frequencies associated with the amino acid this codon encodes
    
                Actual += freqDict[codon]
                Max += max(frequencies)
                Min += min(frequencies)
                Avg += sum(frequencies)/len(frequencies)
    
            #Divide by the window size to get the averages
            Actual = Actual/windowSize
            Max = Max/windowSize
            Min = Min/windowSize
            Avg = Avg/windowSize
            
            if Max-Avg==0:
                percentMax=0
            else:
                percentMax = ((Actual-Avg)/(Max-Avg))*100
            if Avg-Min==0:
                percentMin=0
            else:
                percentMin = ((Avg-Actual)/(Avg-Min))*100
    
            if(percentMax >= 0):
                minMaxValues.append(round(percentMax,2))
            else:
                minMaxValues.append(round(-percentMin,2))
            
    
        #fills in values for codons where window size makes min/max unable to be calculated
        for i in range(int(windowSize/2)):
            minMaxValues.append(0)
        return minMaxValues


def minMax(sequence):
 #   yeastFreq = {'TCA': 0.03991864708764081, 'AAT': 0.25318475701169, 'TGG': 1.0, 'GAT': 0.5201770221162414, 'GAA': 0.8862390592922943, 'TTC': 0.7159360267466637, 'CCG': 0.004733235294022094, 'ACT': 0.5109177475686215, 'GGG': 0.0094131856365968, 'ACG': 0.012460431925551328, 'AGA': 0.7553924587663011, 'TTG': 0.6186602618936331, 'GTC': 0.37106987401182084, 'GCA': 0.061810833121510644, 'TGA': 0.3333333333333333, 'CGT': 0.2117403300720323, 'CAC': 0.5779805928157615, 'CTC': 0.0034381486449209316, 'CGA': 2.429357848105861e-05, 'GCT': 0.6307564359126062, 'ATC': 0.4587067090086894, 'ATA': 0.006829669948812382, 'TTT': 0.2840639732533363, 'TAA': 0.3333333333333333, 'GTG': 0.05014456069837953, 'GCC': 0.2976864674415778, 'GAG': 0.1137609407077057, 'CAT': 0.4220194071842384, 'AAG': 0.7168657641228129, 'AAA': 0.28313423587718717, 'GCG': 0.009746263524305406, 'TCC': 0.18734579350163716, 'GGC': 0.05305331964855917, 'TCT': 0.29262991895195334, 'CCT': 0.197525039410729, 'GTA': 0.020064008635860053, 'AGG': 0.026750929719980346, 'CCA': 0.7716216780595091, 'TAT': 0.2518279573102916, 'ACC': 0.39591756487343627, 'TCG': 0.010858675400038171, 'ATG': 1.0, 'TTA': 0.227024141795565, 'TGC': 0.09136336592913656, 'GTT': 0.5587215566539396, 'CTT': 0.03090535577924939, 'CAG': 0.06548567125740762, 'CCC': 0.02612004723573985, 'ATT': 0.5344636210424982, 'ACA': 0.08070425563239092, 'AAC': 0.74681524298831, 'GGT': 0.9277254784982976, 'AGC': 0.17661704610677711, 'CGG': 2.494122293434685e-05, 'TAG': 0.3333333333333333, 'CGC': 0.0060670466402708454, 'AGT': 0.29262991895195334, 'CTA': 0.0948250514888338, 'CAA': 0.9345143287425924, 'CTG': 0.025147040397797728, 'GGA': 0.009808016216546523, 'TGT': 0.9086366340708634, 'TAC': 0.7481720426897084, 'GAC': 0.4798229778837586}
  #  ecoliFreq = {'TCA': 0.015895609098584198, 'AAT': 0.052646578348652565, 'TGG': 1.0, 'GAT': 0.346006298147288, 'GAA': 0.8167429068458956, 'TTC': 0.8849166752601701, 'CCG': 0.8374133304795861, 'ACT': 0.34154535987431217, 'GGG': 0.015020413126162143, 'ACG': 0.057816421743531035, 'AGA': 7.077535152530984e-05, 'TTG': 0.016363607002198043, 'GTC': 0.09206888115542544, 'GCA': 0.23200897410863494, 'TGA': 0.3333333333333333, 'CGT': 0.7344771136702194, 'CAC': 0.8593343034461259, 'CTC': 0.05441045805693406, 'CGA': 0.0002787741465201721, 'GCT': 0.339918072915499, 'ATC': 0.777485634469149, 'ATA': 0.00017712105716602427, 'TTT': 0.11508332473982981, 'TAA': 0.3333333333333333, 'GTG': 0.23285420120111727, 'GCC': 0.13197587234590208, 'GAG': 0.18325709315410435, 'CAT': 0.14066569655387415, 'AAG': 0.180102682357364, 'AAA': 0.8198973176426361, 'GCG': 0.2960970806299639, 'TCC': 0.26120555798352335, 'GGC': 0.41126982801100653, 'TCT': 0.33629037499154135, 'CCT': 0.05230295955812376, 'GTA': 0.21801911337924054, 'AGG': 1.4110863201276576e-05, 'CCA': 0.10863553736253913, 'TAT': 0.22118859861578222, 'ACC': 0.5878757461439766, 'TCG': 0.03017145622890794, 'ATG': 1.0, 'TTA': 0.006031627625274781, 'TGC': 0.7135656339457533, 'GTT': 0.4570578042642167, 'CTT': 0.023623881846505348, 'CAG': 0.8967245503451191, 'CCC': 0.0016481725997509686, 'ATT': 0.22233724447368503, 'ACA': 0.012762472238180113, 'AAC': 0.9473534216513474, 'GGT': 0.5710345303394535, 'AGC': 0.33629037499154135, 'CGG': 0.00024960802394793387, 'TAG': 0.3333333333333333, 'CGC': 0.2649096179445859, 'AGT': 0.020146626705901727, 'CTA': 0.0015370230908821575, 'CAA': 0.10327544965488093, 'CTG': 0.8980334023782056, 'GGA': 0.0026752285233779633, 'TGT': 0.28643436605424677, 'TAC': 0.7788114013842178, 'GAC': 0.653993701852712}
    
    #This script was written in Python 3 and utilizes the statistics, random, pandas, and matplotlib.pyplot modules
    #Imports happen throughout the script to minimize long pause times and to solve errors that occur when imports happen at the beginning

    freqDict = dict() #dictionary mapping codons to their frequencies
    mapDict = dict() #dictionary mapping codons to amino acid
    aaFreqDict = dict() #dictionary mapping each amino acid to a list of the frequencies of possible codons for that amino acid
    aaMapDict = dict() #dictionary from amino acid to list of codons with frequencies for it (for RRTs)
    
    
    
    mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
               'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
               'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
               'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
               'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
               'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
               'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
               'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
               'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
               'GAC': 'D'}
    
    aaDict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
              'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
              'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
              'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
              'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
              'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
              'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}
    
    
    #The below dictionary contains the usage files for a handful of selected species. The files are taken from HIVE-CUT
    speciesDict = {
            'Escherichia_coli' : {'TTT': 22.38, 'TCT': 8.61, 'TAT': 16.36, 'TGT': 5.19, 'TTC': 16.21, 
                                  'TCC': 8.81, 'TAC': 12.15, 'TGC': 6.34, 'TTA': 13.83, 'TCA': 7.57, 
                                  'TAA': 2.03, 'TGA': 1.04, 'TTG': 13.37, 'TCG': 8.79, 'TAG': 0.25, 
                                  'TGG': 15.21, 'CTT': 11.44, 'CCT': 7.22, 'CAT': 12.84, 'CGT': 20.7,
                                  'CTC': 10.92, 'CCC': 5.56, 'CAC': 9.44, 'CGC': 21.48, 'CTA': 3.93,
                                  'CCA': 8.44, 'CAA': 15.1, 'CGA': 3.67, 'CTG': 52.1, 'CCG': 22.65, 
                                  'CAG': 29.21, 'CGG': 5.72, 'ATT': 30.21, 'ACT': 9.02, 'AAT': 18.26,
                                  'AGT': 9.08, 'ATC': 24.6, 'ACC': 22.88, 'AAC': 21.47, 'AGC': 15.89,
                                  'ATA': 4.88, 'ACA': 7.63, 'AAA': 33.94, 'AGA': 2.43, 'ATG': 27.59, 
                                  'ACG': 14.47, 'AAG': 10.7, 'AGG': 1.48, 'GTT': 18.39, 'GCT': 15.54,
                                  'GAT': 32.43, 'GGT': 24.45, 'GTC': 15.07, 'GCC': 25.45, 'GAC': 19.14,
                                  'GGC': 28.65, 'GTA': 10.97, 'GCA': 20.61, 'GAA': 39.55, 'GGA': 8.44, 
                                  'GTG': 25.9, 'GCG': 32.79, 'GAG': 18.24, 'GGG': 11.29},
    
    
            'Homo_sapien' : {'TTT': 21.42, 'TCT': 16.96, 'TAT': 17.11, 'TGT': 10.99, 'TTC': 23.05, 
                             'TCC': 10.61, 'TAC': 13.48, 'TGC': 8.76, 'TTA': 9.31, 'TCA': 21.14, 
                             'TAA': 1.07, 'TGA': 0.72, 'TTG': 19.58, 'TCG': 12.81, 'TAG': 0.39, 
                             'TGG': 10.68, 'CTT': 20.92, 'CCT': 9.16, 'CAT': 14.3, 'CGT': 11.38, 
                             'CTC': 14.54, 'CCC': 4.46, 'CAC': 9.07, 'CGC': 5.03, 'CTA': 7.69, 
                             'CCA': 27.25, 'CAA': 27.88, 'CGA': 12.5, 'CTG': 12.05, 'CCG': 10.28, 
                             'CAG': 14.88, 'CGG': 4.84, 'ATT': 31.51, 'ACT': 19.3, 'AAT': 30.06, 
                             'AGT': 12.28, 'ATC': 18.5, 'ACC': 10.32, 'AAC': 17.93, 'AGC': 8.32, 
                             'ATA': 9.03, 'ACA': 20.55, 'AAA': 36.4, 'AGA': 15.25, 'ATG': 26.09, 
                             'ACG': 9.22, 'AAG': 25.58, 'AGG': 3.8, 'GTT': 24.14, 'GCT': 22.94, 
                             'GAT': 36.73, 'GGT': 11.12, 'GTC': 13.58, 'GCC': 12.74, 'GAC': 17.27, 
                             'GGC': 6.78, 'GTA': 9.78, 'GCA': 20.35, 'GAA': 41.61, 'GGA': 32.03, 
                             'GTG': 14.53, 'GCG': 8.56, 'GAG': 25.1, 'GGG': 4.32},
        
    
            'Mus_musculus' : {'TTT': 15.94, 'TCT': 17.39, 'TAT': 11.15, 'TGT': 10.68, 'TTC': 18.81,
                              'TCC': 18.32, 'TAC': 14.42, 'TGC': 10.95, 'TTA': 7.29, 'TCA': 13.31, 
                              'TAA': 0.39, 'TGA': 0.76, 'TTG': 13.27, 'TCG': 4.29, 'TAG': 0.35, 
                              'TGG': 11.44, 'CTT': 13.45, 'CCT': 20.06, 'CAT': 11.23, 'CGT': 4.64, 
                              'CTC': 18.83, 'CCC': 18.34, 'CAC': 15.23, 'CGC': 8.6, 'CTA': 8.04, 
                              'CCA': 19.05, 'CAA': 13.11, 'CGA': 6.9, 'CTG': 37.31, 'CCG': 6.11, 
                              'CAG': 36.71, 'CGG': 10.46, 'ATT': 14.66, 'ACT': 13.92, 'AAT': 15.9, 
                              'AGT': 14.11, 'ATC': 20.33, 'ACC': 18.16, 'AAC': 19.75, 'AGC': 20.77,
                              'ATA': 7.28, 'ACA': 16.67, 'AAA': 23.84, 'AGA': 13.03, 'ATG': 21.7, 
                              'ACG': 5.73, 'AAG': 33.95, 'AGG': 12.94, 'GTT': 10.81, 'GCT': 20.19,
                              'GAT': 22.33, 'GGT': 11.09, 'GTC': 14.52, 'GCC': 25.16, 'GAC': 26.3, 
                              'GGC': 19.81, 'GTA': 7.48, 'GCA': 16.8, 'GAA': 30.33, 'GGA': 16.77, 
                              'GTG': 26.58, 'GCG': 5.86, 'GAG': 41.48, 'GGG': 14.91},
    
            'Caenorhabditis_elegans' : {'TTT': 17.06, 'TCT': 16.58, 'TAT': 12.04, 'TGT': 10.54, 'TTC': 17.87, 
                                        'TCC': 17.44, 'TAC': 13.7, 'TGC': 11.15, 'TTA': 8.55, 'TCA': 13.89,
                                        'TAA': 0.46, 'TGA': 0.83, 'TTG': 13.3, 'TCG': 4.18, 'TAG': 0.36, 
                                        'TGG': 11.77, 'CTT': 13.95, 'CCT': 18.88, 'CAT': 11.74, 'CGT': 4.54,
                                        'CTC': 18.06, 'CCC': 19.19, 'CAC': 14.76, 'CGC': 9.06, 'CTA': 7.39, 
                                        'CCA': 18.45, 'CAA': 13.83, 'CGA': 6.36, 'CTG': 36.75, 'CCG': 6.36, 
                                        'CAG': 35.3, 'CGG': 10.88, 'ATT': 16.36, 'ACT': 14.12, 'AAT': 18.16, 
                                        'AGT': 13.72, 'ATC': 18.97, 'ACC': 17.95, 'AAC': 18.36, 'AGC': 19.74, 
                                        'ATA': 7.98, 'ACA': 16.33, 'AAA': 27.15, 'AGA': 13.09, 'ATG': 21.4, 
                                        'ACG': 5.72, 'AAG': 31.89, 'AGG': 12.15, 'GTT': 11.59, 'GCT': 18.77, 
                                        'GAT': 23.68, 'GGT': 10.75, 'GTC': 13.58, 'GCC': 26.18, 'GAC': 24.49, 
                                        'GGC': 20.23, 'GTA': 7.56, 'GCA': 16.89, 'GAA': 33.04, 'GGA': 17.02, 
                                        'GTG': 26.24, 'GCG': 6.26, 'GAG': 39.88, 'GGG': 15.53},
        
        
        'Saccharomyces_cerevisiae' : {'TTT': 26.18, 'TCT': 23.35, 'TAT': 19.05, 'TGT': 7.82, 'TTC': 17.88,
                                      'TCC': 14.07, 'TAC': 14.6, 'TGC': 4.75, 'TTA': 26.33, 'TCA': 19.05, 
                                      'TAA': 0.95, 'TGA': 0.6, 'TTG': 26.5, 'TCG': 8.71, 'TAG': 0.46, 
                                      'TGG': 10.35, 'CTT': 12.27, 'CCT': 13.57, 'CAT': 13.89, 'CGT': 6.26, 
                                      'CTC': 5.52, 'CCC': 6.91, 'CAC': 7.74, 'CGC': 2.63, 'CTA': 13.52, 
                                      'CCA': 17.81, 'CAA': 27.1, 'CGA': 3.1, 'CTG': 10.65, 'CCG': 5.42, 
                                      'CAG': 12.42, 'CGG': 1.82, 'ATT': 30.1, 'ACT': 20.24, 'AAT': 36.61, 
                                      'AGT': 14.6, 'ATC': 16.99, 'ACC': 12.48, 'AAC': 24.8, 'AGC': 9.96,
                                      'ATA': 18.29, 'ACA': 18.18, 'AAA': 42.83, 'AGA': 21.05, 'ATG': 20.68, 
                                      'ACG': 8.15, 'AAG': 30.52, 'AGG': 9.45, 'GTT': 21.47, 'GCT': 20.28,
                                      'GAT': 38.09, 'GGT': 22.59, 'GTC': 11.23, 'GCC': 12.14, 'GAC': 20.39,
                                      'GGC': 9.78, 'GTA': 12.07, 'GCA': 16.26, 'GAA': 45.81, 'GGA': 11.19,
                                      'GTG': 10.72, 'GCG': 6.17, 'GAG': 19.55, 'GGG': 6.06}
    }
    
    
    #Below are the lists of the species contained in the included usage files
    speciesList = ['Escherichia coli', 'Saccharomyces cerevisiae']
    speciesList2 = ['escherichia coli', 'saccharomyces cerevisiae', '1','2']
    
    
    #inputDict allows numbers or names to be entered for ease of use
    inputDict = {'escherichia coli':'Escherichia_coli', '1': 'Escherichia_coli', 'saccharomyces cerevisiae':'Saccharomyces_cerevisiae', '2':'Saccharomyces_cerevisiae'}
    
    #user_filename = input('Name your file:')#Allows user to name their output file - will be saved to the working directory
    user_filename="test";
    
    
    #windowSize = int(input('Sliding Window Length (odd # only): ')) #size of the sliding window for %MinMax
    windowSize=17
    while windowSize%2 != 1:
        windowSize = int(input('ERROR: Please enter an odd numbered window size: '))
    
    #Usage file input handling
#    print("")
#    print("You may either use one of the following preloaded codon frequency tables or input your own.")
#    for i in range(int(len(speciesList))):
#        print(i+1, speciesList[i])
    #usageChoice = input("Would you like to use one of those tables (yes or no): ")
    usageChoice="yes"

    
    usageChoice = usageChoice.lower()
    while usageChoice != "yes" and usageChoice != "no":
        usageChoice = input("Unrecognized input, please respond 'yes' or 'no': ")
        usageChoice = usageChoice.lower()
        
    #If the user wishes to use one of the included input files
    if usageChoice == "yes": 
#        print("")
    #    usageFile = input("Input species for usage file (name or number above): ")
        usageFile="2";
        usageFile = usageFile.lower()
        while usageFile not in speciesList2:
#            print("")
#            print("Unrecognized input, please enter either species name or index number: ")
#            for i in range(int(len(speciesList2)/2)):
#                print(i+1, speciesList[i])
                
            usageFile = input("Input species for usage file (name or corresponding number): ")
    
        freqDict = speciesDict[inputDict[usageFile]]
        
    #Or if they want to input one of their own
    else:
#        print("")
        print("You have the option to manually input your own codon frequencies or to import a HIVE-CUT file.")
        manFile = input("'manually' or 'file': ")
        manFile = manFile.lower()
        while manFile != 'manually' and manFile != 'file':
            print("")
            manFile = input("Unrecognized input, please respond 'manually' or 'file': ")
            manFile = manFile.lower()
            
        #If the user wants to type by hand their own file
        if manFile == 'manually':
            print("")
            print("Please enter codon frequencies in the following format: <Codon><space><frequency><space>")
            print("e.g. ATG .18 CTG .46 GTG .3 TGT .01 ...")
            print("Do not hit enter until every codon has been input")
            print("Frequencies can be in any units (e.g. percent used, frequency per 1000, observed occurences, etc.)")
            inputFreq = input("Frequencies: ")
            inputFreq = inputFreq.upper()
            inputFreq = inputFreq.split()        
            
            #Checks that input "codons" are actually codons and sequence has correct format
            codonFlag = 0
            correctCodons = 0
            while codonFlag == 0:
                correctCodons = 0
                
                while len(inputFreq)%2 !=0:
                    print("")
                    inputFreq = input("Error: Please enter frequencies in the following format: <Codon><space><frequency>")
                    inputFreq = inputFreq.upper()
                    inputFreq = inputFreq.split()
                    
                for i in range(0,len(inputFreq),2):
                    try:
                        mapDict[inputFreq[i]]
                        correctCodons += 1
                    except KeyError:
                        print("The codon " + inputFreq[i] + " was entered incorrectly.")
               
    
                if correctCodons == int(len(inputFreq)/2):
                    codonFlag = 1
                else:
                    print("")
                    print("Only enter codons as triplets comprised of ATCG")
                    inputFreq = input("Please reenter all codon frequencies correcting the above errors: ")
                    inputFreq = inputFreq.upper()
                    inputFreq = inputFreq.split()
                
            #Put input codons in kazusa like format
            newString = ""
            i = 0
            while i < len(inputFreq):
                newString += str(inputFreq[i]) + " " + str(inputFreq[i+1]) + " " + "junk" + " "
                i += 2
            frequenciesFile = []
            frequenciesFile.append(newString)
            
        #If the user wishes to include a file by typing the file path
        else:
            print("")
            print("A standard line for a HIVE-CUT uses NCBI's standard genetic code definition: ")
            print("TTT 26.18 (76390)  TCT 23.35 (68138)  TAT 19.05 (55593)   TGT  7.82 (22826)")
            print('Ensure this is the format of your file')
            print('Copying and pasting from HIVE-CUT into a text file will preserve this format')
            print("The file should be saved as a '.txt' file")
            flag = 0
            while flag == 0:
                filePath = input('Input file path: ')
                try:
                    frequenciesFile = list(open(filePath))
                    flag = 1               
                except FileNotFoundError:
                    print("No file was found at this file path, please try again")
    
        
        
    #Sequence Input Handling    
    #sequence = input('Nucleotide Sequence (start to stop codon): ')
    
    sequence = sequence.upper()
    if sequence[0:3] != "ATG":
        response = input("Sequence does not begin with ATG, would you re-enter the sequence (yes or no): ")
        response = response.lower()
        while response == "yes":
            sequence = input("Please reenter the codon sequence: ")
            if sequence[0:3] != "ATG":
                response = input("Sequence does not begin with ATG, would you like to re-enter the sequence (yes or no): ")
            else:
                response = "no"
    while len(sequence)%3 != 0:
        sequence = input("The entered nucleotide sequence is not divisible by three, please reenter the sequence:")
        sequence = sequence.upper()
    
#    print("Calculating %MinMax...")
    
    #Data Cleaning, creates dictionaries (defined in lines 5-8game length) necessary for %MinMax and harmonization
    if usageChoice == "no":
        for line in frequenciesFile:
            line = line.split()
            i=0
            if len(line)>11:  
                while i < len(line):
                    freqDict[line[i]] = float(line[i+1])
                    aaMapDict[mapDict[line[i]]] = []
                    aaFreqDict[mapDict[line[i]]] = []
                    i+=3
    
        for line in frequenciesFile:
            line = line.split()
            i = 0
            if len(line)>11:
                while i<len(line):
                    aaFreqDict[mapDict[line[i]]].append(float(line[i+1]))
                    aaMapDict[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
                    i+=3
    else: 
        for i in aaDict:#for each amino acid, initialize aaMapDict and aaFreqDict
            aaFreqDict[i] = []
            aaMapDict[i] = []
        for i in aaDict:#for each amino acid
            for j in aaDict[i]:#for each codon in that amino acid
                aaFreqDict[i].append(freqDict[j])
                aaMapDict[i].append(j + " " + str(freqDict[j]))
                         
    #For a given input fasta sequence, break into codons and corresponding amino acids
    codonSeq = []
    
    extras = ""
    for line in sequence:
        line = line.rstrip()
        string = str(extras) + str(line)
        i=0
        j=3
        while j<=len(string):
            codonSeq.append(string[i:j])
            i+=3
            j+=3
        extras = str(string[i:])
#    print (len(codonSeq))
    codonSeq=codonSeq[1:len(codonSeq)-1]
#    print (len(codonSeq))
    
    aaSeq = []
    for codon in codonSeq:
        aaSeq.append(mapDict[codon])

    
    minMaxValues = calculateMinMax(codonSeq, aaFreqDict, freqDict, mapDict, windowSize )
    
    return roundList(minMaxValues)



#==============================================================================
yeastCodonRSCUtable={"TTT"   :0.203,
"TTC"    :1.797,
"TTA"    :0.601,
"TTG"    :5.141,
"CTT"    :0.029,
"CTC"    :0.014,
"CTA"    :0.2,
"CTG"    :0.014,
"ATT"    :1.352,
"ATC"    :1.643,
"ATA"    :0.005,
"ATG"    :1,
"GTT"    :2.161,
"GTC"    :1.796,
"GTA"    :0.004,
"GTG"    :0.039,
"TAT"    :0.132,
"TAC"    :1.868,
"CAT"    :0.394,
"CAC"    :1.606,
"CAA"    :1.987,
"CAG"    :0.013,
"AAT"    :0.1,
"AAC"    :1.9,
"AAA"    :0.237,
"AAG"    :1.763,
"GAT"    :0.713,
"GAC"    :1.287,
"GAA"    :1.968,
"GAG"    :0.032,
"TCT"    :3.359,
"TCC"    :2.327,
"TCA"    :0.122,
"TCG"    :0.017,
"CCT"    :0.179,
"CCC"    :0.036,
"CCA"    :3.776,
"CCG"    :0.009,
"ACT"    :0.899,
"ACC"    :2.063,
"ACA"    :0.025,
"ACG"    :0.013,
"GCT"    :3.005,
"GCC"    :0.948,
"GCA"    :0.044,
"GCG"    :0.004,
"TGT"    :1.857,
"TGC"    :0.143,
"TGG"    :1,
"CGT"    :0.718,
"CGC"    :0.008,
"CGA"    :0.008,
"CGG"    :0.008,
"AGT"    :0.07,
"AGC"    :0.105,
"AGA"    :5.241,
"AGG"    :0.017,
"GGT"    :3.898,
"GGC"    :0.077,
"GGA"    :0.009,
"GGG"    :0.017}

def caiCal(codonList):
    print ("calculating CAI")
    global yeastCodonRSCUtable
    sequence=""
    for codon in codonList:
        sequence+=codon
    global windowSize
    isMovingWindow="yes"
    caiList=[]
    if "yes" in isMovingWindow: 
        usedCodons=windowSize
        for i in range(0,len(codonList)-windowSize+1):
            cai=1
            
            for k in range(i,i+windowSize):
                if codonList[k] in yeastCodonRSCUtable:
                    rscu=float(yeastCodonRSCUtable[codonList[k]])
                    cai*=rscu
                else:
                    print ("codon %s not in map"%codonList[k])
                    usedCodons-=1

            cai=math.pow(cai,1.0/usedCodons)
            caiList.append(cai)
    else:
        for i in range(0,len(codonList),windowSize):
            cai=1
            usedCodons=windowSize
            if i+windowSize<len(codonList):
                for k in range(i,i+windowSize):
                    if codonList[k] in yeastCodonRSCUtable:
                        rscu=yeastCodonRSCUtable[codonList[k]]
                        cai*=rscu
                    else:
                        usedCodons-=1
            cai=math.pow(cai,1/usedCodons)
            caiList.append(cai)     
    return caiList

def calCAIforGene(sequence):
    global yeastCodonRSCUtable
    codonList=loadSequence(sequence)
    usedCodons=0
    cai=1.0
    for codon in codonList:
        if codon in yeastCodonRSCUtable:
            usedCodons+=1
            rscu=yeastCodonRSCUtable[codon]
            cai*=rscu
    cai=math.pow(cai,1.0/usedCodons)
    return cai
    
    
#===============================================================================
import Bio.SeqIO
from tAI import tAI
import pandas as pd


sc_abundances = pd.Series({
  'ACG': 11070, 'AGG': 11070, 'CAG': 11070, 'CGG': 11070, 'CTC': 11070,
  'TCG': 11070, 'ATA': 22140, 'CCT': 22140, 'GAG': 22140, 'GGG': 22140,
  'GTA': 22140, 'GTG': 22140, 'AGC': 33210, 'CTA': 33210, 'GGA': 33210,
  'TCA': 33210, 'ACA': 44280, 'TGC': 44280, 'ATG': 55351, 'GCA': 55351,
  'CGT': 66421, 'TGG': 66421, 'AAA': 77491, 'CAC': 77491, 'TTA': 77491,
  'CAA': 88561, 'TAC': 88561, 'AAC': 110701, 'CCA': 110701, 'TTC': 110701,
  'TTG': 110701, 'ACT': 121771, 'AGA': 121771, 'GCT': 121771, 'TCT': 121771,
  'ATT': 143911, 'AAG': 154982, 'GAA': 154982, 'GTT': 154982, 'GAC': 177122,
  'GGC': 177122
})
tai_abundance = tAI(sc_abundances)

def caltAIForGene(sequence):
    global tai_abundance
    return tai_abundance.calc(sequence)
    









        
def setWindow(inputList,size):
    windowList=[]
    windowSize=size
    cnt=0
    
    while True:
        cnt+=1
        if cnt+windowSize>len(inputList):
            break
        selectedList=inputList[cnt:cnt+windowSize]
        sum=0
        for i in selectedList:
            sum+=float(i)   
        average=sum/len(selectedList)
        windowList.append(average)
    return windowList        

from scipy.stats.mstats import gmean


def calPhiForGene(sequence):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)
    avg=gmean(phiList)
    return avg

def calEtaForGene(sequence):
    codonList=loadSequence(sequence)
    etaList=method2(codonList)
    avg=mean(etaList)
    return avg
    

def calMForGene(sequence):
    codonList=loadSequence(sequence)
    etaList=method3(codonList)
    avg=mean(etaList)
    return avg

#method 1 to 5 is defined above method 5 is minmax,method 6 should be high phi min max, method 7 should be CAI
def results(sequence):
    sequence=parseSequence(sequence)
    codonList=loadSequence(sequence)
    resultLists=[]# store all result lists in a linear list
    method1Result=method1(codonList,"1")
    method2Result=method2(codonList)
    method3Result=method3(codonList)
    method4Result=method4(codonList)
    minMaxResult=minMax(sequence)
    
    
#    print("==")
#    print(len(method1Result))
#    print(len(method2Result))
#    print(len(method3Result))
#    print(len(method4Result))
#    print(len(minMaxResult))
#    print("==")
    resultLists.append(roundList(method1Result))
    resultLists.append(roundList(method2Result))
    resultLists.append(roundList(method3Result))
    resultLists.append(roundList(method4Result))
    resultLists.append(roundList(minMaxResult))
#    print(len(resultLists))
    return resultLists

#translate the results in window

def test():
    sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
    codonList=loadSequence(sequence)
    testResults=method4(codonList)
    print(testResults)
    print (calPhiForGene(sequence))

#
#test()
