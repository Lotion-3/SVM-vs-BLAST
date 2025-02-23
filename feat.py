from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
import time
import timeit


#Feature Extraction
def extractedDictionary(dicName):
    dicName = {'id':[], 'G perc':[], 'C perc':[], 'A perc':[], 'T perc':[], 'GC perc':[], 'Alanine perc': [], 'Arginine perc': [], 'Asparagine perc': [], 'Aspartic Acid perc': [], 'Cysteine perc': [], 'Glutamine perc': [], 'Glutamic Acid perc': [], 'Glycine perc': [], 'Histidine perc': [], 'Isoleucine perc': [], 'Leucine perc': [], 'Lysine perc': [], 'Methionine perc': [], 'Phenylalanine perc': [], 'Proline perc': [], 'Serine perc': [], 'Threonine perc': [], 'Tryptophan perc': [], 'Tyrosine perc': [], 'Valine perc': [], 'Variant':[]}
    return dicName
def extractedDF(splitFile, extractedCsvName):
    startTime=time.time()
    dicName = {'id':[], 'G perc':[], 'C perc':[], 'A perc':[], 'T perc':[], 'GC perc':[], 'Alanine perc': [], 'Arginine perc': [], 'Asparagine perc': [], 'Aspartic Acid perc': [], 'Cysteine perc': [], 'Glutamine perc': [], 'Glutamic Acid perc': [], 'Glycine perc': [], 'Histidine perc': [], 'Isoleucine perc': [], 'Leucine perc': [], 'Lysine perc': [], 'Methionine perc': [], 'Phenylalanine perc': [], 'Proline perc': [], 'Serine perc': [], 'Threonine perc': [], 'Tryptophan perc': [], 'Tyrosine perc': [], 'Valine perc': [], 'Variant':[]}
    idList=[]
    gList=[]
    cList=[]
    aList=[]
    tList=[]
    gcList=[]

    alaList=[]
    argList=[]
    asnList=[]
    aspList=[]
    cysList=[]
    glnList=[]
    gluList=[]
    glyList=[]
    hisList=[]
    ileList=[]
    leuList=[]
    lysList=[]
    metList=[]
    pheList=[]
    proList=[]
    serList=[]
    thrList=[]
    trpList=[]
    tyrList=[]
    valList=[]

    variantList=[]

    i=0
    for sequence in SeqIO.parse(splitFile, "fasta"):
        #print(sequence.id)
        seqLength=(len(sequence))
        seqG=(sequence.count("G"))
        seqC=(sequence.count("C"))
        seqA=(sequence.count("A"))
        seqT=(sequence.count("T"))
        percG=(seqG/seqLength)*1
        percC=(seqC/seqLength)*1
        percA=(seqA/seqLength)*1
        percT=(seqT/seqLength)*1
        percGC=(gc_fraction(sequence))*1
        aminoSeq= ProteinAnalysis(Seq.translate(sequence.seq))

        percAla = aminoSeq.get_amino_acids_percent()['A']*1
        percArg = aminoSeq.get_amino_acids_percent()['R']*1
        percAsn = aminoSeq.get_amino_acids_percent()['N']*1
        percAsp = aminoSeq.get_amino_acids_percent()['D']*1
        percCys = aminoSeq.get_amino_acids_percent()['C']*1
        percGln = aminoSeq.get_amino_acids_percent()['Q']*1
        percGlu = aminoSeq.get_amino_acids_percent()['E']*1
        percGly = aminoSeq.get_amino_acids_percent()['G']*1
        percHis = aminoSeq.get_amino_acids_percent()['H']*1
        percIle = aminoSeq.get_amino_acids_percent()['I']*1
        percLeu = aminoSeq.get_amino_acids_percent()['L']*1
        percLys = aminoSeq.get_amino_acids_percent()['K']*1
        percMet = aminoSeq.get_amino_acids_percent()['M']*1
        percPhe = aminoSeq.get_amino_acids_percent()['F']*1
        percPro = aminoSeq.get_amino_acids_percent()['P']*1
        percSer = aminoSeq.get_amino_acids_percent()['S']*1
        percThr = aminoSeq.get_amino_acids_percent()['T']*1
        percTrp = aminoSeq.get_amino_acids_percent()['W']*1
        percTyr = aminoSeq.get_amino_acids_percent()['Y']*1
        percVal = aminoSeq.get_amino_acids_percent()['V']*1

        #print(percAla)
        #print('Perc G')
        #print(percG)
        #print('Perc C')
        #print(percC)
        #print('Perc A')
        #print(percA)
        #print('Perc T')
        #print(percT)
        #print('Perc GC')
        #print(percGC)
        #print(gc_fraction(sequence))
        idList.append(sequence.id)
        gList.append(percG)
        cList.append(percC)
        aList.append(percA)
        tList.append(percT)
        gcList.append(percGC)

        alaList.append(percAla)
        argList.append(percArg)
        asnList.append(percAsn)
        aspList.append(percAsp)
        cysList.append(percCys)
        glnList.append(percGln)
        gluList.append(percGlu)
        glyList.append(percGly)
        hisList.append(percHis)
        ileList.append(percIle)
        leuList.append(percLeu)
        lysList.append(percLys)
        metList.append(percMet)
        pheList.append(percPhe)
        proList.append(percPro)
        serList.append(percSer)
        thrList.append(percThr)
        trpList.append(percTrp)
        tyrList.append(percTyr)
        valList.append(percVal)

        for subSpec in SeqIO.parse("key.fasta", "fasta"):
            if sequence.id == subSpec.id:
                variantList.append(str(subSpec.seq))
        i=i+1
    
    print(i)
    print(len(gList))
    print(len(cList))
    print(len(aList))
    print(len(tList))
    print(len(gcList))

    print(len(alaList))
    print(len(argList))
    print(len(asnList))
    print(len(aspList))
    print(len(cysList))
    print(len(glnList))
    print(len(gluList))
    print(len(glyList))
    print(len(hisList))
    print(len(ileList))
    print(len(leuList))
    print(len(lysList))
    print(len(metList))
    print(len(pheList))
    print(len(proList))
    print(len(serList))
    print(len(thrList))
    print(len(trpList))
    print(len(tyrList))
    print(len(valList))

    print(len(variantList))
    
    dicName['id']=idList
    dicName['G perc']=gList
    dicName['C perc']=cList
    dicName['A perc']=aList
    dicName['T perc']=tList
    dicName['GC perc']=gcList

    dicName['Alanine perc']=alaList
    dicName['Arginine perc']=argList
    dicName['Asparagine perc']=asnList
    dicName['Aspartic Acid perc']=aspList
    dicName['Cysteine perc']=cysList
    dicName['Glutamine perc']=glnList
    dicName['Glutamic Acid perc']=gluList
    dicName['Glycine perc']=glyList
    dicName['Histidine perc']=hisList
    dicName['Isoleucine perc']=ileList
    dicName['Leucine perc']=leuList
    dicName['Lysine perc']=lysList
    dicName['Methionine perc']=metList
    dicName['Phenylalanine perc']=pheList
    dicName['Proline perc']=proList
    dicName['Serine perc']=serList
    dicName['Threonine perc']=thrList
    dicName['Tryptophan perc']=trpList
    dicName['Tyrosine perc']=tyrList
    dicName['Valine perc']=valList

    dicName['Variant']=variantList
    dfName = pd.DataFrame.from_dict(dicName)
    dfName.to_csv(extractedCsvName, index=False)

    endTime=time.time()
    print(endTime-startTime)
    return(endTime-startTime)

def timer(operation):
    start_time = time.time()
    operation
    end_time = time.time()
    timeTaken = (end_time-start_time)
    print(start_time)
    print(end_time)
    print(timeTaken)
    return timeTaken

def counting(num):
    print(num)

rawDF = pd.read_csv('rawSvmRuns5x10percInc.csv')
print(rawDF)
rawDic = rawDF.to_dict(orient='list')
rawDic["Train Extract Time Taken"] = []
rawDic["Test Extract Time Taken"] = []
trainExtractTimes=[]
testExtractTimes=[]

trainExtractTimes.append((extractedDF("train10(11000)(1).fasta", "train10(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train10(11000)(2).fasta", "train10(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train10(11000)(3).fasta", "train10(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train10(11000)(4).fasta", "train10(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train10(11000)(5).fasta", "train10(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train20(11000)(1).fasta", "train20(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train20(11000)(2).fasta", "train20(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train20(11000)(3).fasta", "train20(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train20(11000)(4).fasta", "train20(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train20(11000)(5).fasta", "train20(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train30(11000)(1).fasta", "train30(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train30(11000)(2).fasta", "train30(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train30(11000)(3).fasta", "train30(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train30(11000)(4).fasta", "train30(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train30(11000)(5).fasta", "train30(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train40(11000)(1).fasta", "train40(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train40(11000)(2).fasta", "train40(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train40(11000)(3).fasta", "train40(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train40(11000)(4).fasta", "train40(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train40(11000)(5).fasta", "train40(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train50(11000)(1).fasta", "train50(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train50(11000)(2).fasta", "train50(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train50(11000)(3).fasta", "train50(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train50(11000)(4).fasta", "train50(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train50(11000)(5).fasta", "train50(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train60(11000)(1).fasta", "train60(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train60(11000)(2).fasta", "train60(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train60(11000)(3).fasta", "train60(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train60(11000)(4).fasta", "train60(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train60(11000)(5).fasta", "train60(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train70(11000)(1).fasta", "train70(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train70(11000)(2).fasta", "train70(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train70(11000)(3).fasta", "train70(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train70(11000)(4).fasta", "train70(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train70(11000)(5).fasta", "train70(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train80(11000)(1).fasta", "train80(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train80(11000)(2).fasta", "train80(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train80(11000)(3).fasta", "train80(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train80(11000)(4).fasta", "train80(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train80(11000)(5).fasta", "train80(11000)(5)(extracted).csv")))

trainExtractTimes.append((extractedDF("train90(11000)(1).fasta", "train90(11000)(1)(extracted).csv")))
trainExtractTimes.append((extractedDF("train90(11000)(2).fasta", "train90(11000)(2)(extracted).csv")))
trainExtractTimes.append((extractedDF("train90(11000)(3).fasta", "train90(11000)(3)(extracted).csv")))
trainExtractTimes.append((extractedDF("train90(11000)(4).fasta", "train90(11000)(4)(extracted).csv")))
trainExtractTimes.append((extractedDF("train90(11000)(5).fasta", "train90(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test90(11000)(1).fasta", "test90(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test90(11000)(2).fasta", "test90(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test90(11000)(3).fasta", "test90(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test90(11000)(4).fasta", "test90(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test90(11000)(5).fasta", "test90(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test80(11000)(1).fasta", "test80(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test80(11000)(2).fasta", "test80(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test80(11000)(3).fasta", "test80(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test80(11000)(4).fasta", "test80(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test80(11000)(5).fasta", "test80(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test70(11000)(1).fasta", "test70(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test70(11000)(2).fasta", "test70(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test70(11000)(3).fasta", "test70(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test70(11000)(4).fasta", "test70(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test70(11000)(5).fasta", "test70(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test60(11000)(1).fasta", "test60(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test60(11000)(2).fasta", "test60(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test60(11000)(3).fasta", "test60(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test60(11000)(4).fasta", "test60(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test60(11000)(5).fasta", "test60(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test50(11000)(1).fasta", "test50(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test50(11000)(2).fasta", "test50(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test50(11000)(3).fasta", "test50(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test50(11000)(4).fasta", "test50(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test50(11000)(5).fasta", "test50(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test40(11000)(1).fasta", "test40(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test40(11000)(2).fasta", "test40(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test40(11000)(3).fasta", "test40(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test40(11000)(4).fasta", "test40(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test40(11000)(5).fasta", "test40(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test30(11000)(1).fasta", "test30(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test30(11000)(2).fasta", "test30(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test30(11000)(3).fasta", "test30(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test30(11000)(4).fasta", "test30(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test30(11000)(5).fasta", "test30(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test20(11000)(1).fasta", "test20(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test20(11000)(2).fasta", "test20(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test20(11000)(3).fasta", "test20(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test20(11000)(4).fasta", "test20(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test20(11000)(5).fasta", "test20(11000)(5)(extracted).csv")))

testExtractTimes.append((extractedDF("test10(11000)(1).fasta", "test10(11000)(1)(extracted).csv")))
testExtractTimes.append((extractedDF("test10(11000)(2).fasta", "test10(11000)(2)(extracted).csv")))
testExtractTimes.append((extractedDF("test10(11000)(3).fasta", "test10(11000)(3)(extracted).csv")))
testExtractTimes.append((extractedDF("test10(11000)(4).fasta", "test10(11000)(4)(extracted).csv")))
testExtractTimes.append((extractedDF("test10(11000)(5).fasta", "test10(11000)(5)(extracted).csv")))

rawDic['Train Extract Time Taken']=trainExtractTimes
rawDic['Test Extract Time Taken']=testExtractTimes
print(rawDic)

rawDF = pd.DataFrame.from_dict(rawDic)
rawDF.to_csv("rawSvmRuns5x10percInc.csv", index=False)

print("extract finished")