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

def timer(operation):
    start_time = time.time()
    operation
    end_time = time.time()
    timeTaken = (start_time-end_time)
    return timeTaken

def svmPredict(trainCSV, testCSV, resultsFile):
    start_time = time.time()
    trainDF = pd.read_csv(trainCSV)
    testDF = pd.read_csv(testCSV)
    xTrain = trainDF.iloc[:, 1:-1]
    yTrain = trainDF.iloc[:, -1:]
    xTest = testDF.iloc[:, 1:-1]
    yTest = testDF.iloc[:, -1:]
    scaler1=StandardScaler()
    xTrainScaled = scaler1.fit_transform(xTrain)
    xTestScaled = scaler1.fit_transform(xTest)
    print(xTrainScaled)
    print(xTestScaled)
    mod = SVC(C=100000000, kernel='rbf', decision_function_shape='ovr')
    mod.fit(xTrainScaled, yTrain.values.ravel())
    svmResults = open(resultsFile, "w")
    predicted = mod.predict(xTestScaled)
    i=0
    print(len(predicted))
    end_time = time.time()
    for i in range(0, (len(predicted))):
        svmResults.write(">"+testDF.iloc[i, 0]+"\n")
        svmResults.write(predicted[i]+"\n")
    return(end_time-start_time)

def accuracyScore(resultsFile):
    k=0
    amountCorrect = 0
    print(k)
    for sequence in SeqIO.parse(resultsFile, "fasta"):
        for keySeq in SeqIO.parse("key.fasta", "fasta"):
            if sequence.id == keySeq.id:
                if sequence.seq==keySeq.seq:
                    amountCorrect=amountCorrect+1
        k = k+1
    print(k)
    accuracy=((amountCorrect/k)*100)
    return accuracy

rawDF = pd.read_csv('rawSvmRuns5x10percInc.csv')
print(rawDF)
rawDic = rawDF.to_dict(orient='list')
rawDic["Prediction Time Taken"] = []
rawDic["Accuracy"] = []
predictionTimes=[]
accuracyList=[]

print("fit start")
predictionTimes.append((svmPredict("train10(11000)(1)(extracted).csv", "test90(11000)(1)(extracted).csv", "shuffled11000(10-90)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train10(11000)(2)(extracted).csv", "test90(11000)(2)(extracted).csv", "shuffled11000(10-90)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train10(11000)(3)(extracted).csv", "test90(11000)(3)(extracted).csv", "shuffled11000(10-90)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train10(11000)(4)(extracted).csv", "test90(11000)(4)(extracted).csv", "shuffled11000(10-90)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train10(11000)(5)(extracted).csv", "test90(11000)(5)(extracted).csv", "shuffled11000(10-90)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train20(11000)(1)(extracted).csv", "test80(11000)(1)(extracted).csv", "shuffled11000(20-80)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train20(11000)(2)(extracted).csv", "test80(11000)(2)(extracted).csv", "shuffled11000(20-80)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train20(11000)(3)(extracted).csv", "test80(11000)(3)(extracted).csv", "shuffled11000(20-80)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train20(11000)(4)(extracted).csv", "test80(11000)(4)(extracted).csv", "shuffled11000(20-80)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train20(11000)(5)(extracted).csv", "test80(11000)(5)(extracted).csv", "shuffled11000(20-80)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train30(11000)(1)(extracted).csv", "test70(11000)(1)(extracted).csv", "shuffled11000(30-70)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train30(11000)(2)(extracted).csv", "test70(11000)(2)(extracted).csv", "shuffled11000(30-70)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train30(11000)(3)(extracted).csv", "test70(11000)(3)(extracted).csv", "shuffled11000(30-70)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train30(11000)(4)(extracted).csv", "test70(11000)(4)(extracted).csv", "shuffled11000(30-70)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train30(11000)(5)(extracted).csv", "test70(11000)(5)(extracted).csv", "shuffled11000(30-70)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train40(11000)(1)(extracted).csv", "test60(11000)(1)(extracted).csv", "shuffled11000(40-60)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train40(11000)(2)(extracted).csv", "test60(11000)(2)(extracted).csv", "shuffled11000(40-60)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train40(11000)(3)(extracted).csv", "test60(11000)(3)(extracted).csv", "shuffled11000(40-60)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train40(11000)(4)(extracted).csv", "test60(11000)(4)(extracted).csv", "shuffled11000(40-60)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train40(11000)(5)(extracted).csv", "test60(11000)(5)(extracted).csv", "shuffled11000(40-60)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train50(11000)(1)(extracted).csv", "test50(11000)(1)(extracted).csv", "shuffled11000(50-50)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train50(11000)(2)(extracted).csv", "test50(11000)(2)(extracted).csv", "shuffled11000(50-50)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train50(11000)(3)(extracted).csv", "test50(11000)(3)(extracted).csv", "shuffled11000(50-50)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train50(11000)(4)(extracted).csv", "test50(11000)(4)(extracted).csv", "shuffled11000(50-50)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train50(11000)(5)(extracted).csv", "test50(11000)(5)(extracted).csv", "shuffled11000(50-50)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train60(11000)(1)(extracted).csv", "test40(11000)(1)(extracted).csv", "shuffled11000(60-40)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train60(11000)(2)(extracted).csv", "test40(11000)(2)(extracted).csv", "shuffled11000(60-40)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train60(11000)(3)(extracted).csv", "test40(11000)(3)(extracted).csv", "shuffled11000(60-40)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train60(11000)(4)(extracted).csv", "test40(11000)(4)(extracted).csv", "shuffled11000(60-40)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train60(11000)(5)(extracted).csv", "test40(11000)(5)(extracted).csv", "shuffled11000(60-40)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train70(11000)(1)(extracted).csv", "test30(11000)(1)(extracted).csv", "shuffled11000(70-30)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train70(11000)(2)(extracted).csv", "test30(11000)(2)(extracted).csv", "shuffled11000(70-30)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train70(11000)(3)(extracted).csv", "test30(11000)(3)(extracted).csv", "shuffled11000(70-30)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train70(11000)(4)(extracted).csv", "test30(11000)(4)(extracted).csv", "shuffled11000(70-30)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train70(11000)(5)(extracted).csv", "test30(11000)(5)(extracted).csv", "shuffled11000(70-30)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train80(11000)(1)(extracted).csv", "test20(11000)(1)(extracted).csv", "shuffled11000(80-20)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train80(11000)(2)(extracted).csv", "test20(11000)(2)(extracted).csv", "shuffled11000(80-20)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train80(11000)(3)(extracted).csv", "test20(11000)(3)(extracted).csv", "shuffled11000(80-20)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train80(11000)(4)(extracted).csv", "test20(11000)(4)(extracted).csv", "shuffled11000(80-20)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train80(11000)(5)(extracted).csv", "test20(11000)(5)(extracted).csv", "shuffled11000(80-20)(5)(SVMresults).fasta")))

predictionTimes.append((svmPredict("train90(11000)(1)(extracted).csv", "test10(11000)(1)(extracted).csv", "shuffled11000(90-10)(1)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train90(11000)(2)(extracted).csv", "test10(11000)(2)(extracted).csv", "shuffled11000(90-10)(2)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train90(11000)(3)(extracted).csv", "test10(11000)(3)(extracted).csv", "shuffled11000(90-10)(3)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train90(11000)(4)(extracted).csv", "test10(11000)(4)(extracted).csv", "shuffled11000(90-10)(4)(SVMresults).fasta")))
predictionTimes.append((svmPredict("train90(11000)(5)(extracted).csv", "test10(11000)(5)(extracted).csv", "shuffled11000(90-10)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(10-90)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(10-90)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(10-90)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(10-90)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(10-90)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(20-80)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(20-80)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(20-80)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(20-80)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(20-80)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(30-70)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(30-70)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(30-70)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(30-70)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(30-70)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(40-60)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(40-60)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(40-60)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(40-60)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(40-60)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(50-50)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(50-50)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(50-50)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(50-50)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(50-50)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(60-40)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(60-40)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(60-40)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(60-40)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(60-40)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(70-30)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(70-30)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(70-30)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(70-30)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(70-30)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(80-20)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(80-20)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(80-20)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(80-20)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(80-20)(5)(SVMresults).fasta")))

accuracyList.append((accuracyScore("shuffled11000(90-10)(1)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(90-10)(2)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(90-10)(3)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(90-10)(4)(SVMresults).fasta")))
accuracyList.append((accuracyScore("shuffled11000(90-10)(5)(SVMresults).fasta")))

rawDic['Prediction Time Taken']=predictionTimes
rawDic['Accuracy']=accuracyList
print(rawDic)

rawDF = pd.DataFrame.from_dict(rawDic)
rawDF.to_csv("rawSvmRuns5x10percInc.csv", index=False)

print("fit end")