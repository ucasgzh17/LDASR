import csv
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import svds
from pylab import *
from sklearn.decomposition import NMF
import Tool

# 定义函数
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        for i in range(len(row)):
            row[i] = float(row[i])
        SaveList.append(row)
    return

def storFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return

DiseaseAndRNABinary = []
ReadMyCsv(DiseaseAndRNABinary, 'DiseaseAndRNABinary.csv')

RNAGaussian = []
ReadMyCsv(RNAGaussian, 'RNAGaussian.csv')
DiseaseGaussian = []
ReadMyCsv(DiseaseGaussian, 'DiseaseGaussian.csv')

DiseaseSimilarityModel1 = []
ReadMyCsv(DiseaseSimilarityModel1, 'DiseaseSimilarityModel1.csv')
DiseaseSimilarityModel2 = []
ReadMyCsv(DiseaseSimilarityModel2, 'DiseaseSimilarityModel2.csv')


DS = []
counter = 0
while counter < len(DiseaseSimilarityModel1):
    row = []
    counter1 = 0
    while counter1 < len(DiseaseSimilarityModel1[counter]):
        if (DiseaseSimilarityModel1[counter][counter1] == 0) & (DiseaseSimilarityModel2[counter][counter1] == 0):
            row.append(DiseaseGaussian[counter][counter1])
        else:
            v = DiseaseSimilarityModel1[counter][counter1] + DiseaseSimilarityModel2[counter][counter1]
            row.append(v/2)
        counter1 = counter1 + 1
    DS.append(row)
    counter = counter + 1
    print(counter)

storFile(DS, 'DS.csv')


PositiveSampleFeature = Tool.PositiveGenerate(DiseaseAndRNABinary, RNAGaussian, DS)
storFile(PositiveSampleFeature, 'PositiveSampleFeature.csv')

NegativeSampleFeature, NegativeSample = Tool.NegativeGenerate(RNAGaussian, DS)
storFile(NegativeSampleFeature, 'NegativeSampleFeature.csv')
storFile(NegativeSample, 'NegativeSample.csv')

