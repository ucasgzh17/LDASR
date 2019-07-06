import os
import time
import numpy as np
import pandas as pd
import csv
import math
import random
from sklearn.model_selection import KFold,StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from RotationForest import RotationForest
from sklearn.ensemble import RandomForestClassifier

from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from MyPlotRoc import plot_ROC
from sklearn.metrics import roc_curve, auc
from scipy import interp
from itertools import cycle
from sklearn.model_selection import KFold,LeaveOneOut,LeavePOut,ShuffleSplit
from RotationForest import RotationForest

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

LOOCVRealAndPrediction = []
ReadMyCsv(LOOCVRealAndPrediction, 'LOOCVRealAndPrediction.csv')


y = []
prediction = []
counter = 0
while counter < len(LOOCVRealAndPrediction):
    y.append(LOOCVRealAndPrediction[counter][0])
    prediction.append(LOOCVRealAndPrediction[counter][1])
    counter = counter + 1


# 画图
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fpr, tpr, thresholds = roc_curve(y, prediction)
print(fpr)
print(len(fpr))
print(tpr)
print(len(tpr))
# 增加零点
fpr = fpr.tolist()
tpr = tpr.tolist()
fpr.insert(0, 0)
tpr.insert(0, 0)
roc_auc = auc(fpr, tpr)
aucs.append(roc_auc)

# from IPython.core.pylabtools import figsize
# figsize(9, 6)

plt.plot(fpr, tpr, lw=2, color='black', alpha=0.8,
         label='ROC LOOCV (AUC = %0.4f)' % (roc_auc))
print(fpr)
print(len(fpr))
print(tpr)
print(len(tpr))
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate',fontsize=13)
plt.ylabel('True Positive Rate',fontsize=13)

# plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
# plt.savefig('LOOCV.svg')
plt.savefig('LOOCV.tif',dpi=300)
plt.show()