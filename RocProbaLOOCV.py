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

# 定义函数
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return

def storFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return


# 读取源文件
SampleFeature = []
ReadMyCsv(SampleFeature, "SampleFeatureAuto.csv")


# SampleLabel
SampleLabel = []
counter = 0
while counter < len(SampleFeature) / 2:
    # Row = []
    # Row.append(1)
    SampleLabel.append(1)
    counter = counter + 1
counter1 = 0
while counter1 < len(SampleFeature) / 2:
    # Row = []
    # Row.append(0)
    SampleLabel.append(0)
    counter1 = counter1 + 1

# 打乱数据集顺序
counter = 0
R = []
while counter < len(SampleFeature):
    R.append(counter)
    counter = counter + 1
random.shuffle(R)

RSampleFeature = []
RSampleLabel = []
counter = 0
while counter < len(SampleFeature):
    RSampleFeature.append(SampleFeature[R[counter]])
    RSampleLabel.append(SampleLabel[R[counter]])
    counter = counter + 1
print('len(RSampleFeature)', len(RSampleFeature))
print('len(RSampleLabel)', len(RSampleLabel))

SampleFeature = []
SampleLabel = []
SampleFeature = RSampleFeature
SampleLabel = RSampleLabel

X = np.array(SampleFeature)
y = np.array(SampleLabel)

# X = SampleFeature
# y = SampleLabel

# 训练
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from MyPlotRoc import plot_ROC
from sklearn.metrics import roc_curve, auc
from scipy import interp
from itertools import cycle
from sklearn.model_selection import KFold,LeaveOneOut,LeavePOut,ShuffleSplit
from RotationForest import RotationForest

loo = LeaveOneOut()

# SplitNum = 5
# cv = StratifiedKFold(n_splits=SplitNum)

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

counter = 0
prediction = []
for train, test in loo.split(X, y):
    print('x')
    model1 = RotationForest()
    model1.fit(X[train], y[train])
    predict = model1.predict_proba(X[test])
    print('!')
    prediction.append(predict[0][1])
    counter = counter + 1
    print(counter)

print(prediction)
# 跑的慢，存一下LOOCV-Prediction
RealAndPrediction = []
counter = 0
while counter < len(prediction):
    pair = []
    pair.append(y[counter])
    pair.append(prediction[counter])
    RealAndPrediction.append(pair)
    counter = counter + 1

storFile(RealAndPrediction, 'RealAndPrediction.csv')

prediction = np.array(prediction)


fpr, tpr, thresholds = roc_curve(y, prediction)


roc_auc = auc(fpr, tpr)
aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, color='b', alpha=0.8,
         label='ROC LOOCV (AUC = %0.2f)' % (roc_auc))


plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate',fontsize=13)
plt.ylabel('True Positive Rate',fontsize=13)
# plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig('LOOCV.png')
plt.show()