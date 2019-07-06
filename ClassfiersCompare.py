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
SampleFeature = []
SampleLabel = []
SampleFeature = RSampleFeature
SampleLabel = RSampleLabel


X = np.array(SampleFeature)
y = np.array(SampleLabel)


# 训练
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from MyPlotRoc import plot_ROC
from sklearn.metrics import roc_curve, auc
from scipy import interp
from itertools import cycle

# RotationForest    ok!
from RotationForest import RotationForest
model1 = RotationForest()

# RandomForestClassifier ok!
from sklearn.ensemble import RandomForestClassifier
model2 = RandomForestClassifier()

# LogisticRegression  ok!
from sklearn.linear_model import LogisticRegression
model3 = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')

from sklearn.naive_bayes import GaussianNB
model4 = GaussianNB()

from sklearn.svm import SVC
model5 = SVC(probability=True)




# k折交叉验证
SplitNum = 5
cv = StratifiedKFold(n_splits=SplitNum)


# RotationForest
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
i = 0

for train, test in cv.split(X, y):
    print('x')
    model1 = RotationForest()
    model1.fit(X[train], y[train])
    y_score1 = model1.predict_proba(X[test])

    print(i)

    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)

    i += 1

# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='RotationForest (AUC = %0.4f)' % (roc_auc), linestyle='-',color='black',
         lw=2, alpha=.8)

# Random Forest
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
i = 0

for train, test in cv.split(X, y):
    print('x')

    model2.fit(X[train], y[train])
    y_score1 = model2.predict_proba(X[test])

    print(i)

    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    i += 1

# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='RandomForest (AUC = %0.4f)' % (roc_auc), linestyle='--',
         lw=2, alpha=.8)


# LogisticRegression
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
i = 0

for train, test in cv.split(X, y):
    print('x')

    model3.fit(X[train], y[train])
    y_score1 = model3.predict_proba(X[test])

    print(i)

    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    # plt.plot(fpr, tpr, lw=1, alpha=0.3,
    #          label='ROC %d (AUC = %0.2f)' % (i, roc_auc))

    i += 1

# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='LogisticRegression (AUC = %0.4f)' % (roc_auc), linestyle='-.',
         lw=2, alpha=.8)


# NaiveBayes
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
i = 0

for train, test in cv.split(X, y):
    print('x')
    model4.fit(X[train], y[train])
    y_score1 = model4.predict_proba(X[test])

    print(i)

    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)

    i += 1

# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='NaiveBayes (AUC = %0.4f)' % (roc_auc), linestyle=':',
         lw=2, alpha=.8)

# SVM
from sklearn.svm import SVC
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
i = 0

for train, test in cv.split(X, y):
    print('x')
    model5 = SVC(degree=5, coef0=0.5, C=1, probability=True)

    model5.fit(X[train], y[train])
    y_score1 = model5.predict_proba(X[test])

    print(i)

    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)

    i += 1

# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='SVM (AUC = %0.4f)' % (roc_auc), linestyle=':', marker='|',
         lw=2, alpha=.8)



# 画标题坐标轴
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate',fontsize=13)
plt.ylabel('True Positive Rate',fontsize=13)
# plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('ClassifierCompare1.tif',dpi=300)
plt.show()