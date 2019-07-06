import matplotlib.pyplot as plt
import numpy as np
import csv
import random
# from pandas_ml import ConfusionMatrix         # 陈博士的confusion matrix！！！
# cm = ConfusionMatrix(y_test, prediction)
# print(cm)
# cm.print_stats()


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

def MyConfusionMatrix(y_real,y_predict):
    from sklearn.metrics import confusion_matrix
    CM = confusion_matrix(y_real, y_predict)
    print(CM)
    CM = CM.tolist()
    TN = CM[0][0]
    FP = CM[0][1]
    FN = CM[1][0]
    TP = CM[1][1]
    print('TN:%d, FP:%d, FN:%d, TP:%d' % (TN, FP, FN, TP))
    Acc = (TN + TP) / (TN + TP + FN + FP)
    Sen = TP / (TP + FN)
    Spec = TN / (TN + FP)
    Prec = TP / (TP + FP)
    MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    # 分母可能出现0，需要讨论待续
    print('Acc:', Acc)
    print('Sen:', Sen)
    print('Spec:', Spec)
    print('Prec:', Prec)
    print('Mcc:', MCC)
    return


# 读取源文件
SampleFeature = []
ReadMyCsv(SampleFeature, "encoded_imgsMany.csv")
print(len(SampleFeature))
print(len(SampleFeature[0]))

# 转换数据类型
counter = 0
while counter < len(SampleFeature):
    counter1 = 0
    while counter1 < len(SampleFeature[counter]):
        SampleFeature[counter][counter1] = float(SampleFeature[counter][counter1])
        counter1 = counter1 + 1
    counter = counter + 1

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
# print(SampleLabel)
# storFile(SampleLabel, 'SampleLabel.csv')      # 保存csv会报错，改变数据类型


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


# 训练
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from scipy import interp
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import train_test_split
x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

model = RandomForestClassifier()
model.fit(x_train, y_train)
prediction = model.predict(x_test)

import math
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

print(classification_report(y_test, prediction))
MyConfusionMatrix(y_test, prediction)



