import matplotlib.pyplot as plt
import numpy as np
import csv
import random


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


def MyEnlarge(x0, y0, width, height, x1, y1, times, mean_fpr, mean_tpr, thickness=1, color = 'blue'):
    # 第一个框的坐标，大小，第二个框的坐标，倍数，全部图像的fpr，tpr传入，粗细，颜色
    def MyFrame(x0, y0, width, height):
        import matplotlib.pyplot as plt
        import numpy as np

        x1 = np.linspace(x0, x0, num=20)  # 生成列的横坐标，横坐标都是x0，纵坐标变化
        y1 = np.linspace(y0, y0, num=20)
        xk = np.linspace(x0, x0 + width, num=20)
        yk = np.linspace(y0, y0 + height, num=20)

        xkn = []
        ykn = []
        counter = 0
        while counter < 20:
            xkn.append(x1[counter] + width)
            ykn.append(y1[counter] + height)
            counter = counter + 1

        plt.plot(x1, yk, color='k', linestyle=':', lw=1, alpha=1)  # 左
        plt.plot(xk, y1, color='k', linestyle=':', lw=1, alpha=1)  # 下
        plt.plot(xkn, yk, color='k', linestyle=':', lw=1, alpha=1)  # 右
        plt.plot(xk, ykn, color='k', linestyle=':', lw=1, alpha=1)  # 上

        return
    # 画虚线框
    width2 = times * width
    height2 = times * height
    MyFrame(x0, y0, width, height)
    MyFrame(x1, y1, width2, height2)

    # 连接两个虚线框
    xp = np.linspace(x0 + width, x1, num=20)
    yp = np.linspace(y0, y1 + height2, num=20)
    plt.plot(xp, yp, color='k', linestyle=':', lw=1, alpha=1)

    # 小虚框内各点坐标
    XDottedLine = []
    YDottedLine = []
    counter = 0
    while counter < len(mean_fpr):
        if mean_fpr[counter] > x0 and mean_fpr[counter] < (x0 + width) and mean_tpr[counter] > y0 and mean_tpr[counter] < (y0 + height):
            XDottedLine.append(mean_fpr[counter])
            YDottedLine.append(mean_tpr[counter])
        counter = counter + 1

    # 画虚线框内的点
    # 把小虚框内的任一点减去小虚框左下角点生成相对坐标，再乘以倍数（4）加大虚框左下角点
    counter = 0
    while counter < len(XDottedLine):
        XDottedLine[counter] = (XDottedLine[counter] - x0) * times + x1
        YDottedLine[counter] = (YDottedLine[counter] - y0) * times + y1
        counter = counter + 1


    plt.plot(XDottedLine, YDottedLine, linestyle='--', color=color, lw=thickness, alpha=1)
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

SplitNum = 5
cv = StratifiedKFold(n_splits=SplitNum)

model2 = RandomForestClassifier()
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 1000)  # ！！！！！！！！！！！！！！！！！插值


colorlist = ['red', 'gold', 'purple', 'green', 'blue', 'black']

i = 0
for train, test in cv.split(X, y):              # https://yq.aliyun.com/ziliao/286979
    print('!')
    model2.fit(X[train], y[train])
    y_score1 = model2.predict_proba(X[test])
    fpr, tpr, thresholds = roc_curve(y[test], y_score1[:, 1])
    p = mean_fpr                                # 横坐标为插值！！！！！！
    q = interp(mean_fpr, fpr, tpr)              # 纵坐标为由很坐标产生的插值！！！！！！
    tprs.append(interp(mean_fpr, fpr, tpr))         # ctrl + b 看源码！！！！！！
    # print(fpr)

    tprs.append(interp(mean_fpr, fpr, tpr))         # 插横坐标，得纵坐标。横坐标mean_fpr，纵坐标tprs[counter]
    tprs[-1][0] = 0.0
    MyEnlarge(0.075, 0.75, 0.1, 0.1, 0.15, 0.25, 4, p, q, 1, colorlist[i])  # MyEnlarge(x0, y0, width, height, x1, y1, times, mean_fpr, mean_tpr)
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    plt.plot(fpr, tpr, lw=1, alpha=0.3, color=colorlist[i],
             label='Fold %d (AUC = %0.4f)' % (i, roc_auc))

    i += 1


# 画均值
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
roc_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, label='Mean (AUC = %0.4f)' % (roc_auc), linestyle='--', lw=2, alpha=.8, color=colorlist[5])


# 画虚线框
MyEnlarge(0.075, 0.75, 0.1, 0.1, 0.15, 0.25, 4, mean_fpr, mean_tpr, 2, colorlist[5])


# 画坐标轴
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])


# 画标题
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()