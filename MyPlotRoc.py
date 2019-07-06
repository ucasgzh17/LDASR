# 传入的为y_prediction, y_real（numpy.array），输出的为x（list）, y（list）, auc

import copy
def plot_ROC(y_prediction, y_real):

    # 把numpy.array转换为list
    prediction1 = y_prediction.tolist()
    y_test1 = y_real.tolist()

    # 生成prediction, y_test的副本
    prediction = prediction1.copy()
    y_test = y_test1.copy()

    # 计算正负label
    pos_num = 0
    neg_num = 0
    counter = 0
    while counter < len(y_test):
        if y_test[counter] == 1:
            pos_num = pos_num + 1
        if y_test[counter] == 0:
            neg_num = neg_num + 1
        counter = counter + 1

    # 生成排序，根据prediction对y_test排序
    PT = []
    counter = 0
    while counter < len(y_test):
        pair = []
        pair.append(prediction[counter])
        pair.append(y_test[counter])
        PT.append(pair)
        counter = counter + 1
    # print('sort before')      # ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    # print(PT)
    counter = 0
    a1 = []
    b1 = []
    while counter < len(PT):
        if PT[counter][0] == 0:
            a1.append(PT[counter])
        if PT[counter][0] == 1:
            b1.append(PT[counter])
        counter = counter + 1
    a1.extend(b1)
    # print('sort after')       # ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    # print(a1)

    # 计算各个点坐标
    auc = 0
    x = [1]
    y = [1]
    counter = 1
    while counter < len(y_test):
        FP = 0
        TP = 0
        counter1 = counter
        while counter1 < len(y_test):
            if a1[counter1][1] == 1:
                TP = TP + 1
            if a1[counter1][1] == 0:
                FP = FP + 1
            counter1 = counter1 + 1
        x.append(FP / neg_num)
        y.append(TP / pos_num)
        auc = auc + (y[counter] + y[counter - 1]) * (x[counter - 1] - x[counter]) / 2
        counter = counter + 1
    x.append(0)
    y.append(0)
    # print(x)  # ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    # print(y)
    m = len(y_test)
    auc = auc + y[m] * x[m] / 2
    return x, y, auc
