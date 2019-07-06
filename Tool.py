import numpy as np
from sklearn.decomposition import NMF
from pylab import *
import matplotlib.pyplot as plt
import csv
import random
from scipy.sparse.linalg import svds
from scipy import sparse
from numpy import *
import numpy as np

# 定义函数
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        for i in range(len(row)):       # 转换数据类型
            row[i] = float(row[i])
        SaveList.append(row)
    return

def ReadMyCsv2(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:          # 注意表头
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return

# 构建AllDiseaseOld
def GenerateAllDisease(LncRNADiseaseAssociation):
    AllDisease = []
    counter1 = 0
    while counter1 < len(LncRNADiseaseAssociation):  # 顺序遍历原始数据，构建AllDisease
        counter2 = 0
        flag = 0
        while counter2 < len(AllDisease):  # 遍历AllDisease
            if LncRNADiseaseAssociation[counter1][1] != AllDisease[counter2]:  # 有新疾病
                counter2 = counter2 + 1
            elif LncRNADiseaseAssociation[counter1][1] == AllDisease[counter2]:  # 没有新疾病，用两个if第二个if会越界
                flag = 1
                counter2 = counter2 + 1
        if flag == 0:
            AllDisease.append(LncRNADiseaseAssociation[counter1][1])
        counter1 = counter1 + 1
    return AllDisease

def GenerateAllRNA(LncRNADiseaseAssociation):
    # 构建AllRNA
    AllRNA = []
    counter1 = 0
    while counter1 < len(LncRNADiseaseAssociation):  # 顺序遍历原始数据，构建AllDisease
        counter2 = 0
        flag = 0
        while counter2 < len(AllRNA):  # 遍历AllDisease
            if LncRNADiseaseAssociation[counter1][0] != AllRNA[counter2]:  # 有新疾病
                counter2 = counter2 + 1
            elif LncRNADiseaseAssociation[counter1][0] == AllRNA[counter2]:  # 没有新疾病，用两个if第二个if会越界
                flag = 1
                break
        if flag == 0:
            AllRNA.append(LncRNADiseaseAssociation[counter1][0])
        counter1 = counter1 + 1
    return AllRNA

# 小写
def LowerData(Data):
    counter = 0
    while counter < len(Data):
        Data[counter][0] = Data[counter][0].lower()
        Data[counter][1] = Data[counter][1].lower()
        counter = counter + 1
    return Data

def MySampleLabel(num):
    SampleLabel = []
    counter = 0
    while counter < num:
        SampleLabel.append(1)
        counter = counter + 1
    counter1 = 0
    while counter1 < num:
        SampleLabel.append(0)
        counter1 = counter1 + 1
    return SampleLabel

# 生成正样本
def PositiveGenerate(DiseaseAndRNABinaryOld, RNAGaussianOld, DiseaseGaussianOld):
    PositiveFeature = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 1:
                pair = []
                pair.extend(RNAGaussianOld[counter1])
                pair.extend(DiseaseGaussianOld[counter])
                PositiveFeature.append(pair)
            counter1 = counter1 + 1
        counter = counter + 1
    return PositiveFeature

# 生成候选负样本
def NegativeCandidateGenerate(DiseaseAndRNABinaryOld, RNAGaussianOld, DiseaseGaussianOld):
    NegativeFeatureAll = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 0:
                pair = []
                pair.extend(RNAGaussianOld[counter1])
                pair.extend(DiseaseGaussianOld[counter])
                NegativeFeatureAll.append(pair)
            counter1 = counter1 + 1
        counter = counter + 1
    return NegativeFeatureAll

def NegativeGenerate(RNAFeatureDAG, DiseaseFeatureDAG):
    # 负样本为全部的disease-rna（328*881）中随机抽取，未在内LncDisease即为负样本
    LncDisease = []
    ReadMyCsv2(LncDisease, 'LncDisease.csv')
    AllDisease = []
    ReadMyCsv2(AllDisease, 'AllDisease.csv')
    AllRNA = []
    ReadMyCsv2(AllRNA, 'AllRNA.csv')
    import random
    NegativeSample = []
    NegativeSampleFeature = []
    counterN = 0
    while counterN < len(LncDisease):  # 随机选出一个疾病rna对
        counterD = random.randint(0, len(AllDisease) - 1)
        counterR = random.randint(0, len(AllRNA) - 1)
        DiseaseAndRnaPair = []
        DiseaseAndRnaPair.append(AllRNA[counterR])
        DiseaseAndRnaPair.append(AllDisease[counterD])
        flag1 = 0
        counter = 0
        while counter < len(LncDisease):
            if DiseaseAndRnaPair == LncDisease[counter]:
                flag1 = 1
                break
            counter = counter + 1
        if flag1 == 1:
            continue
        flag2 = 0
        counter1 = 0
        while counter1 < len(NegativeSample):  # 在已选的负样本中没有，防止重复
            if DiseaseAndRnaPair == NegativeSample[counter1]:
                flag2 = 1
                break
            counter1 = counter1 + 1
        if flag2 == 1:
            continue
        if (flag1 == 0 & flag2 == 0):
            NamePair = []  # 生成对
            NamePair.append(AllRNA[counterR][0])
            NamePair.append(AllDisease[counterD][0])
            NegativeSample.append(NamePair)

            FeaturePair0 = []  # 生成Feature NMF
            FeaturePair0.extend(RNAFeatureDAG[counterR])
            FeaturePair0.extend(DiseaseFeatureDAG[counterD])
            NegativeSampleFeature.append(FeaturePair0)

            counterN = counterN + 1
    return NegativeSampleFeature, NegativeSample

def NegativeGenerateCaseStudy(RNAFeatureDAG, DiseaseFeatureDAG,DN):
    # 负样本为全部的disease-rna（328*881）中随机抽取，未在内LncDisease即为负样本
    LncDisease = []
    ReadMyCsv2(LncDisease, 'LncDisease.csv')
    AllDisease = []
    ReadMyCsv2(AllDisease, 'AllDisease.csv')
    AllRNA = []
    ReadMyCsv2(AllRNA, 'AllRNA.csv')
    import random
    NegativeSample = []
    NegativeSampleFeature = []
    counterN = 0
    while counterN < len(LncDisease):  # 随机选出一个疾病rna对
        counterD = random.randint(0, len(AllDisease) - 1)
        counterR = random.randint(0, len(AllRNA) - 1)
        if AllDisease[counterD][0] == DN:
            continue
        DiseaseAndRnaPair = []
        DiseaseAndRnaPair.append(AllRNA[counterR])
        DiseaseAndRnaPair.append(AllDisease[counterD])
        flag1 = 0
        counter = 0
        while counter < len(LncDisease):
            if DiseaseAndRnaPair == LncDisease[counter]:
                flag1 = 1
                break
            counter = counter + 1
        if flag1 == 1:
            continue
        flag2 = 0
        counter1 = 0
        while counter1 < len(NegativeSample):  # 在已选的负样本中没有，防止重复
            if DiseaseAndRnaPair == NegativeSample[counter1]:
                flag2 = 1
                break
            counter1 = counter1 + 1
        if flag2 == 1:
            continue
        if (flag1 == 0 & flag2 == 0):
            NamePair = []  # 生成对
            NamePair.append(AllRNA[counterR][0])
            NamePair.append(AllDisease[counterD][0])
            NegativeSample.append(NamePair)

            FeaturePair0 = []  # 生成Feature NMF
            FeaturePair0.extend(RNAFeatureDAG[counterR])
            FeaturePair0.extend(DiseaseFeatureDAG[counterD])
            NegativeSampleFeature.append(FeaturePair0)

            counterN = counterN + 1
    return NegativeSampleFeature, NegativeSample

def NegativeGenerateCaseStudy2(RNAFeatureDAG, DiseaseFeatureDAG,DN,num):
    # 负样本为全部的disease-rna（328*881）中随机抽取，未在内LncDisease即为负样本
    LncDisease = []
    ReadMyCsv2(LncDisease, 'LncDisease.csv')
    AllDisease = []
    ReadMyCsv2(AllDisease, 'AllDisease.csv')
    AllRNA = []
    ReadMyCsv2(AllRNA, 'AllRNA.csv')
    import random
    NegativeSample = []
    NegativeSampleFeature = []
    counterN = 0
    while counterN < num:  # 随机选出一个疾病rna对
        counterD = random.randint(0, len(AllDisease) - 1)
        counterR = random.randint(0, len(AllRNA) - 1)
        if AllDisease[counterD][0] == DN:
            continue
        DiseaseAndRnaPair = []
        DiseaseAndRnaPair.append(AllRNA[counterR])
        DiseaseAndRnaPair.append(AllDisease[counterD])
        flag1 = 0
        counter = 0
        while counter < len(LncDisease):
            if DiseaseAndRnaPair == LncDisease[counter]:
                flag1 = 1
                break
            counter = counter + 1
        if flag1 == 1:
            continue
        flag2 = 0
        counter1 = 0
        while counter1 < len(NegativeSample):  # 在已选的负样本中没有，防止重复
            if DiseaseAndRnaPair == NegativeSample[counter1]:
                flag2 = 1
                break
            counter1 = counter1 + 1
        if flag2 == 1:
            continue
        if (flag1 == 0 & flag2 == 0):
            NamePair = []  # 生成对
            NamePair.append(AllRNA[counterR][0])
            NamePair.append(AllDisease[counterD][0])
            NegativeSample.append(NamePair)

            FeaturePair0 = []  # 生成Feature NMF
            FeaturePair0.extend(RNAFeatureDAG[counterR])
            FeaturePair0.extend(DiseaseFeatureDAG[counterD])
            NegativeSampleFeature.append(FeaturePair0)

            counterN = counterN + 1
    return NegativeSampleFeature, NegativeSample

# 利用IsoForest生成强负样本，利用IsoForest取10%异常点，利用打分值选出全局SN
def StrongNegativeGenerate(DiseaseAndRNABinaryOld, RNAGaussianOld, DiseaseGaussianOld,LncRNADiseaseAssociationOld):
    # 生成正样本和所有未标记样本
    print('# 生成正样本和所有未标记样本')
    NegativeFeatureAll = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 0:
                pairFeature = []
                pairFeature.extend(RNAGaussianOld[counter1])
                pairFeature.extend(DiseaseGaussianOld[counter])
                NegativeFeatureAll.append(pairFeature)
            counter1 = counter1 + 1
        counter = counter + 1
    # IsoForest为所有未标记样本赋权值
    print('# IsoForest为所有未标记样本赋权值')
    from sklearn.ensemble import IsolationForest
    clf = IsolationForest(contamination=0.1)
    clf.fit(NegativeFeatureAll)
    scores_pred = clf.decision_function(NegativeFeatureAll)
    # 增加序号标签
    PredictionScoreNum = []
    counter = 0
    while counter < len(scores_pred):
        pair = []
        pair.append(scores_pred[counter])
        pair.append(counter)
        PredictionScoreNum.append(pair)
        counter = counter + 1
    # 选出得分最高的前len(LncRNADiseaseAssociationOld)个作为强负样本
    print('# 选出得分最高的前len(LncRNADiseaseAssociationOld)个作为强负样本')
    SerialNumber = 0
    MaxScoreNum = []
    counter = 0
    while counter < len(LncRNADiseaseAssociationOld):
        max = PredictionScoreNum[0][0]
        counter1 = 0
        while counter1 < len(PredictionScoreNum):
            if max < PredictionScoreNum[counter1][0]:
                max = PredictionScoreNum[counter1][0]
                SerialNumber = counter1
            counter1 = counter1 + 1
        MaxScoreNum.append(PredictionScoreNum[SerialNumber][1])
        del PredictionScoreNum[SerialNumber]
        counter = counter + 1
        # print(counter)
    # 生成负样本NegativeFeature
    print('# 生成负样本NegativeFeature')
    NegativeFeature = []
    counter = 0
    while counter < len(MaxScoreNum):
        NegativeFeature.append(NegativeFeatureAll[MaxScoreNum[counter]])
        counter = counter + 1
    return NegativeFeature, NegativeFeatureAll

# 计算正确的个数，RMSE和MAE
def MyEvaluate(prediction, prediction_proba,TestSample):
    import math
    num = 0
    SumRMSE = 0
    SumMAE = 0
    counter = 0
    while counter < len(prediction):
        SumRMSE = SumRMSE + math.pow((1 - prediction_proba[counter][1]), 2)
        SumMAE = SumMAE + abs(1 - prediction_proba[counter][1])
        if prediction[counter] == 1:
            num = num + 1
        counter = counter + 1
    RMSE = math.sqrt(SumRMSE / len(TestSample))
    MAE = SumMAE / len(TestSample)
    print('TrueNum ?/243: ', num)
    print('RMSE:', RMSE)
    print('MAE:', MAE)
    MyResult = []
    MyResult.append(num)
    MyResult.append(RMSE)
    MyResult.append(MAE)
    return MyResult

# 预测TestSample
def MyPrediction(SampleFeature,SampleLabel,TestSample):
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=100)
    model.fit(SampleFeature, SampleLabel)
    prediction = model.predict(TestSample)
    prediction_proba = model.predict_proba(TestSample)
    print('RandomForestClassifier!')
    result = MyEvaluate(prediction, prediction_proba, TestSample)
    return result

# 预测全部未标记样本的概率值并填补矩阵
def MyPredictionAndMatrixCompletion(SampleFeature,SampleLabel,NegativeFeatureAll,DiseaseAndRNABinaryOld1,DiseaseAndRNABinaryOld2,TestSample):
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=100)
    model.fit(SampleFeature, SampleLabel)
    prediction = model.predict(TestSample)
    prediction_proba = model.predict_proba(TestSample)
    print('RandomForestClassifier!')
    result = MyEvaluate(prediction, prediction_proba, TestSample)
    # 填补矩阵
    prediction_proba_all = model.predict_proba(NegativeFeatureAll)
    num = 0
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld1):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld1[counter]):
            if DiseaseAndRNABinaryOld1[counter][counter1] == 0:
                DiseaseAndRNABinaryOld2[counter][counter1] = prediction_proba_all[num][1]
                num = num + 1
            counter1 = counter1 + 1
        counter = counter + 1
    return DiseaseAndRNABinaryOld2, result

def DiseaseGaussianKernel(DiseaseAndRNABinary):
    # 计算rd
    counter1 = 0
    sum1 = 0
    while counter1 < (len(DiseaseAndRNABinary)):
        counter2 = 0
        while counter2 < (len(DiseaseAndRNABinary[counter1])):
            sum1 = sum1 + pow((DiseaseAndRNABinary[counter1][counter2]), 2)
            counter2 = counter2 + 1
        counter1 = counter1 + 1
    # print('sum1=', sum1)
    Ak = sum1
    Nd = len(DiseaseAndRNABinary)
    rdpie = 0.5
    rd = rdpie * Nd / Ak
    # print('disease rd', rd)
    # 生成DiseaseGaussian
    DiseaseGaussian = []
    counter1 = 0
    while counter1 < len(DiseaseAndRNABinary):  # 计算疾病counter1和counter2之间的similarity
        counter2 = 0
        DiseaseGaussianRow = []
        while counter2 < len(DiseaseAndRNABinary):  # 计算Ai*和Bj*
            AiMinusBj = 0
            sum2 = 0
            counter3 = 0
            AsimilarityB = 0
            while counter3 < len(DiseaseAndRNABinary[counter2]):  # 疾病的每个属性分量
                sum2 = pow((DiseaseAndRNABinary[counter1][counter3] - DiseaseAndRNABinary[counter2][counter3]),2)  # 计算平方
                AiMinusBj = AiMinusBj + sum2
                counter3 = counter3 + 1
            AsimilarityB = math.exp(- (AiMinusBj / rd))
            DiseaseGaussianRow.append(AsimilarityB)
            counter2 = counter2 + 1
        DiseaseGaussian.append(DiseaseGaussianRow)
        counter1 = counter1 + 1
        # print(counter1)
    return DiseaseGaussian

def RNAGaussianKernel(DiseaseAndRNABinary):
    MDiseaseAndRNABinary = np.array(DiseaseAndRNABinary)  # 列表转为矩阵
    RNAAndDiseaseBinary = MDiseaseAndRNABinary.T  # 转置DiseaseAndMiRNABinary
    RNAGaussian = []
    counter1 = 0
    sum1 = 0
    while counter1 < (len(RNAAndDiseaseBinary)):  # rna数量
        counter2 = 0
        while counter2 < (len(RNAAndDiseaseBinary[counter1])):  # disease数量
            sum1 = sum1 + pow((RNAAndDiseaseBinary[counter1][counter2]), 2)
            counter2 = counter2 + 1
        counter1 = counter1 + 1
    # print('sum1=', sum1)
    Ak = sum1
    Nm = len(RNAAndDiseaseBinary)
    rdpie = 0.5
    rd = rdpie * Nm / Ak
    # print('RNA rd', rd)
    # 生成RNAGaussian
    counter1 = 0
    while counter1 < len(RNAAndDiseaseBinary):  # 计算rna counter1和counter2之间的similarity
        counter2 = 0
        RNAGaussianRow = []
        while counter2 < len(RNAAndDiseaseBinary):  # 计算Ai*和Bj*
            AiMinusBj = 0
            sum2 = 0
            counter3 = 0
            AsimilarityB = 0
            while counter3 < len(RNAAndDiseaseBinary[counter2]):  # rna的每个属性分量
                sum2 = pow((RNAAndDiseaseBinary[counter1][counter3] - RNAAndDiseaseBinary[counter2][counter3]),2)  # 计算平方，有问题？？？？？
                AiMinusBj = AiMinusBj + sum2
                counter3 = counter3 + 1
            AsimilarityB = math.exp(- (AiMinusBj / rd))
            RNAGaussianRow.append(AsimilarityB)
            counter2 = counter2 + 1
        RNAGaussian.append(RNAGaussianRow)
        counter1 = counter1 + 1
        # print(counter1)
    return RNAGaussian

# 产生测试样本
def TestSampleFeatureGenerate(LncRNADiseaseAssociationNew, AllDiseaseOld, AllRNAOld, RNAGaussianOld, DiseaseGaussianOld, DiseaseAndRNABinaryOld):
    # 产生Test样本
    ExtraPairNum = []
    ExtraPairName = []
    TestSampleFeature = []
    counter = 0
    while counter < len(LncRNADiseaseAssociationNew):
        rna = LncRNADiseaseAssociationNew[counter][0]
        disease = LncRNADiseaseAssociationNew[counter][1]
        counter1 = 0
        while counter1 < len(AllDiseaseOld):
            if disease == AllDiseaseOld[counter1]:
                counter2 = 0
                while counter2 < len(AllRNAOld):
                    if rna == AllRNAOld[counter2]:
                        if DiseaseAndRNABinaryOld[counter1][counter2] == 0:
                            pairNum = []
                            pairNum.append(counter1)
                            pairNum.append(counter2)
                            ExtraPairNum.append(pairNum)
                            pairName = []
                            pairName.append(AllDiseaseOld[counter1])
                            pairName.append(AllRNAOld[counter2])
                            ExtraPairName.append(pairName)
                            pairFeature = []
                            pairFeature.extend(RNAGaussianOld[counter2])
                            pairFeature.extend(DiseaseGaussianOld[counter1])
                            TestSampleFeature.append(pairFeature)
                        break
                    counter2 = counter2 + 1
                break
            counter1 = counter1 + 1
        counter = counter + 1
    return TestSampleFeature