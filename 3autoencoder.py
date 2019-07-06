import os           # https://blog.csdn.net/marsjhao/article/details/73480859    https://blog.keras.io/building-autoencoders-in-keras.html
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'    # cpu
from keras.layers import Input, Dense
from keras.models import Model
from keras.datasets import mnist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import math
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


# 读数据
PositiveSampleFeature = []
ReadMyCsv(PositiveSampleFeature, "PositiveSampleFeature.csv")
NegativeSampleFeature = []
ReadMyCsv(NegativeSampleFeature, "NegativeSampleFeature.csv")

SampleFeature = []
SampleFeature.extend(PositiveSampleFeature)
SampleFeature.extend(NegativeSampleFeature)


SampleFeature = np.array(SampleFeature)
x = SampleFeature # (3530, 1209)

from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(x, x, test_size=0.2)    # 切分数据集进行训练，用全部数据集x进行“预测”


# 改变数据类型
x_train = x_train.astype('float32') / 1.
x_test = x_test.astype('float32') / 1.


# 变量
encoding_dim = 128
input_img = Input(shape=(len(SampleFeature[0]),))    # 输入维度

# 构建autoencoder
encoded_input = Input(shape=(encoding_dim,))
encoded = Dense(encoding_dim, activation='relu')(input_img)
decoded = Dense(1209, activation='sigmoid')(encoded)

autoencoder = Model(inputs=input_img, outputs=decoded)
decoder_layer = autoencoder.layers[-1]
encoder = Model(inputs=input_img, outputs=encoded)
decoder = Model(inputs=encoded_input, outputs=decoder_layer(encoded_input))

autoencoder.compile(optimizer='adadelta', loss='binary_crossentropy')
autoencoder.fit(x, x, epochs=100, batch_size=128, shuffle=True, validation_data=(x_test, x_test))

# 预测
encoded_imgs = encoder.predict(x)
decoded_imgs = decoder.predict(encoded_imgs)
storFile(encoded_imgs, 'SampleFeatureAuto.csv')
