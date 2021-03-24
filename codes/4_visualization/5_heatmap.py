# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:16:16 2019

@author: Titli
""" 
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import argparse

# =============================================================================
# import seaborn as sns; sns.set(color_codes=True)
# iris = sns.load_dataset("iris")
# species = iris.pop("species")
# print(iris)
# g = sns.clustermap(iris)
# 
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("path", help="enter path to the input file", type=str) # ./../Dataset/Intra_data/
args = parser.parse_args()

data_dir = args.path # input files location

data_location = len(data_dir.split("/"))-1
print(data_dir.split("/"))
print(data_location)

df=pd.read_csv(data_dir+'JaccardDistance.txt', header=0, sep=",")
df=df.drop('0', axis=1)

uniques = df.apply(lambda x: x.nunique())
drop_index = str(uniques[uniques==1].index[0]).strip()
print(drop_index, type(drop_index))
df = df.drop(uniques[uniques==1].index, axis=1)
#print(df)
df1=df.set_index('prot')
df2 = df1.drop(drop_index)
print(df2)
df2.to_csv(data_dir+'df.csv', sep=',')
files=list(df2.columns)
print(len(files), files)

sns.set(font_scale=0.8)
sns_plot = sns.clustermap(df2, method="average",xticklabels=True, yticklabels=True, figsize=(50,50))
sns_plot.savefig(data_dir+'heatmap.jpg')


