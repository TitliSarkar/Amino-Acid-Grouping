# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 18:40:07 2020

@author: titli
"""

## First step towards calculating protrin similarity. This code saves all keys prrsent in the dataset in a file
## Input : .keys files
## Output: a text file with all keys in the sample- 'localFeatureSelection_44PKABC_29_35_allkeys.txt'

import math,glob
import numpy

#saves all the counts for each key 
allKeyCountDict={}

# 44 PKABC
filesList = ['1BX6','1CDK','1CMK','1J3H','1L3R','2C1A','2CPK','3TNQ','3NX8','4UJ1','4YXS','5O5M','5UZK','5VHB','4XW5','1ATP','1MRV','1GZK','1O6K','2JDO','2UW9','2X39','3D0E','3E87','6CCY','4GV1','3OCB','3QKK','4EKK','3CQU','3MV5','4EJN','5KCV','3O96','1GZO','2FK9','1A25','1GMI','1DSY','2NCE','4L1L','5W4S','3GPE','3RDJ']

# =============================================================================
# filesList = []
# for file in glob.glob("*.pdb"):
# 	filesList.append(file.split(".")[0])
# 
# =============================================================================
print (len(filesList),filesList)
n=len(filesList)

for i in filesList:
    print(i)
    f1=open(i+'.keys_29_35','r')
    for j in f1:
        j=j.split()
        key_=j[0].rstrip()
        val_=int(j[1].rstrip())
        if key_ in allKeyCountDict:
            allKeyCountDict[key_].append(val_)
        else:
            allKeyCountDict[key_]=[]
            allKeyCountDict[key_].append(val_)

f2_out=open('localFeatureSelection_44PKABC_29_35_allkeys.txt','w')

for keys_ in allKeyCountDict:
    list_=[]
    list_=allKeyCountDict[keys_]
    numOfMatch=len(list_) #number of proteins that have the key
    numOfGap=n-numOfMatch #number of proteins that DO NOT have the key
    mean_=numpy.mean(list_) #average number of times a key occurs in the list of proteins
    mad_=0.0
   
    for i in range(len(list_)):
        mad_=mad_+math.fabs(list_[i]-mean_)
    mad_=float(mad_)/float(numOfMatch)# calculate MAD
 
    #if ((numOfMatch<=math.ceil(n/5) and mad_<=0.01)):# and (numOfMatch>2)):
    f2_out.writelines([str(keys_),'\t',str(mad_),'\n'])
  
print ("Code End.")
