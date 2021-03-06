# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 18:40:07 2020

@author: titli
"""

## Second step towards calculating protrin similarity. This code finds freq of each key in each protein.
## Input : a text file with all keys in the sample- 'localFeatureSelection_44PKABC_29_35_allkeys.txt'
## Output : a csv file 'localFeatureVect_44PKABC_29_35_allkeys.csv' which contains each row as a protein i, each column as a key j and each cell(i,j) as freq of key j in protein i

import copy,glob

# 44 PKABC
filesList = ['1BX6','1CDK','1CMK','1J3H','1L3R','2C1A','2CPK','3TNQ','3NX8','4UJ1','4YXS','5O5M','5UZK','5VHB','4XW5','1ATP','1MRV','1GZK','1O6K','2JDO','2UW9','2X39','3D0E','3E87','6CCY','4GV1','3OCB','3QKK','4EKK','3CQU','3MV5','4EJN','5KCV','3O96','1GZO','2FK9','1A25','1GMI','1DSY','2NCE','4L1L','5W4S','3GPE','3RDJ']

# =============================================================================
# filesList = []
# for file in glob.glob("*.pdb"):
#     filesList.append(file.split(".")[0])
# =============================================================================
print (len(filesList),filesList)

f1_in=open('localFeatureSelection_44PKABC_29_35_allkeys.txt','r')
f1_out=open('localFeatureVect_44PKABC_29_35_allkeys.csv','w')

keyDict={}
for i in f1_in:
    i=i.split()
    keyDict[i[0].rstrip()]=0
print (len(keyDict))

for i in filesList:
    f2_in=open(i+'.keys_29_35','r')
    keyDict1=copy.deepcopy(keyDict)
    
    for j in f2_in:
    	j=j.split()
    	if j[0].rstrip() in keyDict1:
    	    keyDict1[j[0].rstrip()]=j[1].rstrip()
    for k in keyDict1:
        f1_out.writelines([str(keyDict1[k]),','])
    f1_out.writelines(['0','\n'])
    	
print ("Code End.")
