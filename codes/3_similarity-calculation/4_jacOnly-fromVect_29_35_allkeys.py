# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 18:40:07 2020

@author: titli
"""

## Third and last step towards calculating protrin similarity. This code calculates pairwise similarity between all proteins.
## Input : a csv file 'localFeatureVect_44PKABC_29_35_allkeys.csv' which contains each row as a protein i, each column as a key j and each cell(i,j) as freq of key j in protein i
## Output: a csv file with each cell(m,n) corresponds to Jaccard Similarity between protein m and n 

import os,os.path
import numpy as np

f2_out=open('JaccardSimilarity_44PKABC_29_35_allkeys.txt','w')
fileNames=open('localFeatureVect_44PKABC_29_35_allkeys.csv','r')
dict_={}

num=0
for i in fileNames:
    #print i
    i=i.split(',')
    dict2={}
    numOfCols=len(i)
    #print numOfCols
    for k in range(numOfCols-1):#leave the last entry out as it is the dummy 0 entry for all instances
        dict2[k]=int(i[k].rstrip())
    dict_[num]=dict2
    num+=1
#print (dict_)
print (num)

for i in range(num):

    col=len(dict_[i])
    #print col
    for j in range(num):
        numerator_jac=0
        denomenator_jac=0
        numerator_gen_jac=0
        denomenator_gen_jac=0
        #print i, j
        for k in range(col):
            #print (dict_[i][k], dict_[j][k])
            a=dict_[i][k]
            b=dict_[j][k]
            if a>0:
                a_jac=1
            else:
                a_jac=0
            if b>0:
                b_jac=1
            else:
                b_jac=0
            numerator_jac+=min(a_jac,b_jac)
            denomenator_jac+=max(a_jac,b_jac)
            numerator_gen_jac+=min(a,b)
            denomenator_gen_jac+=max(a,b)
      # print denomenator_jac
        if denomenator_jac==0:
            print (i, j, numerator_jac, denomenator_jac)
        else:
            similarity = float(numerator_gen_jac)/float(denomenator_gen_jac)
            simRounded = np.round(similarity*100,4)

            dist_gen_jac=1.0-similarity
            #print (similarity)
            f2_out.writelines([str(simRounded),'\t'])

    f2_out.writelines(['\n'])
    
print ("Code End.")
