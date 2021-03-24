# -*- coding: utf-8 -*-

import pandas as pd
import multiprocessing
from joblib import Parallel, delayed
import os,glob

# Step 2
#finds common keys among given set of proteins
def find_common_keys(inputProts,ext):
    keys_dict = {} # key=ProtKey, val=Protname
    
    for p in inputProts: 
        with open(p+".keys_"+ext, 'r') as file: # .keys_29_35
            for r in file: # row is list
                if not r:
                    continue
                row = r.split()
                #print (row[0], row[1])
                if row[0] not in keys_dict.keys():   # row[0] = prot, row[1] = key
                    keys_dict[row[0]] = []
                    keys_dict[row[0]].append(p)
                else:
                    keys_dict[row[0]].append(p)

    # remove duplicate proteins coming from multiple instances of a key within same
    for key_, val_ in keys_dict.items():
        val_ = list(set(val_))
        keys_dict[key_]=val_
    for k,v in keys_dict.items():
        if(len(v) == len(inputProts)):
            print(k,len(v))
    
    # remove entries with not common keys
    keys_dict = { k:v for k,v in keys_dict.items() if (len(v)==len(inputProts) and sorted(v)==sorted(inputProts))} 
    
    # find common keys
    common_keys = list(keys_dict.keys())
    print(ext, " #of unq common keys= ",len(common_keys))
    return common_keys

# Step 2: find common key details
def parallel_code_triplets(p, commonKeys,ext):
    print(p," triplets")
    result_file = open(p+'.commonTriplets_'+ext,'w') #.commonTriplets_29_35

    with open(p+'.triplets_'+ext, 'r') as file: # .triplets_29_35
            for row in file:
                if row.split()[0] in commonKeys:
                    result_file.writelines(row)
                    
def parallel_code_keys(p, commonKeys,ext):
    print(p, " keys")
    result_file1 = open(p+'.commonKeys_'+ext,'w') #.commonKeys_29_35

    with open(p+'.keys_'+ext, 'r') as file1: # .keys_29_35
            for row1 in file1:
                if row1.split()[0] in commonKeys:
                    result_file1.writelines(row1)
                    
def find_commonkey_details(inputProts, commonKeys,ext):
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores, verbose=50)(delayed(parallel_code_triplets)(fileName,commonKeys,ext)for fileName in inputProts)
    Parallel(n_jobs=num_cores, verbose=50)(delayed(parallel_code_keys)(fileName,commonKeys,ext)for fileName in inputProts)

    
#func calls
extensions = ['29_35'] 
for ext in extensions:
    fileList=[]
    for i in glob.glob('*.keys_'+ext): #'*.keys_29_35'
        print(i) 
        fileList.append(i.split(".")[0])#.split("/")[5]) #[5]
    print(fileList, len(fileList))
    common_keys = find_common_keys(fileList, ext)  
    print (ext, " #of unq common keys: ", len(common_keys))
    find_commonkey_details(fileList,common_keys,ext)
    
