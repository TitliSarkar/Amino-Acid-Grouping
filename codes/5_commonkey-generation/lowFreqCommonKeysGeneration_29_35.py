# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 18:40:07 2018

@author: Titli
"""
import pandas as pd
import glob

excelfile ='Protein ID for TSR Manuscript 8-24-2018.xlsx' # change for diff DATASET    
df= pd.read_excel(excelfile,sheet_name='PKABC')
group_dict = df.groupby('Name')['PDB IDs'].apply(list).to_dict()
print (group_dict)

# 44 PKABC
files = ['1BX6','1CDK','1CMK','1J3H','1L3R','2C1A','2CPK','3TNQ','3NX8','4UJ1','4YXS','5O5M','5UZK','5VHB','4XW5','1ATP','1MRV','1GZK','1O6K','2JDO','2UW9','2X39','3D0E','3E87','6CCY','4GV1','3OCB','3QKK','4EKK','3CQU','3MV5','4EJN','5KCV','3O96','1GZO','2FK9','1A25','1GMI','1DSY','2NCE','4L1L','5W4S','3GPE','3RDJ']
#files= ['1BX6']

def generate_lowFreqCommonTriplets(files):      
    for prot in files:
        print(prot)
        with open(prot+'.commonkeys_29_35','r') as f:
            lowfreqkeys = [] # first, get low freq keys in a list (freq=1 or freq<=5)
            for line in f:
                if (int(line.split()[1]) == 1): #freq =5 
                    #print (line)
                    lowfreqkeys.append(line.split()[0])
            i=0
            with open(prot+'.commonTriplets_29_35','r') as f: # then only keep keys which has freq=1 and maxdist<=15
                #f_out = open(prot+'.lowFreqCommonTriplets_maxDist15_29_35_freq1','w')
                f_out = open(prot+'.lowFreqCommonTriplets_29_35_freq5','w') # removed maxdist<=15 constraint
                #f_out = open(prot+'.lowFreqCommonTriplets_29_35_maxdist15','w') # removed freq=1 constraint
    
                for line in f:
                    #if (line.split()[0] in lowfreqkeys and float(line.split()[10])<=15):
                    if (line.split()[0] in lowfreqkeys): # removed maxdist<=15 constraint
                    #if (float(line.split()[10])<=15.0):# removed freq=1 constraint
                        #print(line)
                        f_out.write(line)
                        i += 1
            print(prot, len(lowfreqkeys), i)

def add_freq_to_keys(files,ext): # add freq occurence no to each key
    for i in files:
        if i != '1BX6':
            continue
        df = pd.read_csv(i+ext, sep='\t', names = ['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
        #print(df)
        #df['key'] = df['key'].astype(str) + '_' + df.groupby('key')['key'].transform('count').astype(str)
        print(df['key'])
        #df.to_csv(str(i+ext+'_freq_added'), header=None, index=None, sep='\t')
        print (i)       

# generate low freq common triplets from freq_added files        
def key_filtering(files,ext):
    print("----------------------------------------")
    for prot in files:
        df = pd.read_csv(prot+ext, sep='\t', names = ['key','amino0','pos0','amino1','pos1','amino2','pos2','classT1','theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'])
        #print(df)
        df['key_freq'] = df.groupby('key')['key'].transform('count') # correct
        
        df1 = df[(df['key_freq'] == 1)& (df['maxDist']<=15)] #& (df['maxDist']<=12)
        
        df1['key']=df1['key'].astype(str)+'_'+df1['key_freq'].astype(str)
        df2 = df1.drop(columns=['key_freq'])
        print(prot, len(df2))
        #print(df1)
        df2.to_csv(str(prot+ext+'_freq_added_freq1_maxdist15'), index=None, header=None, sep='\t', mode='a')

#add_freq_to_keys(files,ext)
key_filtering(files,'.commonTriplets_29_35')
#add_freq_to_keys(files,'.commonTriplets_29_35_freq1_maxdist12')

# =============================================================================
# for k,v in group_dict.items():
#     fileList = v
#     ext= '.commonTriplets_'+k+'_notOthers_29_35'
#     add_freq_to_keys(fileList,ext)
# =============================================================================
