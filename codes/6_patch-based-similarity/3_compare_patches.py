# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 20:56:46 2020

@author: titli

### comment -  Please ignore columns for jaccard_sim_plus_or_minus 1 and 2 in the files generated for freq, they are not calculated yet
"""

## Code for pairwise similarity calculation between two proteins using allxall patches across two proteins
## Input : all patch files, triplet files
## Output: pairwise patch similarity between two proteins

import pickle
import pandas as pd
import os, time, glob
import numpy as np
from collections import Counter
import argparse

parser = argparse.ArgumentParser()
# input files 
parser.add_argument("prot1", help="enter protein1", type=str) 
parser.add_argument("prot2", help="enter protein2", type=str) 
parser.add_argument("chain", help="enter chain", type=str) 
parser.add_argument("patch_radius", help="enter patch radius", type=int) 
args = parser.parse_args()

prot1 = args.prot1 #'6W75'
prot2 = args.prot2 #'3P9I'
chain = args.chain #'A'  

patch_radius = args.patch_radius #15 # max dist of each neighbour from patch center
patch_neighbours = 300 # max no. of points allowed in patch neighbourhood

data_dir = './../Dataset_'+prot1+chain+'_'+prot2+chain+'_all/'
print(data_dir)
data_location = len(data_dir.split("/"))-1
print(data_dir.split("/"))
print(data_location)

subdir = 'aa_grp_0/'
os.makedirs(data_dir+subdir, exist_ok=True)
files = [x.split("/")[data_location].split(".")[0] for x in glob.glob(data_dir+"*.pdb")]
print(len(files), files)

## compare patches
#prot = pickle.load(open(data_dir+subdir+'1ece_surface.patches', 'rb'))
#print(prot.keys())
#actual_binding = [114,162,238,116,161,27]
#intersection = set(prot.keys()) & set(actual_binding)
#print(intersection)

## claculate local similairy between two patches in terms of patch keys
## jaccard similarity 
def patch_similarity_jaccard(x, y):
   intersection_cardinality = len(set.intersection(*[set(x), set(y)]))
   union_cardinality = len(set.union(*[set(x), set(y)]))
   if union_cardinality!=0: ## avoid division by zero error
       return intersection_cardinality/float(union_cardinality)
   else:
       union_cardinality=1
       return intersection_cardinality/float(union_cardinality)

def patch_similarity_elementwise(x, y): 
    print(x,y)
    intersection = set(x) & set(y)
    print(intersection)
    return 0

def patch_similarity_plusorminus1(x, y): 
    intersection = []
    for i in x: # for each key in list 1
        range_ = [a for a in range(i-1, i+1+1)]
        if set(range_) & set(y) != 0: # if key i in list 1 have a +/-2 in list 2, include i in intersection
            print(i, set(range_) & set(y))    
            intersection.append(i)
            print(intersection)
    print(sorted(intersection), len(intersection))
    return 0

def patch_similarity_plusorminus(x, y, amount): 
    range1 = [a for i in x for a in range(i-amount, i+amount+1)]
    range2 = [b for j in y for b in range(j-amount, j+amount+1)]

    intersection1 = list(set(range1) & set(y))
    intersection2 = list(set(range2) & set(x))

    # percentage of match calculation
    total_keys = list(set(x+y))    
    matched_keys = list(set(intersection1+intersection2))    
    percent_match = (len(matched_keys)/len(total_keys))*100
    
    #print("total keys = ", len(total_keys), "matched keys = ", len(matched_keys), "%match = ", percent_match)
    return percent_match

## calculate similarity for allxall patches for the given pair of proteins
## Input: The patch neighbour residues for each file

# helper: generate keys for each patch 
def get_patch_keys(prot1, prot2):
    print(prot1, prot2)
    
    neigh_res_dict1 = pickle.load(open(data_dir+subdir+prot1+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours), "rb"))
    neigh_res_dict2 = pickle.load(open(data_dir+subdir+prot2+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours), "rb"))
    print(prot1, "residues=", len(neigh_res_dict1), prot2, "residues=", len(neigh_res_dict2))
    
    filetriplets1 = pd.read_csv(data_dir+subdir+prot1+'.triplets_29_35_aa_grp_0', sep='\t', \
                               usecols=[0,2,4,6], \
                               names=['key', 'pos0', 'pos1', 'pos2'])
    filetriplets2 = pd.read_csv(data_dir+subdir+prot2+'.triplets_29_35_aa_grp_0', sep='\t', \
                               usecols=[0,2,4,6], \
                               names=['key', 'pos0', 'pos1', 'pos2'])
    
    # extract keys related to neighbours, for each patch
    neigh_keys_dict1 = {}
    neigh_keys_dict2 = {}
    for (k1,v1) in neigh_res_dict1.items(): # k = patch center
        #print(k1)
        neigh_resi1 = [int(x[0]) for x in v1] # find neighbouting residues
        neigh_resi1.append(k1)  # add patch center to neighbouring residue list
        neigh_keys1 = filetriplets1[(filetriplets1['pos0'].isin(neigh_resi1)) & \
                                  (filetriplets1['pos1'].isin(neigh_resi1)) & \
                                  (filetriplets1['pos2'].isin(neigh_resi1))]['key'].tolist()
        neigh_keys_dict1[k1] = neigh_keys1
        
    for (k2,v2) in neigh_res_dict2.items(): # k = patch center
        #print(k2)
        neigh_resi2 = [int(x[0]) for x in v2] # find neighbouting residues
        neigh_resi2.append(k2)  # add patch center to neighbouring residue list
        neigh_keys2 = filetriplets2[(filetriplets2['pos0'].isin(neigh_resi2)) & \
                                  (filetriplets2['pos1'].isin(neigh_resi2)) & \
                                  (filetriplets2['pos2'].isin(neigh_resi2))]['key'].tolist()
        neigh_keys_dict2[k2] = neigh_keys2
        
    print("key dict len = ", len(neigh_keys_dict1), len(neigh_keys_dict2))
    
    # save neighbour_keys dict
    with open(data_dir+subdir+prot1+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "wb") as f:
        pickle.dump(neigh_keys_dict1, f)
    f.close()
    with open(data_dir+subdir+prot2+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "wb") as f:
        pickle.dump(neigh_keys_dict2, f)
    f.close()


def patch_similarity_calculation(prot1, prot2):
    print(prot1, prot2) 
    # get_patch_keys(prot1, prot2) ## func calls
    
    resi_name_dict1 = pickle.load(open(data_dir+subdir+prot1+'.residuenames', 'rb'))
    resi_name_dict2 = pickle.load(open(data_dir+subdir+prot2+'.residuenames', 'rb'))
    
    # open output file
    if prot1 == prot2:
        sim_out = open(data_dir+subdir+prot1+'_Similarity_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
        dist_out = open(data_dir+subdir+prot1+'_Distance_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
    else:
        sim_out = open(data_dir+subdir+prot1+'_'+prot2+'_Similarity_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
        dist_out = open(data_dir+subdir+prot1+'_'+prot2+'_Distance_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
     
    # write headers
    sim_out.writelines(['prot1', '\t','patchi_id', '\t', 'patchi', '\t', 'prot2', '\t', 'patchj_id', '\t', 'patchj', '\t', 'sim_plus_or_minus_0', '\t', 'sim_plus_or_minus_1', '\t', 'sim_plus_or_minus_2', '\n'])
    dist_out.writelines(['prot1', '\t','patchi_id', '\t', 'patchi', '\t', 'prot2', '\t', 'patchj_id', '\t', 'patchj', '\t', 'dist_plus_or_minus_0', '\t', 'dist_plus_or_minus_1', '\t', 'dist_plus_or_minus_2', '\n'])

    # Load patch keys 
    neigh_keys_dict1 = pickle.load(open(data_dir+subdir+prot1+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "rb"))
    neigh_keys_dict2 = pickle.load(open(data_dir+subdir+prot2+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "rb"))
    print("patch keys = ", len(neigh_keys_dict1), len(neigh_keys_dict2))

    # calculate patch similairty and write to a file    
    for k1,_ in neigh_keys_dict1.items():
        for k2,_ in neigh_keys_dict2.items():
            #print(k1,k2)
            v1 = sorted(neigh_keys_dict1[k1])
            v2 = sorted(neigh_keys_dict2[k2])
            
            sim_pm0_ = np.round(patch_similarity_plusorminus(v1, v2, 0),2) # calculate similarity +/- 0
            dist_pm0_ = np.round(100-sim_pm0_,2) # calculate distance +/- 0
            
            sim_pm1_ = np.round(patch_similarity_plusorminus(v1, v2, 1),2) # calculate similarity +/- 1
            dist_pm1_ = np.round(100-sim_pm1_,2) # calculate distance +/- 1
            
            sim_pm2_ = np.round(patch_similarity_plusorminus(v1, v2, 2),2) # calculate similarity +/- 2
            dist_pm2_ = np.round(100-sim_pm2_,2) # calculate distance +/- 2
            #sim_list.append(sim_)
            
            print(k1, k2, sim_pm0_, dist_pm0_, sim_pm1_, dist_pm1_, sim_pm2_, dist_pm2_)
            sim_out.writelines([prot1, '\t',str(k1), '\t', resi_name_dict1[k1], '\t', prot2, '\t', str(k2), '\t', resi_name_dict2[k2], '\t', str(sim_pm0_), '\t', str(sim_pm1_), '\t', str(sim_pm2_), '\n'])
            dist_out.writelines([prot1, '\t',str(k1), '\t', resi_name_dict1[k1], '\t', prot2, '\t', str(k2), '\t', resi_name_dict2[k2], '\t', str(dist_pm0_), '\t', str(dist_pm1_), '\t', str(dist_pm2_), '\n'])
    
    sim_out.close()
    dist_out.close()
    # end of func 

def patch_similarity_plusorminus_genjaccard(x, y, amount):
    #print(len(x), len(y))
    #print (len(set(x)), len(set(y)))
    
    dict1 = dict(Counter(x))
    dict2 = dict(Counter(y))
    total_unq_keys = list(set(x+y)) # get total unique keys
    #print(dict1, "\n", dict2)
    # Calculate Generalized Jaccard 
    numerator_gen_jac = 0
    denomenator_gen_jac = 0
    for key_ in total_unq_keys:
        a = 0
        b = 0
        if key_ in dict1:
            a = dict1[key_]
        if key_ in dict2:
            b = dict2[key_]
        numerator_gen_jac += min(a,b)
        denomenator_gen_jac += max(a,b)
        #print("key, a ,b   nu, deno= ", key_, a ,b, numerator_gen_jac, denomenator_gen_jac)

    #print("end num, deno = ", numerator_gen_jac, denomenator_gen_jac)
    # print denomenator_jac
    similarity = 0.0
    if denomenator_gen_jac == 0:
        print ("denominator is 0: ", numerator_gen_jac, denomenator_gen_jac, x, y)
    else:
        similarity = np.round(float(numerator_gen_jac)/float(denomenator_gen_jac), 2)
    return similarity
    
## compare patches with considering key frequency in similarity calculation formula
def patch_similarity_calculation_freq(prot1, prot2):
    print(prot1, prot2)
    
    resi_name_dict1 = pickle.load(open(data_dir+subdir+prot1+'.residuenames', 'rb'))
    resi_name_dict2 = pickle.load(open(data_dir+subdir+prot2+'.residuenames', 'rb'))
    
    # open output file
    sim_out = open(data_dir+subdir+prot1+'_'+prot2+'_Similarity_withFreq_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
    dist_out = open(data_dir+subdir+prot1+'_'+prot2+'_Distance_with_Freq_paches_'+chain+'_radius'+str(patch_radius)+'_neighbours'+str(patch_neighbours)+'.txt','w')
    
    # write headers
    sim_out.writelines(['prot1', '\t','patchi_id', '\t', 'patchi', '\t', 'prot2', '\t', 'patchj_id', '\t', 'patchj', '\t', 'sim_genjacc_plus_or_minus_0', '\t', 'sim_genjacc_plus_or_minus_1', '\t', 'sim_genjacc_plus_or_minus_2', '\n'])
    dist_out.writelines(['prot1', '\t','patchi_id', '\t', 'patchi', '\t', 'prot2', '\t', 'patchj_id', '\t', 'patchj', '\t', 'dist_genjacc_plus_or_minus_0', '\t', 'dist_genjacc_plus_or_minus_1', '\t', 'dist_genjacc_plus_or_minus_2', '\n'])
    
    # Load patch keys 
    neigh_keys_dict1 = pickle.load(open(data_dir+subdir+prot1+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "rb"))
    neigh_keys_dict2 = pickle.load(open(data_dir+subdir+prot2+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours)+".patchKeys", "rb"))
    print(" each prot has residues = ", len(neigh_keys_dict1), len(neigh_keys_dict2))
    for k1,_ in neigh_keys_dict1.items(): # k1 =  patch center from prot 1
        for k2,_ in neigh_keys_dict2.items(): # k2 = patch center from prot 2
            v1 = sorted(neigh_keys_dict1[k1]) # v1 =  key list in patch 1
            v2 = sorted(neigh_keys_dict2[k2]) # v2 = key list in patch 2
# =============================================================================
#             if k1==128 and k2==128:
#                 print("k1, v1 = ", k1, v1)
#                 print("k2, v2 = ", k2, v2)
#                 sim_genjac_0 = patch_similarity_plusorminus_genjaccard(v1, v2, 0)*100
#                 print(sim_genjac_0)
#                 break
# =============================================================================
            ##calculate patch similairty and write to a file    
            ## Generalized Jaccard (with freq)
            sim_genjac_0 = patch_similarity_plusorminus_genjaccard(v1, v2, 0)
            dist_gen_jac_0 = np.round(100.0-sim_genjac_0, 2)
            
            sim_genjac_1 = patch_similarity_plusorminus_genjaccard(v1, v2, 1)
            dist_gen_jac_1 = np.round(100.0-sim_genjac_1, 2)
            
            sim_genjac_2 = patch_similarity_plusorminus_genjaccard(v1, v2, 2)
            dist_gen_jac_2 = np.round(100.0-sim_genjac_2, 2)
            
            print(k1, k2, sim_genjac_0, dist_gen_jac_0, sim_genjac_1, dist_gen_jac_1, sim_genjac_2, dist_gen_jac_2)
            sim_out.writelines([prot1, '\t',str(k1), '\t', resi_name_dict1[k1], '\t', prot2, '\t', str(k2), '\t', resi_name_dict2[k2], '\t', str(sim_genjac_0), '\t', str(sim_genjac_1), '\t', str(sim_genjac_2), '\n'])
            dist_out.writelines([prot1, '\t',str(k1), '\t', resi_name_dict1[k1], '\t', prot2, '\t', str(k2), '\t', resi_name_dict2[k2], '\t', str(dist_gen_jac_0), '\t', str(dist_gen_jac_1), '\t', str(dist_gen_jac_2), '\n'])            
    sim_out.close()
    dist_out.close()
    # end of func 
    
    
## func call
start = time.time()

get_patch_keys(files[0], files[1]) ## func calls for creating file with patches with keys 

patch_similarity_calculation(files[0], files[1])
patch_similarity_calculation_freq(files[0], files[1])

print("Time taken to calculate allxall patch similarity = ", np.round((time.time()-start))/60, 4)

