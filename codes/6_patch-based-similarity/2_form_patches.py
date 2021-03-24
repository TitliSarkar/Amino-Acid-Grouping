# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 18:40:07 2020

@author: titli
"""

## Code for forming patches centering each c-alpha in protein chain, with predefined patch radiues and patch neighbours
## Input : .pdb, .residues and .triplets files
## Output: temporary files for neighbouring residues for each patch center

import pickle
import os, glob
import argparse
from scipy.spatial import distance
import numpy as np

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

subdir ='aa_grp_0/'
os.makedirs(data_dir+subdir, exist_ok=True)
files= [x.split("/")[data_location].split(".")[0] for x in glob.glob(data_dir+"*.pdb")]
print(len(files), files)

## form patches
for fileName in files:
    #if fileName=='2kau':
        #continue
    print (fileName)
    inFile=open(data_dir+fileName+'.pdb','r')
    xCord={}
    yCord={}
    zCord={}
    seq_number={}
    residue={}
    counter=0
    
    for i in inFile:
        if ((i[0:6].rstrip()=="END")):
            break
        if(i[0:4].rstrip())=="ATOM" and (i[13:15].rstrip())=="CA" and i[21:22].strip()==chain:
            #print (i)
            xCord[counter]=(float(i[30:38].strip()))
            yCord[counter]=(float(i[38:46].strip()))
            zCord[counter]=(float(i[46:54].strip()))
            seq_number[counter]=str(i[22:27].strip())
            residue[int(seq_number[counter])]=str(i[17:20].strip())
            counter+=1
    #print(residue)
    with open(data_dir+subdir+fileName+".residuenames", "wb") as f:
        pickle.dump(residue, f)
    f.close()
    
    # for each residue in chain, find neighbours with distance < 12 
    neigh_dict={}
    visited=[] # store visited nodes
    for c1,seq1 in seq_number.items():
        #print(c1,seq, xCord[c1], yCord[c1], zCord[c1])
        neigh_list=[]
        for c2,seq2 in seq_number.items():
            dist=distance.euclidean((xCord[c1], yCord[c1], zCord[c1]), (xCord[c2], yCord[c2], zCord[c2]))
            dist=np.round(dist,4)
            if dist !=0 and dist<=patch_radius:
                #print(c1, c2, dist)
                neigh_list.append((int(seq2),float(dist)))
        neigh_list.sort(key=lambda x:x[1]) # neighbours are sorted by distance
        neigh_dict[int(seq1)]=neigh_list[:patch_neighbours] # pick n nearest neighbours
    
    # print each patch and its neighbours
    for k,v in neigh_dict.items():
        #print(k, len(v))
        if k==1:
            print(k, len(v), v)
        
    with open(data_dir+subdir+fileName+".patches_radius"+str(patch_radius)+"_neigh"+str(patch_neighbours), "wb") as f:
        pickle.dump(neigh_dict, f)
    f.close()


            
