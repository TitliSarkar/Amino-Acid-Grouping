# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:16:16 2019

@author: Titli
""" 
import os
import plotly
import plotly.plotly as py
import plotly.figure_factory as ff
import numpy as np

from scipy.cluster.hierarchy import linkage

#path = '/home/C00222141/Desktop/134Proteins/'
#os.chdir(path)
print(os.getcwd())
print(plotly.__version__)

# before doing this, loging to plotty @ https://plot.ly with email: xxx@xxx.edu username:xxx pwd:xxx and generate API_key
plotly.tools.set_credentials_file(username='Titli',api_key='DYo3p5aNbKfX5cHX6CFo')


X = np.loadtxt('JaccardSimilarity_44PKABC_29_35_allkeys.txt')
print(type(X))

#44 Proteins
names=['1BX6_PKA','1CDK_PKA','1CMK_PKA','1J3H_PKA','1L3R_PKA','2C1A_PKA','2CPK_PKA','3TNQ_PKA','3NX8_PKA','4UJ1_PKA','4YXS_PKA','5O5M_PKA','5UZK_PKA','5VHB_PKA','4XW5_PKA','1ATP_PKA','1MRV_PKB','1GZK_PKB','1O6K_PKB','2JDO_PKB','2UW9_PKB','2X39_PKB','3D0E_PKB','3E87_PKB','6CCY_PKB','4GV1_PKB','3OCB_PKB','3QKK_PKB','4EKK_PKB','3CQU_PKB','3MV5_PKB','4EJN_PKB','5KCV_PKB','3O96_PKB','1GZO_PKB','2FK9_PKC','1A25_PKC','1GMI_PKC','1DSY_PKC','2NCE_PKC','4L1L_PKC','5W4S_PKC','3GPE_PKC','3RDJ_PKC']
print(len(names))

figure = ff.create_dendrogram(X, orientation='left', labels=names,linkagefun=lambda x: linkage(X, 'ward', metric='euclidean'))
py.plot(figure, filename='dendrogram_with_labels')
