import urllib.request    
from urllib.request import urlopen  
import multiprocessing
from joblib import Parallel, delayed

# 44 PKABC
fileList = ['1BX6','1CDK','1CMK','1J3H','1L3R','2C1A','2CPK','3TNQ','3NX8','4UJ1','4YXS','5O5M','5UZK','5VHB','4XW5','1ATP','1MRV','1GZK','1O6K','2JDO','2UW9','2X39','3D0E','3E87','6CCY','4GV1','3OCB','3QKK','4EKK','3CQU','3MV5','4EJN','5KCV','3O96','1GZO','2FK9','1A25','1GMI','1DSY','2NCE','4L1L','5W4S','3GPE','3RDJ']
print (len(fileList))
def parallelcode(fname):
    print (fname)
    
    url = urllib.request.urlretrieve('http://files.rcsb.org/download/'+fname+'.pdb', fname+'.pdb')
    try:
        urlopen(url)
    except:
        print (fname+'not found')
        pass

def getInputs():
    #num_cores = multiprocessing.cpu_count()
    #Parallel(n_jobs=num_cores, verbose=50)(delayed(parallelcode)(fname)for fname in fileList)

    for fname in fileList:
        print (fname)
    
    url = urllib.request.urlretrieve('http://files.rcsb.org/download/'+fname+'.pdb', fname+'.pdb')
    try:
        urlopen(url)
    except:
        print (fname+'not found')
        pass
getInputs()

print ("code end.")


