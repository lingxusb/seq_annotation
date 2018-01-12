import sys
import numpy as np
import scipy.io
import matplotlib.pylab as plt

#find TSS before start codon
def fore_TSS(expression,BHI_R1f,pos_TIS,i):
    r= list(pos_TIS[i,])
    for j in range(150):
            if expression[i] > 2 and BHI_R1f[pos_TIS[i,0]-j,1] == 0 and \
            (np.mean(BHI_R1f[pos_TIS[i,0]-j-20:pos_TIS[i,0]-j,1])<expression[i]/10 or\
             np.mean(BHI_R1f[pos_TIS[i,0]-j-8:pos_TIS[i,0]-j,1])==0):
                r[0] = pos_TIS[i,0]-j
                break
    return r

#find TSS after start codon
def after_TSS(expression,BHI_R1f,pos_TIS,i):
    r = 0
    if expression[i]>2 and pos_TIS[i,0] == pos_TIS[i,0] and BHI_R1f[pos_TIS[i,0],1] == 0:
            for k in range(50):
                if BHI_R1f[pos_TIS[i,0]+k,1] != 0:
                    r = k-1
                    break
    return r

#find terminator after end codon
def after_TAA(expression,BHI_R1f,pos_TIS,i):
    r = 0
    for j in range(150):
            if expression[i] > 2 and BHI_R1f[pos_TIS[i,1]+j,1] == 0 and \
            np.mean(BHI_R1f[pos_TIS[i,1]+j:pos_TIS[i,1]+j+8,1]) == 0:
                r = j-1
                break
    return r
    
#refine transcript start and end sites
def refine_TSS(BHI_R1f, pos_TIS):
    BHI_R1f = BHI_R1f
    expression = np.zeros((pos_TIS.shape[0],1))
    pos_TS1 = np.zeros(pos_TIS.shape)
    for i in range(pos_TIS.shape[0]):
        expression[i] = np.mean(BHI_R1f[pos_TIS[i,0]:pos_TIS[i,1],1])
        
        #refine transcript start site
        pos_TS1[i] = np.array(fore_TSS(expression,BHI_R1f,pos_TIS,i))
        pos_TS1[i,0] += after_TSS(expression,BHI_R1f,pos_TIS,i)
        
        #refine transcript end site
        pos_TS1[i,1] += after_TAA(expression,BHI_R1f,pos_TIS,i)        

    return expression,pos_TS1

#search for orphan promoters, defined at least 150bp away from gene coding regions.
def pro_ophan(counts):
    posl = []
    strenl = []
    pos = [0,0]
    stren = 0
    i=50
    while i < len(counts)-50:
        i += 1
        if counts[i]!=0 and (np.mean(counts[i:i+5])-counts[i])>=3*np.mean(counts[i-50:i])+1:
            pos[0] = i-1
            pos[1] = len(counts)
            stren = np.mean(counts[i:i+50])
            posl.append(list(pos))
            strenl.append(stren)
            i = i + 500
    return posl, strenl
    
def forphan(BHI_R1f,pos_TS1):
    exp = []
    sites  = []
    for i in range(len(pos_TS1)-1):
        pos, stren = pro_ophan(BHI_R1f[pos_TS1[i,1]:pos_TS1[i+1,0],1])
        if len(pos) > 0:
            exp = exp + stren
            sites= sites + (pos+pos_TS1[i,1]).tolist()
    return np.resize(np.array(exp),(len(exp),1)), np.array(sites)

#search for internal promoters
def pro_internal(counts):
    posl = []
    strenl = []
    pos = [0,0]
    stren = 0
    i = 50
    while i < len(counts)-50:
        i += 1
        if counts[i]!=0 and (np.mean(counts[i:i+5])-counts[i])>=3*np.mean(counts[i-50:i])+1:
            pos[0] = i-1
            pos[1] = len(counts)
            stren = np.mean(counts[i:i+50])
            posl.append(list(pos))
            strenl.append(stren)
            i = i + 500
    return posl, strenl
    
def finternal(BHI_R1f,pos_TS1):
    exp = []
    sites  = []
    for i in range(len(pos_TS1)):
        pos, stren = pro_internal(BHI_R1f[pos_TS1[i,0]:pos_TS1[i,1],1])
        if len(pos) > 0:
            exp= exp + stren
            sites= sites+(pos+pos_TS1[i,0]).tolist()
    return np.resize(np.array(exp),(len(exp),1)), np.array(sites)

#load RNA seq data(position based) for forward strand of chromsome 1
BHI_R1f = np.genfromtxt('D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\BHI\\RNA-seq\\wig\\CP016345.1_f.wig')
#load position gene coding region for forward strand of chr1 (start and end position)
pos_TIS1 = np.loadtxt('Data\\pos_chr1.txt')
pos_TIS1 = pos_TIS1.astype(int)
print pos_TIS1[0:10,]

#find transcription start/end sites near gene annotation
exp1,sites1 = refine_TSS(BHI_R1f,pos_TIS1)
#find transcription start/end sites between gene coding regions
exp2,sites2 = forphan(BHI_R1f,pos_TIS1)
#find TSS in gene coding regions
exp3,sites3 = finternal(BHI_R1f,pos_TIS1)
