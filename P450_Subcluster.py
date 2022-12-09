#%%
import os
cpass=os.getcwd()
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import seaborn as sns
import argparse

print("------------")
now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Start:"+str(sys.argv))
print("Ver.3.0.0_(2022_Dec.)")
print("by Hiroki Onoda NUSR")
print("------------")

# Argument & Help
parser = argparse.ArgumentParser(description="The program classify the protein sequences")
parser.add_argument('-n', default='5', help='nuber of amino acid sequence')
parser.add_argument('-r', default='5', help='nuber of recyclee')
parser.add_argument('-e', default='10', help='nuber of exponation')
parser.add_argument('-o', default='auto', help='input XXXXX if you want to name the output [XXXXX.csv]; ex) plddt (def.)')
parser.add_argument('-c', default='tab10_r', help='color of heat map ex) tab10_r (def.)')
parser.add_argument('-a', default='XX', help='Prefix of seqeunce name')
parser.add_argument('-b', default='00', help='Number of subcluster')
parser.add_argument('-m', default='2', help='minimum subcluster No.')
args = parser.parse_args()
# Argument & Help

#l----------------------------------------------l
#l           Valiable information               l
#l----------------------------------------------l

# Output file name
oname=str(args.o)
print("Output .png file name: ***"+str(oname)+"_***.png")
print("It can be changed by [-o "+oname+"]")
print("")


# Color
cname=str(args.c)
print("Output .png file color: "+str(cname))
print("It can be changed by [-c "+cname+"]")
print("")

#l-----------------------------------------------l
#l           Valiable information                l
#l-----------------------------------------------l

##########-Import-##########
##--Input sequence number
ssub=1
numb=int(str(args.n))
##--
##--MKDIR
os.mkdir("TSV")
os.mkdir("TSV/rank")
os.mkdir("TSV/class")
##
os.mkdir("CSV")
os.mkdir("CSV/rank")
os.mkdir("CSV/class")
os.mkdir("CSV/newclass")
##
os.mkdir("PNG")
os.mkdir("PNG/rearange")
os.mkdir("PNG/sc_rearange")
os.mkdir("PNG/chop")
os.mkdir("PNG/ncls")
os.mkdir("PNG/nrcls")
##--
##########-Import-##########



##########-Define-##########
##--Make matrix from tsv
def init():
    global cls1
    global newr
    ##--Initial setting
    cls0=np.zeros((numb,numb),dtype=float)
    inir=np.arange(numb,dtype=float)
    ##--
    
    print("Import:hits_sorted.tsv")
    with open("../hits_sorted.tsv", 'r') as go:
        for line in go:
            linen=line.split("	")
            lineo=linen[0].split("_")
            linez=linen[2].split("_")
            line0=int(lineo[1].strip())-1
            line2=int(linez[1].strip())-1
            line4=linen[4]
            if (line0)%ssub==w:
                if (line2)%ssub==w:
                    cls0[int((line0)/ssub),int((line2)/ssub)] = line4
                    pass
                pass
            pass
        pass
    cls1=cls0/100
    
    #--Decrease tje effect of slef effect
    for t in range(numb):
        cls1[t,t] = 0.15
        pass
    ##--
    
    ##--Save input matrix
    np.set_printoptions(threshold=np.inf)
    np.savetxt("TSV/class/"+oname+"_"+str(w)+"_i_class.tsv", cls1)
    ##--
    
    ##--rearange matrix by large cluster
    scor0=np.sum(cls1,axis=0)
    rank0=list(np.argsort(scor0))
    rankz=np.zeros(numb,dtype=float)
    symr0=np.zeros((numb,numb),dtype=float)
    for i in range(numb):
        symr0[i:i+1, :]=cls1[int(rank0[i]):int(rank0[i])+1, :]
        rankz[i]=int(rank0[i])+1
        pass
    cls1=np.zeros((numb,numb),dtype=float)
    for i in range(numb):
        cls1[:, i:i+1]=symr0[:, int(rank0[i]):int(rank0[i])+1]
        pass
    newr=rankz
    np.savetxt("CSV/rank/"+oname+"_"+str(w)+"_rank_0.csv", newr, delimiter=",")
    ##--
    np.savetxt("TSV/class/"+oname+"_"+str(w)+"_0_classt.tsv", cls1)
    
    ##--Save initial matrix
    plt.figure(figsize=(15, 12))
    ax = sns.heatmap(cls1, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
    plt.savefig("PNG/rearange_"+oname+"_"+str(w)+"_0.png")
    plt.savefig("PNG/rearange/rearange_"+oname+"_"+str(w)+"_0.png")
    plt.clf()
    plt.close()
    ##--
    
    ##--
    now = time.ctime()
    cnvtime = time.strptime(now)
    print(time.strftime("%Y/%m/%d", cnvtime))
    print("Finish_0:"+time.strftime("%H:%M:%S", cnvtime))
    pass
##--


##########-Define-##########
##--Make cruster
def clst():
    global cou
    global cls1
    global newr
    cou+=1
    if dup>1:
        clsd=np.power(cls1,dup)
        pass
    else:
        clsd=cls1
        pass
    
    ##--New next rank calicuration
    prob=clsd/np.full((numb,numb),np.sum(clsd,axis=1)).T
    mean=prob*np.full((numb,numb),np.array(range(numb))+1)
    scor=np.sum(mean,axis=1)
    rank=list(np.argsort(scor))
    
    ##--Matrix rearrangement
    newrz=np.zeros(numb,dtype=float)
    symr=np.zeros((numb,numb),dtype=float)
    for i in range(numb):
        symr[i:i+1, :]=cls1[int(rank[i]):int(rank[i])+1, :]
        newrz[i]=newr[int(rank[i])]
        pass
    cls1=np.zeros((numb,numb),dtype=float)
    for i in range(numb):
        cls1[:, i:i+1]=symr[:, int(rank[i]):int(rank[i])+1]
        pass
    newr=newrz
    np.savetxt("CSV/rank/"+oname+"_"+str(w)+"_rank_"+str(cou)+".csv", newr, delimiter=",")
    ##--
    
    ##--Save current matrix
    plt.figure(figsize=(15, 12))
    ax = sns.heatmap(cls1, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
    plt.savefig("PNG/rearange/rearange_"+oname+"_"+str(w)+"_"+str(cou)+".png")
    plt.clf()
    plt.close()
    ##

    now = time.ctime()
    cnvtime = time.strptime(now)
    print("Finish_"+str(cou)+":"+time.strftime("%H:%M:%S", cnvtime))
    ##--
    pass
##--


##########-Define-##########
##--Make class
def ssclst():
    global cls1
    global newr
    nop=5
    rann=np.zeros(numb)
    for d in range(numb):
        ccl=0
        maxp=0
        sann=0
        for k in range(nop):
            sann+=(k+1)*(k+1)
            cullen=np.mean(cls1[d:d+1,int(numb*(sann-(k+1)*(k+1))/55):int(numb*(sann)/55)])
            if cullen > 0.1:
                if maxp < cullen:
                    maxp=cullen
                    ccl=k
                    pass
                pass
            pass
        rann[d]=ccl
        pass
    ##--
    for k in range(nop):
        coo=0
        clnb=""
        axxa=np.count_nonzero(rann==k)
        cls2=np.zeros((axxa,numb), dtype=float)
        for d in range(numb):
            if rann[d]==k:
                cls2[coo:coo+1,:]=cls1[d:d+1,:]
                clnb+=str(int(newr[d])+1)+","
                coo+=1
                pass
            pass
        with open("CSV/class/c_"+str(nop)+"_"+str(w)+"_"+str(k)+"_class.csv", "w") as po:
            po.write(clnb)
            pass
        plt.figure(figsize=(25, 5))
        ax = sns.heatmap(cls2, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
        plt.savefig("PNG/chop/c_"+str(nop)+"_"+str(w)+"_"+str(k)+"_class.png")
        plt.clf()
        plt.close()
        pass
    pass
##--


##########-Run----##########
for w in range(ssub):
    init()
    cou=0
    for h in range(int(str(args.e))):
        dup=h+1
        for j in range(int(str(args.r))):
            clst()
            np.savetxt("TSV/class/"+oname+"_"+str(w)+"_"+str(cou)+"_classt.tsv", cls1)
            pass            
        plt.figure(figsize=(15, 12))
        ax = sns.heatmap(cls1, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
        plt.savefig("PNG/rearange_"+oname+"_"+str(w)+"_"+str(cou)+".png")
        plt.clf()
        plt.close()
        pass
    plt.figure(figsize=(15, 12))
    ax = sns.heatmap(cls1, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
    plt.savefig("rearange_"+oname+"_"+str(w)+"_"+str(cou)+".png")
    plt.savefig("../../../subclass_"+oname+"_"+str(args.a)+"_"+str(args.b)+"_"+str(cou)+".png")
    plt.clf()
    plt.close()
    ssclst()
    pass
##########-Run----##########
import subprocess
import csv
maxm=np.zeros(int(args.n),dtype=float)
numm=np.zeros(int(args.n),dtype=float)
for w in range(ssub):
    for ten in range (5):
        csv_file = open("CSV/class/c_5_"+str(w)+"_"+str(ten)+"_class.csv")
        csv_reader=csv.reader(csv_file)
        data_list=[]
        for row in csv_reader:
            if not row == []:
                roz=[]
                for n in row[:-1]:
                    roz.append(int(n)*ssub+1+w-ssub)
                    pass
                symj=np.zeros((len(roz),int(args.n)),dtype=float)
                symr=np.zeros((len(roz),len(roz)),dtype=float)
                count=0
                with open("../hits_sorted.tsv", 'r') as go:
                    for line in go:
                        linen=line.split("	")
                        lineo=linen[0].split("_")
                        line0=int(lineo[1].strip())
                        if line0 in roz:
                            linez=linen[2].split("_")
                            line2=int(linez[1].strip())
                            line4=linen[4]
                            symj[roz.index(line0),line2-1] = line4
                            if line2 in roz:
                                symr[roz.index(line0),roz.index(line2)] = line4
                            pass
                        pass
                    pass
                pass
            plt.figure(figsize=(15, 12))
            ax = sns.heatmap(symr, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=100)
            plt.savefig("PNG/ncls/newcls_5_"+str(w)+"_"+str(ten)+"_class.png")
            plt.clf()
            plt.close()
            
            ## Reclust
            avem=np.average(symj, axis=0)
            for jack in range(int(args.n)):
                if maxm[jack] < avem[jack]:
                    maxm[jack]=avem[jack]
                    numm[jack]=w*5+1+ten
                    pass
                pass
            pass
        pass
    pass
###
print(numm)
###


###
for k in range(5):
    coo=0
    clnb=""
    axxa=np.count_nonzero(numm==k+1)
    cls2=np.zeros((axxa,int(args.n)), dtype=float)
    for d in range(int(args.n)):
        if numm[d]==k+1:
            clnb+=str(d)+","
            coo+=1
            pass
        pass
    with open("CSV/newclass/rc_"+str(5*ssub)+"_"+str(k)+"_class.csv", "w") as po:
        po.write(clnb)
        pass
    pass
###
##
now = time.ctime()
cnvtime = time.strptime(now)
print("Finish_CLS_NB:"+time.strftime("%H:%M:%S", cnvtime))
##--

##########-Define-##########
def fig5():
    global symr
    cls1=symr/100
    #--Decrease tje effect of slef effect
    for t in range(int(len(roz))):
        cls1[t,t] = 0.15
        pass
    ##--rearange matrix by large cluster
    scor0=np.sum(cls1,axis=0)
    rank0=list(np.argsort(scor0))
    rankz=np.zeros(int(len(roz)),dtype=float)
    symr0=np.zeros((int(len(roz)),int(len(roz))),dtype=float)
    for i in range(int(len(roz))):
        symr0[i:i+1, :]=cls1[int(rank0[i]):int(rank0[i])+1, :]
        rankz[i]=int(rank0[i])+1
        pass
    cls1=np.zeros((int(len(roz)),int(len(roz))),dtype=float)
    for i in range(int(len(roz))):
        cls1[:, i:i+1]=symr0[:, int(rank0[i]):int(rank0[i])+1]
        pass
    ##--Save initial matrix
    plt.figure(figsize=(15, 12))
    ax = sns.heatmap(cls1, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=1)
    plt.savefig("PNG/sc_rearange/rearange_"+oname+"_"+str(args.a)+"_"+str(args.b)+"_5_"+str(ten)+".png")
    plt.clf()
    plt.close()
    ##--
    
    ##--
    now = time.ctime()
    cnvtime = time.strptime(now)
    print(time.strftime("%Y/%m/%d", cnvtime))
    print("Finish_0:"+time.strftime("%H:%M:%S", cnvtime))
    pass
##--


maxm=np.zeros(int(args.n),dtype=float)
numm=np.zeros(int(args.n),dtype=float)
for w in range(ssub):
    for ten in range (5):
        csv_file = open("CSV/newclass/rc_"+str(5*ssub)+"_"+str(ten)+"_class.csv")
        csv_reader=csv.reader(csv_file)
        data_list=[]
        for row in csv_reader:
            if not row == []:
                roz=[]
                for n in row[:-1]:
                    roz.append(int(n))
                    pass
                symj=np.zeros((len(roz),int(args.n)),dtype=float)
                symr=np.zeros((len(roz),len(roz)),dtype=float)
                count=0
                with open("../hits_sorted.tsv", 'r') as go:
                    for line in go:
                        linen=line.split("	")
                        lineo=linen[0].split("_")
                        line0=int(lineo[1].strip())
                        if line0 in roz:
                            linez=linen[2].split("_")
                            line2=int(linez[1].strip())
                            line4=linen[4]
                            symj[roz.index(line0),line2-1] = line4
                            if line2 in roz:
                                symr[roz.index(line0),roz.index(line2)] = line4
                            pass
                        pass
                    pass
                pass
            plt.figure(figsize=(15, 12))
            ax = sns.heatmap(symr, linewidth=0, cmap=cname, square=True, xticklabels=100, yticklabels=100, vmin=0, vmax=100)
            plt.savefig("PNG/nrcls/newrcls_"+str(5*ssub)+"_"+str(ten+1)+"_reclass.png")
            plt.clf()
            plt.close()
            #
            fig5()
            #
            ## Reclust
            avem=np.average(symj, axis=0)
            for jack in range(int(args.n)):
                if maxm[jack] < avem[jack]:
                    maxm[jack]=avem[jack]
                    numm[jack]=1+ten
                    pass
                pass
            pass
        pass
    pass
###
print(numm)
###

cocon=np.zeros(5,dtype=float)
for t in range(6):
    os.mkdir(str(t))
for pety in range(int(args.n)):
    cocon[int(numm[pety])-1]+=1
    papa=int(cocon[int(numm[pety])-1])
    os.symlink("../../../../"+str(args.b)+"/"+str(args.a)+str(args.b)+"_"+"0"*(7-len(str(pety+1)))+str(pety+1)+".faa", str(int(numm[pety]))+"/"+str(args.a)+str(args.b)+"-"+str(int(numm[pety]))+"_"+"0"*(7-len(str(papa)))+str(papa)+".faa")
    pass
