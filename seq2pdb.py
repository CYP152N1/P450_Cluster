#%%
from urllib import request
import os
import argparse
import time
import sys


print("------------")
now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Start:"+str(sys.argv))
print("Ver.1.0_(2022_Dec.)")
print("by Hiroki Onoda NUSR")
print("------------")


parser = argparse.ArgumentParser(description="The program is AF_pdb scraiper from fasta file")
parser.add_argument('-i', default='XXNN', help='subfig')
args = parser.parse_args()

#l----------------------------------------------l
#l           Valiable information               l
#l----------------------------------------------l

# Output file name
fname=str(args.i)

#l-----------------------------------------------l
#l           Valiable information                l
#l-----------------------------------------------l


##--MKDIR
os.mkdir(fname+"_AF")
os.mkdir(fname+"_AF/V1")
os.mkdir(fname+"_AF/V2")
os.mkdir(fname+"_AF/V3")
os.mkdir(fname+"_AF/V4")
os.mkdir(fname+"_AF/L")
##

files=os.listdir("./"+fname)
files=sorted(files)
logg=""
for f in files:
    logg=""
    print(f)
    cou=0
    with open (fname+"/"+f) as fa:
        for line in fa:
            if line[0:1]==">":
                try:
                    source = request.urlopen("https://alphafold.ebi.ac.uk/files/AF-"+line[1:].rstrip()+"-F1-model_v1.pdb") 
                    data=source.read()
                    text=data.decode()
                    with open (fname+"_AF/V1/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open (fname+"_AF/L/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open(fname+"_AF/V1.txt", 'w') as j:
                        j.write(line)
                        j.close()
                        pass
                    logg=f[:-4]+","+line[1:].rstrip()+",1"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                    cou+=1
                except:
                    logg=f[:-4]+","+line[1:].rstrip()+",0"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                pass
            pass
        pass
    pass
    logg=""
    with open (fname+"/"+f) as fa:
        for line in fa:
            if line[0:1]==">":
                try:
                    source = request.urlopen("https://alphafold.ebi.ac.uk/files/AF-"+line[1:].rstrip()+"-F1-model_v2.pdb") 
                    data=source.read()
                    text=data.decode()
                    with open (fname+"_AF/V2/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open (fname+"_AF/L/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open(fname+"_AF/V2.txt", 'a') as j:
                        j.write(line)
                        j.close()
                        pass
                    logg=",1"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                    cou+=1
                except:
                    logg=",0"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                pass
            pass
        pass
    logg=""
    with open (fname+"/"+f) as fa:
        for line in fa:
            if line[0:1]==">":
                try:
                    source = request.urlopen("https://alphafold.ebi.ac.uk/files/AF-"+line[1:].rstrip()+"-F1-model_v3.pdb") 
                    data=source.read()
                    text=data.decode()
                    with open (fname+"_AF/V3/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open (fname+"_AF/L/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open(fname+"_AF/V3.txt", 'a') as j:
                        j.write(line)
                        j.close()
                        pass
                    logg=",1"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                    cou+=1
                except:
                    logg=",0"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                pass
            pass
        pass
    logg=""
    with open (fname+"/"+f) as fa:
        for line in fa:
            if line[0:1]==">":
                try:
                    source = request.urlopen("https://alphafold.ebi.ac.uk/files/AF-"+line[1:].rstrip()+"-F1-model_v4.pdb") 
                    data=source.read()
                    text=data.decode()
                    with open (fname+"_AF/V4/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open (fname+"_AF/L/"+f[:-4]+".pdb", mode="w") as go:
                        go.write(text)
                        go.close()
                        pass
                    with open(fname+"_AF/V4.txt", 'a') as j:
                        j.write(line)
                        j.close()
                        pass
                    logg=",1"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                    cou+=1
                except:
                    logg=",0"
                    with open(fname+"_AF/AF.csv", 'a') as j:
                        j.write(logg)
                        j.close()
                        pass
                    pass
                pass
            pass
        pass
    logg=","+str(cou)+"\n"
    with open(fname+"_AF/AF.csv", 'a') as j:
        j.write(logg)
        j.close()
        pass
    pass


