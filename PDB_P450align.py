#%%
import os
cpass=os.getcwd()
import jax
import math
import sys
import pickle

import matplotlib.pyplot as plt
import numpy as np

#---------------import-------------------
import argparse
import time
#---------------arg----------------------
parser = argparse.ArgumentParser(description="The runner can accept arguments")
parser.add_argument('-i', default='PDB/test.pdb', help='relative fasta path of imput sequence ; ex) fasta/test.fasta (def.)')
parser.add_argument('-jo', "--jobname", default=0,                     help='input XXXXX if you want to name the output directory = predition_XXXXX_<hash>; ex) XXXX of XXXX.fasta (def.)')
parser.add_argument('-stp',"--no_save_tmp_pdb", action='store_false',  help='dont save tmp.pdb')
parser.add_argument('-ali',"--save_aligned_pdb", action='store_true',  help='save aligned.pdb at first chain (A chain)')
parser.add_argument('-alc', default="1",                               help='Change chain number from first chain to other chain defalt: 1')
parser.add_argument('-al1', default=0,                                 help='Center the mid point of al1 and al2, the point (al1) shift on x-axis')
parser.add_argument('-al2', default=0,                                 help='Center the mid point of al1 and al2, the point (al2) shift on x-asis and opozit site of al1')
parser.add_argument('-al3', default=0,                                 help='the point (al3) shift on xy-plane')
args = parser.parse_args()


##--MKDIR

if not os.path.isdir("o"):
    os.mkdir("o")
    os.mkdir("o/L")
    pass
##

####################
def CheckP():
    global w4
    Aa9=""
    Aa8=""
    Aa7=""
    Aa6=""
    Aa5=""
    Aa4=""
    Aa3=""
    Aa2=""
    Aa1=""
    Aa0=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line [0:4]=="ATOM":
                if line[12:16].strip()=="CA":
                    #print(line[17:20].strip())
                    Aa9=str(Aa8)
                    Aa8=str(Aa7)
                    Aa7=str(Aa6)
                    Aa6=str(Aa5)
                    Aa5=str(Aa4)
                    Aa4=str(Aa3)
                    Aa3=str(Aa2)
                    Aa2=str(Aa1)
                    Aa1=str(Aa0)
                    Aa0=line[17:20].strip()
                    if Aa2=="CYS":
                        if Aa0=="GLY" or Aa0=="ALA":
                            if Aa9=="PHE" or Aa9=="TRP" or Aa9=="TYR":
                                w4=int(float(line[22:26].strip()))
                                o4=w4-9
                                t4=w4-2
                                print(Aa9+","+Aa2+","+Aa0+":"+line[22:26].strip())
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    pass

def ProxiP():
    lines=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line [0:4]=="ATOM":
                if float(line[46:54].strip())<4:
                    lines+=line
                    pass
                pass
            elif line [0:4]=="ATOM" or line [0:6]=="HETATM":
                if line[17:20].strip()=="HEM":
                    lines+=line
                    pass
                pass
            pass
        pass
    with open(cpass+"/p/"+args.i[:-4]+"_ploxi.pdb", 'w') as f:
        f.write(lines)
        pass
    pass

def HemeaH():
    lines=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line [0:4]=="ATOM" or line [0:6]=="HETATM":
                if line[17:20].strip()=="HEM":
                    lines+=line
                    pass
                pass
            pass
        pass
    with open(cpass+"/h/"+args.i[:-4]+"_heme.pdb", 'w') as f:
        f.write(lines)
        pass
    pass

####################
def ArgHis():
    lines=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line [0:4]=="ATOM":
                if float(line[46:54].strip())<8 and float(line[30:38].strip())<-4 and float(line[38:46].strip())<0 and float(line[30:38].strip())>-15 and float(line[38:46].strip())>-15 and float(line[46:54].strip())>-15:
                    lines+=line
                    pass
                elif int(float(line[22:26].strip()))>w4-11 and int(float(line[22:26].strip()))<w4+2:
                    lines+=line
                    pass
                pass
            elif line [0:4]=="ATOM" or line [0:6]=="HETATM":
                if line[17:20].strip()=="HEM":
                    lines+=line
                    pass
                pass
            pass
        pass
    with open(cpass+"/RH/"+args.i[:-4]+"_ploxi.pdb", 'w') as f:
        f.write(lines)
        pass
    pass

def CheckC():
    Aa9=""
    Aa8=""
    Aa7=""
    Aa6=""
    Aa5=""
    Aa4=""
    Aa3=""
    Aa2=""
    Aa1=""
    Aa0=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line [0:4]=="ATOM":
                if line[12:16].strip()=="CA":
                    #print(line[17:20].strip())
                    Aa9=str(Aa8)
                    Aa8=str(Aa7)
                    Aa7=str(Aa6)
                    Aa6=str(Aa5)
                    Aa5=str(Aa4)
                    Aa4=str(Aa3)
                    Aa3=str(Aa2)
                    Aa2=str(Aa1)
                    Aa1=str(Aa0)
                    Aa0=line[17:20].strip()                    
                    if Aa2=="CYS":
                        if Aa0=="GLY" or Aa0=="ALA":
                            if Aa9=="PHE" or Aa9=="TRP" or Aa9=="TYR":
                                w4=int(float(line[22:26].strip()))
                                o4=w4-9
                                t4=w4-2
                                print(Aa9+","+Aa2+","+Aa0+":"+line[22:26].strip())
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    print(str(o4)+","+str(t4)+","+str(w4))

    novch={"1":"A", "2":"B", "3":"C", "4":"D", "5":"E", "6":"F", "7":"G", "8":"H", "9":"I", "10":"J", "11":"K", "12":"L", "13":"M", "14":"N", "15":"O", "16":"P", "17":"Q", "18":"R", "19":"S", "20":"T", "21":"U", 22:"V", "23":"W", "24":"X", "25":"Y", "26":"Z"}
    cn=novch[args.alc]

    print("align residues were set to ("+ str(o4)+","+ str(w4)+","+ str(t4)+")")

    xyz1=np.zeros((4,3),dtype=float)
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line[21:22]==cn:
                if line[0:4]=="ATOM":
                    if o4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line[:-2])
                            xyz1[0,0] = float(line[30:38].strip())
                            xyz1[0,1] = float(line[38:46].strip())
                            xyz1[0,2] = float(line[46:54].strip())

                    elif w4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line[:-2])
                            xyz1[1,0] = float(line[30:38].strip())
                            xyz1[1,1] = float(line[38:46].strip())
                            xyz1[1,2] = float(line[46:54].strip())

                    elif t4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line[:-2])
                            xyz1[2,0] = float(line[30:38].strip())
                            xyz1[2,1] = float(line[38:46].strip())
                            xyz1[2,2] = float(line[46:54].strip())

    print("------------")
    xyz1[3,0] = xyz1[2,0]
    xyz1[3,1] = xyz1[2,1]
    xyz1[3,2] = xyz1[2,2]
    xyzm=np.zeros((1,3),dtype=float)
    xyzm[0,0] = xyz1[3,0]
    xyzm[0,1] = xyz1[3,1]
    xyzm[0,2] = xyz1[3,2]
    print("Initial xyz")
    print(xyz1)
    #print(xyzm)
    xyz2=xyz1-xyzm
    #print("Center xyz")
    #print(xyz2)
    #print("------------")
    # rotation z fit oy and ty to zero
    if xyz2[0,1]>0:
        pmy=-1

    else:
        pmy=1

    distaxy=math.sqrt(math.pow(xyz2[0,0], 2)+math.pow(xyz2[0,1], 2))
    radax=pmy*math.acos(xyz2[0,0]/distaxy)
    rot1=np.zeros((3,3),dtype=float)
    rot1[0,0]=math.cos(radax)
    rot1[1,0]=-math.sin(radax)
    rot1[0,1]=math.sin(radax)
    rot1[1,1]=math.cos(radax)
    rot1[2,2]=1
    xyz3=np.dot(xyz2,rot1)
    print("rotation z to fit 1y and 3y to zero")
    #print(xyz3)
    #print("------------")
    # rotation y fit ox to zero
    if xyz3[0,2]>0:
        pmz=-1

    else:
        pmz=1

    dista1=math.sqrt(math.pow(xyz3[0,0], 2)+math.pow(xyz3[0,1], 2)+math.pow(xyz3[0,2], 2))
    radaxy=pmz*math.acos(xyz3[0,0]/dista1)
    rot2=np.zeros((3,3),dtype=float)
    rot2[0,0]=math.cos(radaxy)
    rot2[2,0]=-math.sin(radaxy)
    rot2[0,2]=math.sin(radaxy)
    rot2[2,2]=math.cos(radaxy)
    rot2[1,1]=1
    xyz4=np.dot(xyz3,rot2)
    print("rotation y to fit 1z and 3z to zero")
    #print(xyz4)
    #print("------------")
    # rotation x fit wz to zero
    if xyz4[1,2]>0:
        pmyz=-1

    else:
        pmyz=1

    distbbyz=math.sqrt(math.pow(xyz4[1,1], 2)+math.pow(xyz4[1,2], 2))
    radbbyz=pmyz*math.acos(xyz4[1,1]/distbbyz)
    rot3=np.zeros((3,3),dtype=float)
    rot3[1,1]=math.cos(radbbyz)
    rot3[2,1]=-math.sin(radbbyz)
    rot3[1,2]=math.sin(radbbyz)
    rot3[2,2]=math.cos(radbbyz)
    rot3[0,0]=1
    xyz5=np.dot(xyz4,rot3)
    print("rotation x to fit 2z to zero")
    #print(xyz5)
    #print("------------")
    # Check all shift & roation 
    xyz6=np.dot(np.dot(np.dot(xyz1-xyzm,rot1),rot2),rot3)
    print("Check all shift & roation")
    print(xyz6)
    print("------------")
    xyz=np.zeros((1,3),dtype=float)
    liii=""
    with open(cpass+"/"+args.i, 'r') as go:
        for line in go:
            if line[0:4]=="ATOM":
                xyz[::1]=0
                xyz[0,0] = float(line[30:38].strip())
                xyz[0,1] = float(line[38:46].strip())
                xyz[0,2] = float(line[46:54].strip())
                xyzn=np.dot(np.dot(np.dot(xyz-xyzm,rot1),rot2),rot3)
                qxxl=" "*(8-len(str(round(xyzn[0,0],4))))+str(round(xyzn[0,0],4))
                qyyl=" "*(8-len(str(round(xyzn[0,1],4))))+str(round(xyzn[0,1],4))
                qzzl=" "*(8-len(str(round(xyzn[0,2],4))))+str(round(xyzn[0,2],4))
                liii+=line[0:30]+qxxl+qyyl+qzzl+line[54:]
            elif line[0:6]=="HETATM":
                xyz[::1]=0
                xyz[0,0] = float(line[30:38].strip())
                xyz[0,1] = float(line[38:46].strip())
                xyz[0,2] = float(line[46:54].strip())
                xyzn=np.dot(np.dot(np.dot(xyz-xyzm,rot1),rot2),rot3)
                qxxl=" "*(8-len(str(round(xyzn[0,0],4))))+str(round(xyzn[0,0],4))
                qyyl=" "*(8-len(str(round(xyzn[0,1],4))))+str(round(xyzn[0,1],4))
                qzzl=" "*(8-len(str(round(xyzn[0,2],4))))+str(round(xyzn[0,2],4))
                liii+=line[0:30]+qxxl+qyyl+qzzl+line[54:]
            else:
                liii+=line

    with open(cpass+"/o/"+args.i[:-4]+"_align.pdb", 'w') as f:
        f.write(liii)
        pass
    pass
####################

CheckP()
#ArgHis()
#HemeaH()
#ProxiP()
CheckC()
#AlignA()
