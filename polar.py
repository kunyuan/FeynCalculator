import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob, os, re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size=12

rs=1
Beta=10
# kF=(9.0*np.pi/4.0)**(1.0/3.0)/rs #3D
kF=np.sqrt(2.0)/rs #2D
# Bubble=0.08871  # 3D, Beta=0.5, rs=1
# Bubble=0.09720  #3D, Beta=10, rs=1
# Bubble=0.11635  #2D, Beta=0.5, rs=1
Bubble=0.15916  #2D, Beta=10, rs=1

ScanOrder=[1,2]
# ScanOrder=[2]
Index={}
Index[1]=[1,]
Index[2]=[1,]
Index[3]=[1,2,3,4,5]
DataAll={}
Data={}
Normalization=1


folder="./Beta"+str(Beta)+"_rs1.0/Data/"
# os.chdir(folder)
files=os.listdir(folder)
for order in ScanOrder:
    Num=0
    data0=None
    for f in files:
        if re.match("Diag"+str(order)+"_[0-9]+.dat", f):
            print f
            Num+=1
            d=np.loadtxt(folder+f)
            if data0 is None:
                data0=d
            else:
                data0[:,1:]+=d[:,1:]

    print "Found {0} files.".format(Num)
    data0[:,1:]/=Num

    if order==1:
        Normalization=data0[0,1]/Bubble
    data0[:,1]/=Normalization
    data0[:,1]*=(-1)**(order-1)

    # print data0

    DataAll[order]=np.array(data0)

    for i in Index[order]:
        Data[i]=[]
        Num=0
        data=None
        # for f in glob.glob("Diag"+str(order)+"_*_"+str(i)+".dat"):
        for f in files:
            if re.match("Diag"+str(order)+"_[0-9]+_"+str(i)+".dat", f):
                print f
                Num+=1
                d=np.loadtxt(folder+f)
                if data is None:
                    data=d
                else:
                    data[:,1:]+=d[:,1:]
        print "Found {0} files.".format(Num)
        data[:,1:]/=Num
        data[:,1]/=Normalization
        data[:,1]*=(-1)**(order-1)
        print data
        Data[i].append(np.array(data))

# Data[2][0][:,1]/=8.0*np.pi/(1.0+Data[2][0][:,0]**2)
# Data[2][0][:,1]=-(Data[2][0][:,1])**0.5
        

def ErrorPlot(p, data, color, marker, label=None, size=4):
    data[:,0]/=kF
    # data[:,1]-=data[0,1]
    p.plot(data[:,0],data[:,1],marker=marker,c=color, label=label,lw=1, markeredgecolor="None", linestyle="--", markersize=size)
    # p.errorbar(data[:,0],data[:,1], yerr=data[:,2], c=color, ecolor=color, capsize=0, linestyle="None")
    # p.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], alpha=0.5, facecolor=color, edgecolor=color)

w=1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)

ErrorPlot(ax, DataAll[1], 'k', 's', "Order 1")
ErrorPlot(ax, DataAll[2], 'r', 'o', "Order 2")
# ErrorPlot(ax, DataAll[3], 'g', 'o', "Order 3")
# ErrorPlot(ax, Data[2][1], 'b', 'o', "Order 2-Order 2")
# ErrorPlot(ax, Data[0], 'r', 's', "Diag 1")
# ErrorPlot(ax, Data[1], 'b', 's', "Diag 2")
# ErrorPlot(ax, Data[2], 'g', 's', "Diag 3")
# ErrorPlot(ax, Data[3], 'm', 's', "Diag 4")
# ErrorPlot(ax, Data[4], 'c', '*', "Diag 5")

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")


x=np.arange(0,0.2,0.001)
y=0.5*x**w

# ax.plot(x,y,'k-', lw=2)

ax.set_xlim([0.0, 3.0])
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.01, 0.0])
ax.set_xlabel("$q/k_F$", size=size)
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("polarization.pdf")
plt.show()

