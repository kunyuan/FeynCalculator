import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob, os
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size=12

rs=1
kF=(9.0*np.pi/4.0)**(1.0/3.0)/rs #3D
# kF=np.sqrt(2.0)/rs #2D

# Index=[1,2,3,4]
# Index=[1,2,3,4,5]
Index=[1,]
Data=[]
Normalization=0

CurrOrder=2

folder="./Beta10.0_rs1.0/Data/"
os.chdir(folder)
Num=0
data0=None
for f in glob.glob("Diag"+str(CurrOrder)+"_"+"?.dat"):
    print f
    Num+=1
    d=np.loadtxt(f)
    if data0 is None:
        data0=d
    else:
        data0[:,1:]+=d[:,1:]

data0[:,1:]/=Num
# Normalization=abs(sum(data0[:,1]))
Normalization=data0[0,1]
data0[:,1]/=Normalization

for i in Index:
    Num=0
    data=None
    for f in glob.glob("Diag"+str(CurrOrder)+"_?_"+str(i)+".dat"):
        print i, f
        Num+=1
        d=np.loadtxt(f)
        if data is None:
            data=d
        else:
            data[:,1:]+=d[:,1:]
    data[:,1:]/=Num
    data[:,1]/=Normalization
    Data.append(data)

# Data[3][:,1]*=2

# folder="../Beta10.0_rs1.0_1"
# os.chdir(folder)
# Num=0
# data2=None
# for f in glob.glob("*.dat"):
    # Num+=1
    # d=np.loadtxt(f)
    # if data2 is None:
        # data2=d
    # else:
        # data2[:,1:]+=d[:,1:]

# data2[:,1:]/=Num
# data2[:,1]/=abs(sum(data2[:,1]))

# folder="../Beta20.0_rs1.0"
# os.chdir(folder)
# Num=0
# data3=None
# for f in glob.glob("*.dat"):
    # Num+=1
    # d=np.loadtxt(f)
    # if data3 is None:
        # data3=d
    # else:
        # data3[:,1:]+=d[:,1:]

# data3[:,1:]/=Num
# data3[:,1]/=abs(sum(data3[:,1]))

# exit(0)
# L14=np.loadtxt("./L128_1.5056/Int1D_2.dat")

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

# ErrorPlot(ax, L1, 'g', '>', "$L=16, t=0.01$")
# ErrorPlot(ax, L2, 'm', '<', "$L=32, t=0.01$")

# ErrorPlot(ax, L0, 'y', '*', "$BK, \gamma=0.3332052$")
# ErrorPlot(ax, L1, 'b', '*', "$t=0.5$")
# ErrorPlot(ax, L2, 'r', 's', "$t=0.4$")
# ErrorPlot(ax, L2, 'r', 's', "$t=0.4$")
# ErrorPlot(ax, data3, 'g', 's', "Beta=10_1")
# ErrorPlot(ax, data2, 'k', 's', "Beta=10_0")
# ErrorPlot(ax, data3, 'k', 's', "$T=E_F/40$")
# ErrorPlot(ax, data2, 'b', 's', "Beta=10_1")

ErrorPlot(ax, data0, 'k', 's', "Total")
ErrorPlot(ax, Data[0], 'r', 's', "Diag 1")
# ErrorPlot(ax, Data[1], 'b', 's', "Diag 2")
# ErrorPlot(ax, Data[2], 'g', 's', "Diag 3")
# ErrorPlot(ax, Data[3], 'm', 's', "Diag 4")
# ErrorPlot(ax, Data[4], 'c', '*', "Diag 5")

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")
# ErrorPlot(ax, L4, 'm', '<', "$t=1.0$")
# ErrorPlot(ax, L5, 'y', 'v', "$BK, \gamma=2.0$")
# ErrorPlot(ax, L6, 'b', 'o', "$t=10.0$")
# ErrorPlot(ax, L7, 'c', 'o', "$trap center$")
# ErrorPlot(ax, L8, 'k', 'v', "S=1 BK")
# ErrorPlot(ax, DI2, 'c', 'o', "S=1/2 2BK, $d=2$")
# ErrorPlot(ax, DI4, 'b', '<', "S=1/2 2BK, $d=4$")
# ErrorPlot(ax, DI8, 'r', 'v', "S=1/2 2BK, $d=8$")
# ErrorPlot(ax, DI16, 'm', 's', "S=1/2 2BK, $d=16$")


x=np.arange(0,0.2,0.001)
# y=0.5-0.5*x**w
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
ax.set_ylabel("normalized $P(\omega=0, q)-P(\omega=0, q=0)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

os.chdir("../")
# plt.savefig("polarization.pdf")
plt.show()

