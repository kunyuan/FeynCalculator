import numpy as np

rs=1.0
Lambda=1.0
Beta=25

##############   3D    ##################################
kF=(9.0*np.pi/4.0)**(1.0/3.0)/rs #3D
###### Bare Green's function    #########################
# Bubble=0.08871  # 3D, Beta=0.5, rs=1
# Bubble=0.0971916  #3D, Beta=10, rs=1
# Bubble=0.0971613  #3D, T=0.04Ef, rs=1
# Bubble= 0.097226 # 3D, zero temperature, rs=1
###### Fock dressed Green's function ###################
Bubble, Density, dMu2=0.088883, 0.2387, -0.2699 #3D, Beta=0.1, rs=1, Lambda=1.0

##############   2D    ##################################
###### Bare Green's function    #########################
# kF=np.sqrt(2.0)/rs #2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble=0.15916  #2D, Beta=10, rs=1


folder="./Beta{0}_rs{1}_lambda{2}_freq/".format(Beta, rs, Lambda) 
RawData=np.loadtxt(folder+"output.dat")
Order={}
OrderList=[1,2,3]
for order in OrderList:
    Order[order]=RawData[RawData[:,0]==order]

Norm=Order[1][:,1]/Bubble
for order in OrderList:
    Order[order][:,1]/=Norm

for order in OrderList:
    print "Order {0}: {1}+-{2}".format(order, Order[order][:,1].mean(), Order[order][:,1].std()/np.sqrt(len(Order[order][:,1])))
    print Order[order]
