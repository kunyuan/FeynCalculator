#!/usr/bin/python
import os

# FileName="pyrochlo.dat"
os.system("rm output.dat")
os.system("cat ./Beta25_rs1.0_lambda1.0_freq/output.dat>> output.dat")
FileName="output.dat"
Head="#para"
Parameters=[1,2]
# DataRank={1:"Epot", 4:"E2", 3: "Energy", 17:"M_Eq",26:"Binder", 28:"Binder_Eq", 8:"wind"}
DataRank={1:"Order1", 2:"Order2", 3:"Order3"}

############   Read data from file and cut them into blocks   ####################

f=open(FileName,"r")
DataRankList=DataRank.keys()
print DataRankList
print DataRank.values()
#Divide the data file into several data blocks
Blocks=[elem.strip() for elem in f.read().split(Head)]
Blocks=filter(lambda x:x!="",Blocks)
#print Blocks[0]
f.close()

IndexMap=[]
alllines=Blocks[0].split("\n")
for line in alllines[1:]:
    try:
        if(int(line.split()[0]) in DataRankList):
            IndexMap.append(int(line.split()[0]))
    except:
        pass

#############  Extract parameters and variables from data blocks    ###############

ParaList=[]
WeightList=[]
DataList=[]
for block in Blocks:
    #extract the parameters list from all the data blocks
    blocklist=block.split("\n")
    #print blocklist[0].split()[-1]
    ParaList.append(tuple([float(blocklist[0].split()[i-1]) for i in Parameters]))
    #WeightList.append(float(blocklist[0].split()[-1]))
    WeightList.append(1.0)
    DataListTemp=[]
    #extract the datas asked by DataRank from all the data blocks
    
    for line in blocklist[1:]:
        try:
            if(int(line.split()[0]) in DataRankList):
                DataListTemp.append([float(line.split()[1])])
        except:
            pass
    DataList.append(DataListTemp)
# print ParaList
ParaSet=sorted(list(set(ParaList)))
# print DataList

#############   Merge data blocks  and create DataMap     ##################

#DataMap={Parameters:Data}
DataMap={}

for paras in ParaSet:
    Index=[]
    DataValue=[]
    datas=[]
    for i in range(0,len(ParaList)):
        if(ParaList[i]==paras):
            Index.append(i)
    WeightTemp=[WeightList[j] for j in Index]
    DataTemp=[DataList[j] for j in Index]
    for i in range(0,len(DataTemp[0])):
        mean=0.0
        std=0.0
        for j in range(0,len(DataTemp)):
            mean+=DataTemp[j][i][0]
            std+=DataTemp[j][i][1]**2*WeightTemp[j]
        mean=mean/len(DataTemp)
        std=(std/len(DataTemp)/sum(WeightTemp))**0.5
        datas.append([mean,std])
    # print datas,DataTemp
    DataMap[paras]=datas
    # print DataMap


###########    Output DataMap as the format             ###################
#os.remove("./*.ext")

Switch=True   # Classify data by variables
#Switch=False  #Classify data by parameters and variables


SubParameter=[]
for elem in ParaSet:
    SubParameter.append(elem[0])
SubParameter=list(set(SubParameter))

os.system("rm *.ext*.dat")
for para in SubParameter:
    SubDataMap={}
    ParaInFileList=[]
    for elem in ParaSet:
        if(Switch or elem[0]==para):
            #if(elem[0]==elem[1]):     #add some filter here
                paras=(elem[0],elem[1])
                SubDataMap[paras]=DataMap[elem]

##    for i in range(0,len(DataRank.values())):
##        if(not Switch):
##            OutputName=DataRank.values()[i]+"_"+str(para)+".ext.dat"
##        else:
##            OutputName=DataRank.values()[i]+".ext"
##        f=open(OutputName,"w")
##        for elem in sorted(SubDataMap.keys()):
##            OutStr="   ".join([str(e) for e in elem])+"   "+str(SubDataMap[elem][i][0])+"   "+str(SubDataMap[elem][i][1])+"\n"
##            f.write(OutStr)
##        f.close()
    # print IndexMap
    for i in range(0,len(IndexMap)):
        if(not Switch):
            OutputName=DataRank[IndexMap[i]]+"_"+str(para)+".ext.dat"
        else:
            OutputName=DataRank[IndexMap[i]]+".extcom.dat"
        f=open(OutputName,"w")
        for elem in sorted(SubDataMap.keys()):
            # if i==2:
                # print DataRank[IndexMap[i]], elem, SubDataMap[elem][i]
            # print elem, DataRank[IndexMap[i]], IndexMap[i]
            print SubDataMap
            OutStr="   ".join([str(e) for e in elem])+"   "+str(SubDataMap[elem][i][0])+"   "+str(SubDataMap[elem][i][1])+"\n"
            # print OutStr
            f.write(OutStr)
        f.close()

