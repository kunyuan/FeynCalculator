#!/usr/bin/python
import random
import os, sys
IsCluster=False
sourcedir="./"
execute="bubble"
compiler="ifort"
# compiler="gfortran"
filelist=os.listdir(sourcedir)
sourcename=[elem for elem in filelist if elem[0:6]=="bubble" and elem[-3:]=="f90"]
print sourcename
sourcename.sort()
sourcename=sourcename[-1]
tothomedir=os.getcwd()
inlist=open(tothomedir+"/inlist","r")

i=0
NBlock=1024
#isload=int(inlist.readline().split(":")[1])
seed=-int(random.random()*100000)

for eachline in inlist:
        i+=1
        para=eachline.split()

        if len(para)==0:
            print "All submitted!"
            break

        homedir=os.getcwd()+"/Beta"+str(para[0])+"_rs"+str(para[1])
        if(os.path.exists(homedir)!=True):
            os.system("mkdir "+homedir)

        os.system("cp DiagPolar*.txt "+homedir)
        os.system("mkdir "+homedir+"/Data")

        os.system(compiler+" "+sourcedir+"/"+sourcename+" -O3 -o "+homedir+"/"+execute)
        infilepath=homedir+"/infile"
        outfilepath=homedir+"/outfile"
        jobfilepath=homedir+"/jobfile"
        if(os.path.exists(infilepath)!=True):
            os.system("mkdir "+infilepath)
        if(os.path.exists(outfilepath)!=True):
            os.system("mkdir "+outfilepath)
        if(os.path.exists(jobfilepath)!=True):
            os.system("mkdir "+jobfilepath)

        for j in range(int(para[-1])):
            seed-=1
            infile="_in"+str(i)+"_"+str(j)
            outfile="_out"+str(i)+"_"+str(j)
            jobfile="_job"+str(i)+"_"+str(j)+".sh"
            f=open(infilepath+"/"+infile,"w")
            item=para[0:-1]
            item.append(str(seed))
            # item.append(str(NBlock))
            item.append(str(j))
            stri=" ".join(item)
            f.write(stri)
            f.close()

            if IsCluster==False:
                os.chdir(homedir)
                os.system("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile+" &")
                os.chdir("..")
                #os.system("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile+" &")
                #os.system("./"+execute+" < "+infilepath+"/"+infile)
            else:
                with open(jobfilepath+"/"+jobfile, "w") as fjob:
                    fjob.write("#!/bin/sh\n"+"#PBS -N "+jobfile+"\n")
                    fjob.write("#PBS -o "+homedir+"/Output\n")
                    fjob.write("#PBS -e "+homedir+"/Error\n")
                    fjob.write("#PBS -l walltime=2000:00:00\n")
                    fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
                    fjob.write("cd "+homedir+"\n")
                    fjob.write("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile)

                os.system("qsub "+jobfilepath+ "/"+jobfile)
                os.system("rm "+jobfilepath+ "/"+jobfile)
print("Jobs manage daemon is ended")
sys.exit(0)
