#!/usr/bin/env python
import sys
def nsf(num, n=1):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n-1)).format(num)
    return float(numstr)
newton_its=0
with open(sys.argv[1],'r') as f:
   for line in f.readlines():
       if "Iteration Statistics" in line:
           print "Average Iteration Statistics\n==================================="
       elif "out_" in line:
           if newton_its > 0:
#               tot=avg
               avg/=float(newton_its)
               avg11/=float(tot)
               print "average ksp its ",avg," over ", newton_its," newton iterations"
               print "max ksp its ",max_its," over ", newton_its," newton iterations"
               print "failed ksp its ",fail," over ", newton_its," newton iterations"
               print "average A11 ksp its ",avg11," over ", tot," ksp iterations"
               print "max A11 ksp its ",max_its11," over ", tot," ksp iterations"
               print "failed A11 ksp its ",fail11," over ", tot," ksp iterations\n"
               print "time taken ",time
               print "flops taken ",flops
               print "KSP solve percentage time ",ksppt
               print "KSP A11 solve percentage time ",ksp11pt
               print "KSP S solve percentage time ",kspSpt
               print "PCApply percentage time ",pcapply
               print "PCSetUp solve percentage time ",pcsetup,"\n"
               orig_stdout = sys.stdout
               w = open('latex_results.txt', 'a')
               sys.stdout = w
               print fname
               print "& "+str(nsf(avg,4))+"/"+str(max_its)+" & "+str(nsf(avg11,4))+"/"+str(max_its11)
               print "& "+str(time)+" & "+str(ksppt)+"\% & "+str(ksp11pt)+"/"+str(kspSpt)+"\%"
               sys.stdout = orig_stdout
               w.close()
               print line[:-1]
           else:
               print line[:-1]
           fnum=line[-7:-5]
           fname=line[17:-5]
           avg=0.0
           max_its = 0
           newton_its=0
           avg11=0.0
           max_its11=0
           tot=0
           fail=0
           fail11=0
       elif "rans2p_ solve converged" in line:
           words = line.split()
           for i,w in enumerate(words):
               if w == "iterations":
                   avg+=float(words[i+1])
                   newton_its+=1
                   if int(words[i+1]) > max_its:
                       max_its = int(words[i+1])
       elif "rans2p_ solve did not converge" in line:
           fail+=1
       elif "rans2p_fieldsplit_velocity_ solve converged" in line:
           words = line.split()
           for i,w in enumerate(words):
               if w == "iterations":
                   avg11+=float(words[i+1])
                   tot+=1
                   if int(words[i+1]) > max_its11:
                       max_its11 = int(words[i+1])
       elif "rans2p_fieldsplit_velocity_ did not converge" in line:
           fail11+=1
       elif "Time (sec):" in line:
           words = line.split()
           for i,w in enumerate(words):
               if w == "(sec):":
                   time = float(words[i+1])
                   if time >= 1000:
                       time = int(time)
       elif "Flops:" in line:
           words = line.split()
           for i,w in enumerate(words):
               if w == "Flops:":
                   flops = int(float((words[i+1])))
       elif "PCSetUp" in line:
           words = line.split()
           pcsetup = int(words[10])
       elif "PCApply" in line:
           words = line.split()
           pcapply = int(words[10])
       elif "KSPSolve_FS_0" in line:
           words = line.split()
           ksp11pt = int(words[10])
       elif "KSPSolve_FS_Schu" in line:
           words = line.split()
           kspSpt = int(words[10])
       elif "KSPSolve " in line:
           words = line.split()
           ksppt = int(words[10])
   if newton_its > 0:
#       tot=avg
       avg/=float(newton_its)
       avg11/=float(tot)
       print "average ksp its ",avg," over ", newton_its," newton iterations"
       print "max ksp its ",max_its," over ", newton_its," newton iterations"
       print "failed ksp its ",fail," over ", newton_its," newton iterations"
       print "average A11 ksp its ",avg11," over ", tot," ksp iterations"
       print "max A11 ksp its ",max_its11," over ", tot," ksp iterations"
       print "failed A11 ksp its ",fail11," over ", tot," ksp iterations\n"
       print "time taken ",time
       print "flops taken ",flops
       print "KSP solve percentage time ",ksppt
       print "KSP A11 solve percentage time ",ksp11pt
       print "KSP S solve percentage time ",kspSpt
       print "PCApply percentage time ",pcapply
       print "PCSetUp solve percentage time ",pcsetup
       orig_stdout = sys.stdout
       w = open('latex_results.txt', 'a')
       sys.stdout = w
       print fname
       print "& "+str(nsf(avg,4))+"/"+str(max_its)+" & "+str(nsf(avg11,4))+"/"+str(max_its11)
       print "& "+str(time)+" & "+str(ksppt)+"\% & "+str(ksp11pt)+"/"+str(kspSpt)+"\%"
       sys.stdout = orig_stdout
       w.close()
