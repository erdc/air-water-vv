#!/usr/bin/env python
import sys
newton_its=0.0
with open(sys.argv[1],'r') as f:
   for line in f.readlines():
       if "Iteration Statistics" in line:
           print "Average Iteration Statistics\n==================================="
       elif "out_" in line:
           if newton_its > 0.0:
               tot=avg
               avg/=newton_its
               avg11/=tot
               print "average ksp its ",avg," over ", newton_its," newton iterations"
               print "max ksp its ",max_its," over ", newton_its," newton iterations"
               print "average A11 ksp its ",avg11," over ", tot," ksp iterations"
               print "max A11 ksp its ",max_its11," over ", tot," ksp iterations\n"
               print line[:-1]
           else:
               print line[:-1]
           avg=0.0
           max_its = 0
           newton_its=0.0
           avg11=0.0
           max_its11=0.0
       elif "rans2p_ solve" in line:
           words = line.split()
           for i,w in enumerate(words):
               if w == "iterations":
                   avg+=float(words[i+1])  
                   newton_its+=1.0
                   if int(words[i+1]) > max_its:
                       max_its = int(words[i+1])
       elif "rans2p_fieldsplit_velocity_ solve" in line:
           words = line.split()
           for i,w in enumerate(words):
              if w == "iterations":
                  avg11+=float(words[i+1])
                  if int(words[i+1]) > max_its11:
                      max_its11 = int(words[i+1])
   if newton_its > 0.0:
       tot=avg
       avg/=newton_its
       avg11/=tot
       print "average ksp its ",avg," over ", newton_its," newton iterations"
       print "max ksp its ",max_its," over ", newton_its," newton iterations"
       print "average A11 ksp its ",avg11," over ", tot," ksp iterations"
       print "max A11 ksp its ",max_its11," over ", tot," ksp iterations"
