import AnalysisTools as AT
import os
import numpy as np
print "Reading generation probes"





T = 1.94
H = 0.025
depth = 1.
L = 5.
folder = "../output"
os.chdir(folder)
dataW = AT.readProbeFile("pressure_gaugeArray.csv")

print dataW[1]
Z= -depth + dataW[1][0][1]

Nwaves = 10

Tend = dataW[2][-1]
Tstart = Tend-Nwaves*T

print Tstart,Tend,Z

Ht = 0
Hi = 0
bf = 2.

zc =[]
for dd in range(0,len(dataW[3][0,:])):
    dat = AT.zeroCrossing(dataW[2],dataW[3][:,dd],Tstart, Tend,minfreq=1/(bf*T),maxfreq=(bf/T))
    dat[1]=AT.pressureToHeight(dat[1],Z,depth,L,998.2,9.81)
    zc.append(dat[1]/0.025)
    print "x= ",dataW[1][dd][0], " H= ", dat[1]

