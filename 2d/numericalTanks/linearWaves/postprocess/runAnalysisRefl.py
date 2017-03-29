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

Nwaves = 3

Tend = dataW[2][-1]
Tstart = Tend-Nwaves*T

print Tstart,Tend,Z

Ht = 0
Hi = 0
bf = 1.2

zc =[]
for dd in range(0,len(dataW[3][0,:])):
    dat = AT.zeroCrossing(dataW[2],dataW[3][:,dd],Tstart, Tend,minfreq=1/(bf*T),maxfreq=(bf/T))
    print "[period,pressure] = ", dat
    dat[1]=AT.pressureToHeight(dat[1],Z,depth,L,998.2,9.81)
    print "[period,height] = ", dat

    zc.append(dat)

dx_array = dataW[1][1][0]-dataW[1][0][0]          
print dx_array
Narray =int(round(L/6./dx_array))
RR = []
HH=[]
for dd in range(0,len(dataW[3][0,:])-2*Narray):
    dd1 = dd + Narray
    dd2 = dd1 + Narray
    H1 = zc[dd][1]
    H2 = zc[dd1][1]
    H3 = zc[dd2][1]
    HH.append(AT.reflStat(H1,H2,H3,Narray*dx_array,L)[0])
    RR.append(AT.reflStat(H1,H2,H3,Narray*dx_array,L)[2])

print np.mean(HH[21:])
print np.mean(RR[21:])
