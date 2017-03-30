import AnalysisTools as AT
import os
import numpy as np
print "Reading generation probes"





T = 3. 
H = 0.15
depth = 1.
L = 8.74
folder = "../output"
os.chdir(folder)
dataW = AT.readProbeFile("pressure_gaugeArray.csv")

dataW[1] = dataW[1][::2]
print dataW[1]

#Z= -depth + dataW[1][0][1]

Nwaves = 10

Tend = dataW[2][-1]
Tstart = Tend-Nwaves*T

print Tstart,Tend

Ht = 0
Hi = 0
bf = 2.5

zc =[]
for dd in range(0,len(dataW[3][0,:])):
    dat = AT.zeroCrossing(dataW[2],dataW[3][:,dd],Tstart, Tend,minfreq=1/(2*T),maxfreq=(bf/T))
#    dat[1]=AT.pressureToHeight(dat[1],Z,depth,L,998.2,9.81)

    zc.append(dat)

dx_array = dataW[1][1][0]-dataW[1][0][0]          
#print dx_array
Narray =int(round(L/6/dx_array))
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

print HH[:]
print np.mean(RR[:])
