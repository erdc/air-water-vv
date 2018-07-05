from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import collections as cll
import csv
import os
import matplotlib.pyplot as plt
import random_waves as rw
import numpy as np
from AnalysisTools import signalFilter,zeroCrossing,reflStat


#####################################################################################

## Reading probes into the file
folder = "RWFast"
os.chdir(folder)
file_vof = 'column_gauges.csv'
file_p = 'pressure_gaugeArray.csv'
file_s='../series.txt'
#eta_v= np.loadtxt(file_s)
def readProbeFile(filename):
    with open (filename, 'rb') as csvfile:
        data = np.loadtxt(csvfile, delimiter=",",skiprows=1)
        time = data[:,0]
        data = data[:,1:]
        csvfile.seek(0)
        header = csvfile.readline()
        header = header.replace("time","")
        header = header.replace("[","")
        header = header.replace("]","")
        header = header.replace(","," ")
        header = header.split()
        probeType = []
        probex = []
        probey = []
        probez = []        
        for ii in range(0,len(header),4):
            probeType.append(header[ii])
            probex.append(float(header[ii+1]))
            probey.append(float(header[ii+2]))
            probez.append(float(header[ii+3]))
        probeCoord = list(zip(np.array(probex),np.array(probey),np.array(probez)))
        datalist = [probeType,probeCoord,time,data]
        return datalist

data_vof = readProbeFile(file_vof)
data_p = readProbeFile(file_p)

#####################################################################################
# Eta for RandomWavesFast
time = data_vof[2]
vof = data_vof[3]
tank_dim = rw.opts.tank_dim
waterLevel = rw.opts.water_level
eta_fast = np.array(tank_dim[1]-vof[:,0]-waterLevel)
#####################################################################################
# Eta for RandomWaves
Tp = rw.Tp
Hs = rw.Hs
mwl = rw.mwl
depth =  rw.depth
waveDir = np.array(rw.waveDir)
g = np.array([0,-9.81,0.])#(rw.g)
N =  rw.N
bandFactor = rw.bandFactor
spectName =  rw.spectName
phi = rw.phi
wave_ref = rw.wt.RandomWaves(Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params=None,phi=phi,fast=True)
eta_bc= rw.tank.BC['x-'].vof_dirichlet.uOfXT
zin = np.linspace(old_div(rw.opts.he,2.),rw.tank_dim[1]-old_div(rw.opts.he,2),old_div(rw.opts.tank_dim[1],rw.opts.he))
######################################################################################

Tstart =0.# rw.Tstart
Tend = rw.Tend
x0 = rw.x0
Lgen = np.array([5., 0., 0.])
#wave_fast = wt.RandomWavesFast(Tstart,Tend,x0,Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,
#                               spectral_params=None,phi=phi,Lgen=Lgen,Nwaves=30,Nfreq=64)
X = np.array([0., 0., 0.])
eta_ref = []
eta_bca = [ ]
#eta_fast = []
for i in range(0,len(time)):
    eta_ref.append(wave_ref.eta(X,time[i]))
    aa = 0.
    for ii in range(len(zin)):
        xx =np.array([x0[0],zin[ii],x0[1]])
        aa+=rw.opts.he*eta_bc(xx,time[i])
    eta_bca.append(tank_dim[1]-waterLevel -aa)
        
        
    #eta_fast.append(wave_fast.eta(X,time[i]))
eta_ref = np.array(eta_ref)
eta_bca = np.array(eta_bca)
fp = old_div(1.,Tp)
minf = old_div(fp,bandFactor)
maxf = bandFactor*fp
time_int = np.linspace(0,time[-1],len(time))
eta_fast = np.interp(time_int,time,eta_fast)
eta_fast = signalFilter(time,eta_fast,minf, maxf, 1.1*maxf, 0.9*minf)

#eta_fast = np.array(eta_fast)

#####################################################################################

# Plotting the probes
plt.figure(num='eta')
plt.plot(time_int, eta_fast, 'b', label='End of RZ')
plt.plot(time, eta_ref, 'r--', label='RandomWaves (calculated)')
#plt.plot(eta_v[:,0], eta_v[:,1], 'y--', label='RandomWaves (printed)')
plt.plot(time, eta_bca, 'g:', label='setUnsteady (bc module)')
plt.legend(loc='best')
plt.xlabel('time [sec]')
plt.ylabel('eta [m]')
plt.ylim(-0.08,0.08)
plt.xlim(0.,20.)
plt.suptitle('Surface elevation against time for the random waves')
plt.grid()
plt.show()
plt.savefig('eta_RW.png')
