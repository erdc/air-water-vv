import numpy as np
import collections as cll
import csv
import os
import matplotlib.pyplot as plt
#import random_waves as rw
from proteus import WaveTools as wt
import math

#####################################################################################

## Reading probes into the file
folder = "output"
os.chdir(folder)
file_vof = 'column_gauges.csv'
file_p = 'pressure_gaugeArray.csv'

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
        probeCoord = zip(np.array(probex),np.array(probey),np.array(probez))
        datalist = [probeType,probeCoord,time,data]
        return datalist

data_vof = readProbeFile(file_vof)
data_p = readProbeFile(file_p)

#####################################################################################

# Eta for RandomWavesFast
time = data_vof[2]
vof = data_vof[3]
tank_dim = [15., 1.5]   # rw.opts.tank_dim
waterLevel = 1.   # rw.opts.water_level
eta_fast = np.array(tank_dim[1]-vof[:,0]-waterLevel)

# Eta for RandomWaves
Tp = 1.94  # rw.Tp
Hs = 0.05   # rw.Hs
mwl = 1.0   # rw.mwl
depth = 1.0 # rw.depth
waveDir = np.array([1.,0.,0.])   # rw.waveDir
g = np.array([0.,-9.81,0.])
N = 2000   # rw.N
bandFactor = 2.0   # rw.bandFactor
spectName = 'JONSWAP'   # rw.spectName
phi = np.loadtxt('phases.txt')
wave_ref = wt.RandomWaves(Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params=None,phi=phi,fast=True)
Tstart = 0.
Tend = 120.
x0 = np.array([0., 0. ,0. ])
Lgen = np.array([5., 0., 0.])
#wave_fast = wt.RandomWavesFast(Tstart,Tend,x0,Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,
#                               spectral_params=None,phi=phi,Lgen=Lgen,Nwaves=30,Nfreq=64)
X = np.array([0., 0., 0.])
eta_ref = []
#eta_fast = []
for i in range(0,len(time)):
    eta_ref.append(wave_ref.eta(X,time[i]))
    #eta_fast.append(wave_fast.eta(X,time[i]))
eta_ref = np.array(eta_ref)
#eta_fast = np.array(eta_fast)

# Pressure
rho = 998.2
a = wave_ref.ai
k = wave_ref.ki
omega = wave_ref.omega
z = -0.5
x = 0.
P_ref = []
for j in range(len(time)):
    S = 0.
    for i in range(N):
        S += rho*abs(g[1])*a[i]*math.cos(k[i]*x-omega[i]*time[j]+phi[i])*math.cosh(k[i]*(depth+z))/math.cosh(k[i]*depth)
    P_ref.append(S)
P_fast = data_p[3][:,0]

P_ref = np.array(P_ref)-np.mean(np.array(P_ref))
P_fast = np.array(P_fast)-np.mean(np.array(P_fast))

#####################################################################################

# Plotting the probes
plt.figure(num='eta')
plt.plot(time, eta_fast, 'b', label='RandomWavesFast')
plt.plot(time, eta_ref, 'r--', label='RandomWaves')
plt.legend(loc='best')
plt.xlabel('time [sec]')
plt.ylabel('eta [m]')
plt.suptitle('Surface elevation against time for the random waves')
plt.grid()
plt.show()
plt.savefig('eta_RW.png')

plt.figure(num='p')
plt.plot(time, P_fast, 'b', label='RandomWavesFast')
plt.plot(time, P_ref, 'r--', label='RandomWaves')
plt.legend(loc='best')
plt.xlabel('time [sec]')
plt.ylabel('pressure [Pa]')
plt.suptitle('Pressure against time for the random waves')
plt.grid()
plt.show()
plt.savefig('pressure_RW.png')




