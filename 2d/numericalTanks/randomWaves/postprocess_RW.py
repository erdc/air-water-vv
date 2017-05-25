import collections as cll
import csv
import os
import matplotlib.pyplot as plt
import random_waves as rw
import numpy as np

#####################################################################################

## Reading probes into the file
folder = "output"
os.chdir(folder)
file_vof = 'column_gauges.csv'
file_p = 'pressure_gaugeArray.csv'
file_s='../series.txt'
eta_v= np.loadtxt(file_s)

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
tank_dim = rw.opts.tank_dim
waterLevel = rw.opts.water_level
eta_fast = np.array(tank_dim[1]-vof[:,0]-waterLevel)

# Eta for RandomWaves
Tp = rw.Tp
Hs = rw.Hs
mwl = rw.mwl
depth = rw.depth
waveDir = np.array(rw.waveDir)
g = np.array(rw.g)
N = rw.N
bandFactor = rw.bandFactor
spectName =  rw.spectName
phi = rw.phi
wave_ref = rw.wt.RandomWaves(Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params=None,phi=phi,fast=True)
eta_bc = rw.tank.BC['x-'].vof_dirichlet.uOfXT
zin = np.linspace(rw.he/2,rw.tank_dim[1]-rw.he/2,rw.tank_dim[1]/rw.he)

Tstart = rw.Tstart
Tend = rw.Tend
x0 = rw.x0
Lgen = np.array([0., 0., 0.])

X = np.array([0., 0., 0.])
eta_ref = []
eta_bca = []

for i in range(0,len(time)):
    eta_ref.append(wave_ref.eta(X,time[i]))
    aa = 0.
    for ii in range(len(zin)):
        xx = np.array([x0[0],zin[ii],x0[1]])
        aa += rw.he*eta_bc(xx,time[i])
    eta_bca.append(tank_dim[1]-waterLevel -aa)

eta_ref = np.array(eta_ref)
eta_bca = np.array(eta_bca)

#####################################################################################

# Plotting the probes
plt.figure(num='eta')
plt.plot(time, eta_fast, 'b', label='End of RZ')
plt.plot(time, eta_ref, 'r--', label='RandomWaves (calculated)')
plt.plot(eta_v[:,0], eta_v[:,1], 'y--', label='RandomWaves (printed)')
plt.plot(time, eta_bca, 'g:', label='setUnsteady (bc module)')
plt.legend(loc='best')
plt.xlabel('time [sec]')
plt.ylabel('eta [m]')
plt.ylim(-0.08,0.08)
plt.suptitle('Surface elevation against time for the random waves')
plt.grid()
plt.show()
plt.savefig('eta_RW.png')

#####################################################################################

# Validation of the results

S = 0.
c = 0.
for i in range(len(time)):
    c += 1.
    S += (eta_fast[i]-eta_ref[i])**2
err = np.sqrt(S/c)
err = 100*err/(rw.opts.Hs)
val = open('validation_eta_RW.txt', 'w')
val.write('Surface elevation against time for the random waves'+'\n')
val.write('Gauges taken between 0s and 30s'+'\n')
val.write('Average error (%) between the theoretical function and the simulation:'+'\n')
val.write(str(err))
val.close()


