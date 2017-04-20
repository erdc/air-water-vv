import numpy as np
import collections as cll
import csv
import os
import matplotlib.pyplot as plt
import linear_waves as lw
from proteus import WaveTools as wt

#####################################################################################

## Reading probes into the file
folder = "output"
os.chdir(folder)
file_vof = 'column_gauges.csv'

def readProbeFile(filename):
    with open (filename, 'rb') as csvfile:
        data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
        time=data[:,0]
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

#####################################################################################

# Exctracting probes
time = data_vof[2]
vof = data_vof[3]
eta_num = []
tank_dim = lw.opts.tank_dim
waterLevel = lw.opts.water_level
i_mid = len(vof[0])/2-1
for i in range(0, len(vof)):
    eta_num.append(tank_dim[1]-vof[i][i_mid]-waterLevel)
eta_num = np.array(eta_num)

# Theoretical eta
x = np.array(data_vof[1][2*i_mid])
wave = lw.wave
eta_th = []
for i in range(0,len(time)):
    eta_th.append(wave.eta(x,time[i]))
 
#####################################################################################

# Plotting the probes
plt.figure(num='eta')
plt.plot(time, eta_num, 'b', label='numerical')
plt.plot(time, eta_th, 'r--', label='theoretical')
plt.legend(loc='upper right')
plt.xlabel('time [sec]')
plt.ylabel('eta [m]')
plt.xlim((0.,30.))
plt.ylim((-0.1,0.1))
plt.suptitle('Surface elevation against time in the middle of the tank.')
plt.grid()
plt.show()
plt.savefig('eta_LW.png')

#####################################################################################

# Validation of the result
S = 0.
c = 0.
istart = np.where(time>=6.)[0][0]
iend = np.where(time>=18.)[0][0]
for i in range(istart,iend):
    c = c + 1.
    S = S + (eta_th[i]-eta_num[i])**2
err = np.sqrt(S/c)
err = 100*err/(lw.opts.wave_height+waterLevel)
val = open('validation_eta_NLW.txt', 'w')
val.write('Eta in the middle of the tank.'+'\n')
val.write('Gauges taken between 6s and 18s'+'\n')
val.write('Average error (%) between the theoretical function and the simulation:'+'\n')
val.write(str(err))
val.close()
