from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
#import DingemansWaveShoaling as dws

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
        probeCoord = list(zip(np.array(probex),np.array(probey),np.array(probez)))
        datalist = [probeType,probeCoord,time,data]
        return datalist

data_vof = readProbeFile(file_vof)

#####################################################################################

# Exctracting probes
time = data_vof[2]
vof = data_vof[3]

ETA = []
tank_dim = [58., 1.26] #dws.opts.tank_dim
waterLevel = 0.86 #dws.opts.water_level
gauge_x = [6.26, 10.26, 12.66, 23.26, 27.26, 29.26, 31.26, 33.66, 36.86, 40.26, 44.26] #dws.opts.gauge_x

for j in range(old_div(len(data_vof[1]),2)):
    eta = []
    for i in range(len(vof)):
        eta.append(tank_dim[1]-vof[:,j][i]-waterLevel)
    ETA.append(eta)
ETA = np.array(ETA)

#####################################################################################

# Plotting the probes
fig = plt.figure(figsize=(25,15))
ax = ['' for x in range(old_div(len(data_vof[1]),2))]
for i in range(old_div(len(data_vof[1]),2)):
    ax[i] = fig.add_subplot(4,3,i+1)
    ax[i].plot(time, ETA[i], 'r')
    ax[i].set_ylim([-0.02,0.02])
    ax[i].tick_params(labelsize=10)
    ax[i].set_title('Eta [m] against time [sec] at x='+str(gauge_x[i]), color='b', fontsize=12)
    ax[i].grid()
plt.tight_layout() 
plt.savefig('eta.png')   
#plt.show()
