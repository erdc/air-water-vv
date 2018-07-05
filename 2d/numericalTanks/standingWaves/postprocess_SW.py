from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
import numpy as np
import math
import csv
import os
import standing_waves as sw

#####################################################################################

# Reading probes into the file
folder = "output"
os.chdir(folder)
file_p = 'pressure_gaugeArray.csv'

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

data_p = readProbeFile(file_p)

#####################################################################################

# Exctracting probes
T = sw.opts.wave_period #1.94
H = sw.opts.wave_height #0.025
depth = sw.opts.water_level #1.
L = sw.opts.wavelength #5.
time = data_p[2]
pressure = data_p[3][:,-1]
Z = -depth + data_p[1][0][1]
Nwaves = old_div((sw.opts.tank_dim[0]+sw.opts.tank_sponge[0]+sw.opts.tank_sponge[1]),L)

# Calculating the height with the pressure

def pressureToHeight(data,Z,depth,wavelength,rho,g):
    k = 2*math.pi/wavelength
    Kp = rho*g*np.cosh(k*(depth+Z))/np.cosh(k*depth)
    return old_div(data,Kp)

Tend = time[-1]
Tstart = Tend-Nwaves*T
Tend_period = Tstart+T
istart = np.where(time>=Tstart)[0][0]
iend_period = np.where(time>Tend_period)[0][0]
p_period = pressure[istart:iend_period]
p_0 = max(p_period)-min(p_period)
Hr = pressureToHeight(p_0,Z,depth,L,998.2,9.81)

#####################################################################################

# Validation of the result
err = 100*abs(2*H-Hr)/(2*H)
val = open('validation_SW.txt', 'w')
val.write('Wave height for standing waves.'+'\n')
val.write('Gauge taken at the right side of the tank.'+'\n')
val.write('Theory'+'\t'+'Simulation'+'\t'+'Error (%)'+'\n')
val.write(str(2*H)+'\t'+str(Hr)+'\t'+str(err))
val.close()


