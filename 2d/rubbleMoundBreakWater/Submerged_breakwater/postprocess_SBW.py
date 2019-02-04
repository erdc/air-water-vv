import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import submerged_breakwater as sbw
from AnalysisTools import zeroCrossing

#####################################################################################

## Reading probes into the file
folder = "output"
os.chdir(folder)
file_vof = 'line_integral_gauges_1.csv'

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
ETA = []
tank_dim = sbw.tank_dim #[13.7815, 0.75] 
waterLevel = sbw.waterLevel #0.315
gauge_x = []
for k in range(len(sbw.columnLines1)):
    gauge_x.append(sbw.columnLines1[k][0][0])
for j in range(len(data_vof[1])/2):
    eta = []
    for i in range(len(vof)):
        eta.append(tank_dim[1]-vof[:,j][i]-waterLevel)
    ETA.append(eta)
ETA = np.array(ETA)

#####################################################################################

# Plotting the probes
fig = plt.figure(figsize=(25,15))
ax = ['' for x in range(len(data_vof[1])/2)]
for i in range(len(data_vof[1])/2):
    ax[i] = fig.add_subplot(6,4,i+1)
    ax[i].plot(time, ETA[i], 'r')
    ax[i].set_ylim([-0.06,0.08])
    ax[i].tick_params(labelsize=10)
    ax[i].set_title('Eta [m] against time [sec] at x='+str(gauge_x[i]), color='b', fontsize=12)
    ax[i].grid()
plt.tight_layout() 
plt.savefig('eta.png')   
#plt.show()

#####################################################################################


#print time[300:]

time1 = time[900:]
ETA1 = ETA[1]
ETA2 = ETA1[900:]

ETA3 = ETA[16]
ETA4 = ETA3[900:]


#a = len(ETA)
# Transmission coefficient
zc = []
#for i in range(len(ETA)):
zc=zeroCrossing(time1,ETA2)
zc1=zeroCrossing(time1,ETA4)
zc = np.array(zc)
zc1 = np.array(zc1)
K = zc1[1]/zc[1]
#K = np.mean(zc1[1])/sbw.opts.wave_height
print 'Transmission coefficient'+'\t'+'='+'\t'+str(K)
print zc
print zc1
#print ETA[17,7]
#print len(time)
