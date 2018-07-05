from builtins import str
from builtins import zip
from builtins import range
from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv
import os
import matplotlib.pyplot as plt

#####################################################################################

## Reading probes into the file
folder = "../output"
os.chdir(folder)
filename='combined_column_gauge.csv'

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

#####################################################################################

# Extracts the datas from the function readProbeFile  
datalist = readProbeFile(filename)
time = datalist[2]

# Calculates the time-average discharge over the crest
U = []
for i in range(0,len(datalist[3])):
    U.append(np.mean(datalist[3][i]))
U = np.array(U)
Q = U*0.25

#####################################################################################

# Plotting the probes
plt.plot(time, Q)
plt.xlabel('time [sec]')
plt.ylabel('Q [m^2/s]')
plt.suptitle('Time-averaged discharge under the gate of the sluice gate')
plt.ylim((0.6,1.4))
plt.xlim((-1.0,30.0))
plt.grid(True)
plt.savefig('CrumpWeir_discharge.png')
plt.show()

#####################################################################################

# Validation of the result
Q_th = 1.037 #Theoretical discharge between 20 s and 30 s
T = time.tolist()
T_20 = T.index(20.0)
T_30 = T.index(30.0)
T_20_to_30 = []
for i in range(T_20,T_30):
    T_20_to_30.append(Q[i])
Q_pr = np.mean(T_20_to_30) #Discharge between 20 s and 30 s obtained with PROTEUS
err = 100*abs(Q_th-Q_pr)/Q_th
val = open('validation_discharge_sluice.txt', 'w')
val.write('Gauges taken under the gate.'+'\n')
val.write('Average discharge between 20 s and 30 s.'+'\n')
val.write('Theory'+'\t'+'Simulation'+'\t'+'Error (%)'+'\n')
val.write(str(Q_th)+'\t'+str(Q_pr)+'\t'+str(err))
val.close()



