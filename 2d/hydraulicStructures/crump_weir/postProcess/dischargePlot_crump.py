from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv
import os
import matplotlib.pyplot as plt

#####################################################################################

## Reading probes into the file
#folder = "../output"
folder = "../results_delete"
os.chdir(folder)
file_vof='column_gauge.csv'
file_u='u_over_crest.csv'

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
        header =  header.replace(","," ")
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
data_u = readProbeFile(file_u)

#####################################################################################

# Solves the problem if velocity and vof have not the same number of time steps
if len(data_u[2])==len(data_vof[2]) or len(data_u[2])<len(data_vof[2]):
    time = data_u[2]
else:
    time = data_vof[2]

# Extracts the datas from the function readProbeFile   
U = []
VOF = []
water_level = []
for i in range(0,len(data_u[3])):
    U.append(data_u[3][i])
for i in range(0,len(data_vof[3])):
    VOF.append(data_vof[3][i][16])
    water_level.append(3.75-data_vof[3][i][16])
U = np.array(U)   
VOF = np.array(VOF)
water_level = np.array(water_level)
u_coord = np.array(data_u[1])

# Changes the order of the velocity
tmp = []
u = []
for i in range(0,len(U[0])):
    for j in range(0,len(U)):
        tmp.append(U[j][i])
    u.append(tmp)
    tmp = []
u = np.array(u)

# Creates an average velocity of the water at each time step over the crest
vel = []
for i in range(0,len(time)):
    for j in range(0,len(u_coord)):
        if u_coord[j][1]<water_level[i]:
            tmp.append(u[j][i])
    vel.append(np.mean(tmp))
    tmp = []
vel = np.array(vel)

# Calculates the time-average discharge over the crest
section = []
Q = []
for i in range(0,len(time)):
    section.append(water_level[i]-0.5)
    Q.append(section[i]*vel[i])

#####################################################################################

## Plotting the probes
fig = plt.plot(time, Q)
plt.xlabel('time [sec]')
plt.ylabel('Q [m^2/s]')
plt.suptitle('Flow discharge at the crest of the crump weir')
#plt.ylim((0.,3.5))
#plt.xlim((0.,30.))
plt.grid(True)
plt.savefig('SluiceGate_discharge.png')
plt.show()

#####################################################################################

# Validation of the result
Q_th = 2.0175 #Theoretical discharge between 20 s and 30 s
T = time.tolist()
T_20 = T.index(0.1)
T_30 = T.index(0.8)
T_20_to_30 = []
for i in range(T_20,T_30):
    T_20_to_30.append(Q[i])
Q_pr = np.mean(T_20_to_30) #Discharge between 20 s and 30 s obtained with PROTEUS
err = 100*abs(Q_th-Q_pr)/Q_th
val = open('validation_discharge_crump.txt', 'w')
val.write('Gauges taken over the crest of the crump weir.'+'\n')
val.write('Average discharge between 20 s and 30 s.'+'\n')
val.write('Theory'+'\t'+'Simulation'+'\t'+'Error (%)'+'\n')
val.write(str(Q_th)+'\t'+str(Q_pr)+'\t'+str(err))
val.close()









