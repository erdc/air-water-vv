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
file_pressurePoint = 'pressure_PointGauge.csv'
file_pressureLine = 'pressure_LineGauge.csv'
file_vof = 'vof_LineIntegralGauge.csv'

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

data_p_point = readProbeFile(file_pressurePoint)
data_p_line = readProbeFile(file_pressureLine)
data_vof = readProbeFile(file_vof)

#####################################################################################

# Exctracting probes
time = data_p_point[2]
P_point = data_p_point[3]
P_line = data_p_line[3]
vof = data_vof[3]
water_level = []
for i in range(0,len(vof)):
    water_level.append(1.8-vof[i])
Y = []
for i in range(0,len(data_p_line[1])):
    Y.append(data_p_line[1][i][1])
Y = np.array(Y)

#####################################################################################

# Definition of the theoretical pressure under and over water
rho_w = 998.2
rho_a = 1.205
g = 9.81
h = 0.6
H = 1.8
def p(y):
    if np.any(y<=0.6):
        return rho_a*g*(H-h) + rho_w*g*(h-y)
    else:
        return rho_a*g*(H-y)

#####################################################################################

# Plotting the probes

plt.figure(num='Water level')
plt.plot(time, water_level)
plt.xlabel('time [sec]')
plt.ylabel('Water level in the middle of the tank [m]')
plt.ylim((0.,1.8))
plt.savefig('water_level_QW.png')

plt.figure(num='Pressure point')
plt.plot(time, P_point)
plt.xlabel('time [sec]')
plt.ylabel('Pressure at (3.22, 0.6) [Pa]')
plt.ylim((2850,2980))
plt.savefig('pressure_point_QW.png')

plt.figure(num='Pressure line')
plt.plot(P_line[-1], Y, 'ro', label='Numerical')
plt.plot(p(Y), Y, 'b', label='Theoretical')
plt.legend(loc='upper right')
plt.xlabel('Pressure [Pa]')
plt.ylabel('Y position [m]')
plt.savefig('pressure_line_QW.png')

plt.show()

#####################################################################################

# Validation of the result

water_level_th = 0.6
wl = water_level
water_level_num = wl[-1][0]
err_wl = 100*abs(water_level_th-water_level_num)/water_level_th
val = open('validation_WaterLevel_QW.txt', 'w')
val.write('Water level in the middle of the tank.'+'\n')
val.write('Gauge taken after 1s.' + '\n')
val.write('Theory'+'\t'+'Simulation'+'\t'+'Error (%)'+'\n')
val.write(str(water_level_th)+'\t'+str(water_level_num)+'\t'+str(err_wl))
val.close()

press_point_th = p(0.3)
press_point_num = P_point[-1][0]
err_pp = 100*abs(press_point_th-press_point_num)/press_point_th
val1 = open('validation_PressurePoint_QW.txt', 'w')
val1.write('Pressure at the point (x,y) = (3.22, 0.3)'+'\n')
val1.write('Gauge taken after 1s.'+'\n')
val1.write('Theory'+'\t'+'\t'+'Simulation'+'\t'+'Error (%)'+'\n')
val1.write(str(press_point_th)+'\t'+str(press_point_num)+'\t'+str(err_pp))
val1.close()

S = 0.
for i in range(0,len(Y)-4): # Ignores the 4 last points beside the water surface
    S = S + abs(p(Y[i])-P_line[-1][i])/p(Y[i])
err_pl = 100*S/(len(Y)-4)
val2 = open('validation_PressureLine_QW.txt', 'w')
val2.write('Pressure under the water at a column in the middle of the tank.'+'\n')
val2.write('Gauges taken after 1s.'+'\n')
val2.write('Average error (%) between the theoretical function and the simulation:'+'\n')
val2.write(str(err_pl))
val2.close()







