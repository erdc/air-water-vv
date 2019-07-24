from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Put relative path below
filename='p_u_gauge.csv'

# Reading file
with open (filename, 'rb') as csvfile:
    data=csv.reader(csvfile, delimiter=",")
    a=[]
    time=[]
    probes=[]
    nRows=0

    for row in data:
# Time steps
        if nRows!=0:
            time.append(float(row[0]))

# Probes location ######################
        if nRows==0:                   #
            for i in row:              #
                if i!= '      time':   #
                    i=float(i[14:24])  #
                    probes.append(i)   #
########################################

        row2=[]
        for j in row:
            if j!= '      time' and nRows>0.:
                j=float(j)
                row2.append(j)
        a.append(row2)
        nRows+=1

#####################################################################################

### ---- PRESSURE ---- #
##    x_p = 1 #[1:39]
##    pressure = []
##    for k in range(1,nRows):
##        pressure.append(a[k][x_p])  
### Plot pressure in time   
##    plt.plot(time,pressure)
##    plt.xlabel('time [sec]')    
##    plt.ylabel('pressure [Pa]')
##    plt.suptitle('Pressure against time')
##    plt.grid(True)
##    plt.show()
##    savefig('Pressure_at_left.png')

# ---- VELOCITY ---- #
    y = 0.0
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for j in range(40,78)
        u = []
        for k in range(1,nRows):
            u.append(a[k][j])
        ax.plot(time, u, y)
        y = y+0.021875
    plt.xlabel('time [sec]')    
    plt.ylabel('velocity [m/s]')
    plt.zlabel('height [m]')
    plt.suptitle('Velocity')
    plt.grid(True)
    plt.show()
    savefig('Velocities_at_left.png')
    
#####################################################################################

### Print an output file to validate the results
##    maxPressureCal = max(pressure2)
##    maxPressureRef = 999.53
##    err = 100*abs(maxPressureRef-maxPressureCal)/maxPressureRef
##    val = open('validation_pressure.txt', 'w')
##    val.write('Only for gauges taken at (x,y)=(0.5,0.5)'+'\n')
##    val.write('Maximum pressure:'+'\n')
##    val.write('Reference'+'\t'+'Simulation'+'\t'+'\t'+'Error'+'\n')
##    val.write(str(maxPressureRef)+'\t'+str(maxPressureCal)+'\t'+str(err))
##    val.close()
    
#####################################################################################

### Print an output file
##    info = open('probes.txt','w')
##    string1=str(probes)
##    string2=string1.replace('[',' ')
##    string3=string2.replace(']',' ')   
##    string4=string3.replace(',','\t') 
##    info.write(str('x')+'\t'+string4+'\n')
##    for j in range(1,nRows):
##        string5=str(a[j])
##        string6=string5.replace('[',' ')
##        string7=string6.replace(']',' ')   
##        string8=string7.replace(',','\t')   
##        info.write(string8+'\n')
##    info.close()




