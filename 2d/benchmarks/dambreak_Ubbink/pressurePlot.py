from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put relative path below
filename='pressureGauge.csv'

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

# Choose which probes to plot  
    x = 1
    pressure2=[]
    for k in range(1,nRows):
        pressure2.append(a[k][x])    
# Plot pressure in time
    import matplotlib.pyplot as plt
    plt.plot(time,pressure2)
    plt.xlabel('time step [sec]')    
    plt.ylabel('pressure [Pa]')
    plt.suptitle('Pressure against time at (x,y)=(0.292,0.04)')
    plt.grid(True)
    plt.show()
    savefig('pressure_in_time.png')

#####################################################################################

# Print an output of pressure
    maxPressure = max(pressure2)
    s = 0
    for i in range(1,len(pressure2)):
        s = s+pressure2[i]
    averagePressure = old_div(s,len(pressure2))
    val = open('validation.txt', 'w')
    val.write('Only for gauges taken at (x,y)=(0.292,0.04)'+'\n')
    val.write('Maximum pressure [Pa]'+'\t'+'Average pressure [Pa]'+'\n')
    val.write(str(maxPressure)+'\t'+'\t'+'\t'+str(averagePressure))
    val.close()

#####################################################################################

# Print an output file
    info = open('probes.txt','w')
    string1=str(probes)
    string2=string1.replace('[',' ')
    string3=string2.replace(']',' ')   
    string4=string3.replace(',','\t') 
    info.write(str('x')+'\t'+string4+'\n')
    for j in range(1,nRows):
        string5=str(a[j])
        string6=string5.replace('[',' ')
        string7=string6.replace(']',' ')   
        string8=string7.replace(',','\t')   
        info.write(string8+'\n')
    info.close()




