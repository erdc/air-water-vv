from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put relative path below
filename='point_gauges.csv'

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
# Choose which time step to plot  
    print('Number of rows : '+ str(nRows))
    line = int(raw_input('Enter which line (time step) to plot : '))
    pressure=a[line-1][1:]  
 # Plot pressure in space
    import matplotlib.pyplot as plt
    plt.plot(probes,pressure)
    plt.xlabel('probes location [m]')    
    plt.ylabel('pressure[Pa]')
    plt.show()
    savefig('pressure_in_space_%d.png' %(line))

# Choose which probes to plot  
    print('Number of probes : '+ str(len(probes)))
    x = int(raw_input('Enter which probes to plot : '))
    pressure2=a[x-1]  
 # Plot pressure in time
    import matplotlib.pyplot as plt
    plt.plot(time,pressure2)
    plt.xlabel('time step [sec]')    
    plt.ylabel('pressure [Pa]')
    plt.show()
    savefig('pressure_in_space_%d.png' %(x))
