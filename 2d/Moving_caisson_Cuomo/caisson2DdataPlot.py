import numpy as np
from scipy import *
from pylab import *
import collections as cll
import matplotlib.pyplot as plt

# Put relative path below 
print 'Data plot script' 
filename= 'caisson2D.csv' #raw_input('Type filename.csv : ')       

# Reading file 
with open (filename, 'rb') as csvfile:
    data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
    time=data[:,0]
    data = data[:,1:]
    csvfile.seek(0)
    header = csvfile.readline()
    datalist = [time,data]
print ' Time and data separeted'

# Choosing which data print out
print 'Choose which column plot in time'
index=int(raw_input('Type index [first index is 0 and is NOT time] : ')) 
label=raw_input('Type which label is it and units (like position [m]) :')
x=np.zeros(len(time))
for i in range(len(time)):
    x[i]=data[i][index]

# Plot results 
plt.plot(time,x)
plt.xlabel('time [s]')    
plt.ylabel(label)
plt.show()
savefig('dataPlot_'+label+'.png')


