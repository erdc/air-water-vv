import numpy as np
from scipy import *
from pylab import *
import collections as cll
import matplotlib.pyplot as plt
#import tank as tank

# Put relative path below 
print '################'
print 'Data plot script' 
print '################'
print ''

filename= 'caisson2D.csv' #raw_input('Type filename.csv : ')       

# Reading file 
with open (filename, 'rb') as csvfile:
    data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
    time=data[:,0]
    data = data[:,1:]
    csvfile.seek(0)
    header = csvfile.readline()
    datalist = [time,data]

print '################'
print 'Time and data separeted'
print '################'
print ''

# Choosing which data print out
xPos, yPos, zRot, Fx, Fy, Mz, uxPl, uxEl, uy = np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time))
for j in range(len(time)):
    xPos[j]=data[j][0]
    yPos[j]=data[j][1]
    zRot[j]=data[j][5]
    Fx[j]=data[j][6]
    Fy[j]=data[j][7]
    Mz[j]=data[j][11]
    uxEl[j]=data[j][-3]
    uxPl[j]=data[j][-2]
    uy[j]=data[j][-1]

# Post-processing
xInit=xPos[0]
xPos1=(xPos-xInit)*1000
yInit=yPos[0]
yPos1=(yPos-yInit)*1000

mu=0.5
width=0.4
Kx=541553.2
Ky=582633.7

Fx1=Fx*width
Fy1=abs(mu*width*(Ky*uy))
Fy2=(Kx*uxEl)*width

# Plot results
print '################'
print 'Plotting charts'
print '################'
print ''


plt.figure(1)
plt.plot(time,xPos1, 'k')
plt.xlabel('time [s]')    
plt.ylabel('Position x-axis [mm]')
plt.tight_layout()
plt.grid(True)
savefig('dataPlotX.png')

plt.figure(2)
plt.plot(time,yPos1, 'k')
plt.xlabel('time [s]')    
plt.ylabel('Position y-axis [mm]')
plt.tight_layout()
plt.grid(True)
savefig('dataPlotY.png')

plt.figure(3)
plt.plot(time,zRot, 'k')
plt.xlabel('time [s]')    
plt.ylabel('Rotation z-axis [rad]')
plt.tight_layout()
plt.grid(True)
savefig('dataPlotZrot.png')

plt.figure(4)
lineFx, =plt.plot(time,Fx1, 'b', label='Horizontal load')
lineFy, =plt.plot(time,Fy1, 'r', label='Frictional force')
plt.xlabel('time [s]')    
plt.ylabel('Loadings [N]')
plt.tight_layout()
plt.grid(True)
plt.legend(handles=[lineFx, lineFy], ncol=1, loc=2)
savefig('dataPlotLoadings.png')

plt.figure(5)
lineFy, =plt.plot(time,Fy1, 'r', label='Frictional force')
lineFy2, =plt.plot(time,Fy2, 'g', label='Soil reaction')
plt.xlabel('time [s]')    
plt.ylabel('Reactions [N]')
plt.tight_layout()
plt.grid(True)
plt.legend(handles=[lineFy, lineFy2], ncol=1, loc=2)
savefig('dataPlotReactions(elastic).png')



plt.figure(6)
plt.plot(time,Mz, 'k')
plt.xlabel('time [s]')    
plt.ylabel('Moment z-axis [N.m]')
plt.tight_layout()
plt.grid(True)
savefig('dataPlotMomentZ.png')



