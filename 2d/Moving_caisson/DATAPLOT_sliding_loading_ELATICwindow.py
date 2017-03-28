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


#   Start time of analysis
stime = raw_input('Enter start time for analysis (press Enter for default 0): ')
try:
    stime = float(stime)    
    ja = min(find(time>=stime))
except:
    stime = 0
    ja = 0
print('Start time entered: '+str(stime))
#   end time of analysis. Preferable to put integer number of periods
etime = raw_input('Enter end time for analysis (preferably an integer number of periods, press Enter for default end time): ')
try:
    etime = float(etime)    
    jb = min(find(time>=etime))
except:
    etime = time[-1]
    jb = len(time)-1
print('End time entered: '+str(etime))


# Choosing which data print out
xPos, yPos, zRot, Fx, Fy, Mz, uxPl, uxEl, uy = [], [], [], [], [], [], [], [], [] 
times = []
for j in range(len(time)):
    if time[j]>stime and time[j]<etime:
        times.append(time[j]-stime)
        xPos.append(data[j][0])
        yPos.append(data[j][1])
        zRot.append(data[j][5])
        Fx.append(data[j][6])
        Fy.append(data[j][7])
        Mz.append(data[j][11])
        uxEl.append(data[j][-3])
        uxPl.append(data[j][-2])
        uy.append(data[j][-1])

# Post-processing
time=times
xInit=xPos[0]
xPos1=(xPos-xInit)*1000
yInit=yPos[0]
yPos1=(yPos-yInit)*1000

mu=0.5
width=0.4
Kx=541553.2
Ky=582633.7

Fx1, Fy1, Fy2, xTrasl, xPl = [], [], [], [], []
for iii in range(len(time)):
    Fx1.append(Fx[iii]*width)
    Fy1.append(abs(mu*width*(Ky*uy[iii])))
    Fy2.append((Kx*uxEl[iii])*width)
    xTrasl.append((uxEl[iii]+uxPl[iii]-uxPl[0])*1000.)
    xPl.append((uxPl[iii]-uxPl[0])*1000)

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



plt.figure(7)
lineXT, = plt.plot(time,xTrasl, 'k', label='Total x-axis displacement')
lineXP, = plt.plot(time,xPl, 'b', label='Plastic x-axis displacement')
plt.xlabel('time [s]')    
plt.ylabel('X-axis displacements [mm]')
plt.tight_layout()
plt.grid(True)
plt.legend(handles=[lineXT, lineXP], ncol=1, loc=2)
savefig('dataPlotxT-xP.png')



