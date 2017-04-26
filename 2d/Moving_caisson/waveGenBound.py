import numpy as np
from scipy import *
from pylab import *
import collections as cll
import matplotlib.pyplot as plt
import tank as ta

# Put relative path below 
print '################'
print 'Wave generation at the boundary script' 
print '################'
print ''

tStart = 0.0
tEnd = 20.0
dt = 0.01
nT = int( (tEnd - tStart)/dt )
eta = np.zeros(nT)
time = np.linspace(tStart, tEnd, nT, dtype=float)
point = np.array([0.,0.,0.])
for i in range(nT):
    t = time[i]
    eta[i] = ta.waveinput.eta(point,t)

# Plot results
print '################'
print 'Plotting chart of eta function in time'
print '################'
print ''

plt.figure(1)
plt.plot(time, eta, 'k')
plt.xlabel('time [s]')    
plt.ylabel('eta [m]')
plt.tight_layout()
plt.grid(True)
savefig('etaGenBound.png')

