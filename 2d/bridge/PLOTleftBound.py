import numpy as np
from scipy import *
from pylab import *
import collections as cll
import matplotlib.pyplot as plt
import tank as ta

# Put relative path below 
print '################'
print 'Plotting left boundary script' 
print '################'
print ''

t = raw_input('Type the timestep to plot : ')
bottom = 0.0
top = ta.opts.th
dy = ta.he
nE = (top - bottom) / dy
y = np.linspace(bottom,top,nE)
field = []
for i in y:
    field.append(ta.tank.BC['x-'].k_dirichlet.uOfXT((0.,i,0.),t))
field = np.array(field)

# Plot results
print '################'
print 'Plotting chart of eta function in time'
print '################'
print ''

plt.figure(1)
plt.plot(field, y, 'k')
plt.xlabel('field variable [s]')    
plt.ylabel('y [m]')
plt.tight_layout()
plt.grid(True)
savefig('LeftBound.png')

