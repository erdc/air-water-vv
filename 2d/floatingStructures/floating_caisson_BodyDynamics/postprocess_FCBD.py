import numpy as np
import matplotlib.pyplot as plt
import math
import AnalysisTools as at

## Reading probes into the file
probes = 'caisson2D.csv'        
        
datalist = at.readProbeFile(probes)
probeType = datalist[0]
time = datalist[1]
data = datalist[2]
rz = data[:,5]

alpha = []
for i in range(0,len(rz)):
    alpha.append(rz[i]*180/math.pi)
alpha = np.array(alpha)

it = np.where(time>2.5)[0][0]
period = at.zeroCrossing(time[:it],alpha[:it],up=False)[0]

period_ref = 0.93
err = 100*abs(period_ref-period)/abs(period_ref)
val = open('validation_FCBD.txt', 'w')
val.write('Period for the rotation angle'+'\n')
val.write('Theory'+'\t'+'Simulation'+'\t'+'\t'+'Error (%)'+'\n')
val.write(str(period_ref)+'\t'+str(period)+'\t'+str(err))
val.close()

plt.plot(time,alpha)
plt.xlabel('time [sec]')
plt.ylabel('rotation angle [deg]')
plt.suptitle('Rotation angle for the floating caisson 2D (body dynamics)')
plt.xlim((0.,4.5))
plt.grid(True)
plt.savefig('rotation_angle.png')
plt.show()
