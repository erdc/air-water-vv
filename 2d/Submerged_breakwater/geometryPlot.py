from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put tankVertices and caissonVertices arrays here

#tankVertices
tv=np.array([[0.0, 0.0],
 [1.393, 0.0],
 [4.179, 0.0],
 [4.679, 0.25],
 [4.929, 0.25],
 [5.304, 0.0],
 [10.1795, 0.0],
 [12.9655, 0.0],
 [12.9655, 0.75],
 [10.1795, 0.75],
 [1.393, 0.75],
 [0.0, 0.75]],
                         )
nt=len(tv)

#caissonVertices
cv=np.array( [ 0.,   0.  ])
nc=len(cv)


xt=[]
yt=[]
xc=[]
yc=[]

for i in range(nt):
    xt.append(tv[i][0])
    yt.append(tv[i][1])
xt.append(tv[0][0])
yt.append(tv[0][1])
"""
for j in range(nc):
    xc.append(cv[j][0])
    yc.append(cv[j][1])
xc.append(cv[0][0])
yc.append(cv[0][1])
"""

# Plot geometry
import matplotlib.pyplot as plt
plt.plot(xt,yt)
#plt.plot(xc,yc)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.suptitle('geometry')
plt.show()
savefig('geometry.png')





