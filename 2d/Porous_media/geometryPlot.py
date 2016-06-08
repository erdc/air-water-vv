from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put tankVertices and caissonVertices arrays here

#tankVertices
tv=np.array([[0.0, 0.0],
 [8.9830376539957904, 0.0],
 [17.966075307991581, 0.0],
 [18.932075307991582, 0.644],
 [20.072075307991582, 0.644],
 [20.072075307991582, 0.244],
 [21.072075307991582, 0.244],
 [21.560075307991582, 0.244],
 [22.048075307991581, 0.0],
 [31.031112961987372, 0.0],
 [40.014150615983162, 0.0],
 [40.014150615983162, 2.0],
 [31.031112961987372, 2.0],
 [8.9830376539957904, 2.0],
 [0.0, 2.0]],
                         )
nt=len(tv)

#caissonVertices
cv=np.array( [[ 20.07207531,   0.244     ],
       [ 21.07207531,   0.244     ],
       [ 21.07207531,   1.364     ],
       [ 20.07207531,   1.364     ]],
                        )
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

for j in range(nc):
    xc.append(cv[j][0])
    yc.append(cv[j][1])
xc.append(cv[0][0])
yc.append(cv[0][1])


# Plot geometry
import matplotlib.pyplot as plt
plt.plot(xt,yt)
plt.plot(xc,yc)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.suptitle('geometry')
plt.show()
savefig('geometry.png')





