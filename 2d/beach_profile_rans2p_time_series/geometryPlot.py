from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put tankVertices and caissonVertices arrays here

#tankVertices
tv=np.array([[0.0, 0.0],
 [0.38992935246580784, 0.0],
 [1.1697880573974235, 0.0],
 [1.3197880573974234, 0.075],
 [2.0197880573974234, 0.075],
 [2.1697880573974233, 0.0],
 [2.949646762329039, 0.0],
 [3.7295054672606547, 0.0],
 [3.7295054672606547, 1.0],
 [2.949646762329039, 1.0],
 [0.38992935246580784, 1.0],
 [0.0, 1.0],],
                         )
nt=len(tv)

#caissonVertices
cv=np.array( [[ 1.44978806,  0.075     ],
       [ 1.88978806,  0.075     ],
       [ 1.88978806,  0.475     ],
       [ 1.44978806,  0.475     ]],
                        )
nc=len(cv)


xt=[]
yt=[]
xc=[]
yc=[]

for i in range(nt):
    xt.append(tv[i][0])
    yt.append(tv[i][1])
#xt.append(tv[0][0])
#yt.append(tv[0][1])

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





