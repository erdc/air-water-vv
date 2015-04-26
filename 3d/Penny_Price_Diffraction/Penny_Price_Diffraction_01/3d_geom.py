###############INPUT
import tank3D
#=================================

# Vertex coordinates
vx=[]
vy=[]
vz=[]          
for coord in vertices:
    vx.append(coord[0])
    vy.append(coord[1])
    vz.append(coord[2])
# Define which faces to show: FaceX, FaceY, FaceZ, allFaces
setface = "FaceX"

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


#def zz(x,y):
    
# Figure
fig = plt.figure()
ax = fig.gca(projection='3d')
#Plotting vertices
sc1=ax.scatter(vx, vy, vz)
#Plotting line segments
for seg in segments:
    l1=ax.plot([vx[seg[0]],vx[seg[1]]],[vy[seg[0]],vy[seg[1]]],[vz[seg[0]],vz[seg[1]]])

# Plotting facets
for face in facets:
    #Coordinates of face vertices
    fx = []
    fy = []
    fz = []
    #Counter for determining if it is a X,Y,Z facet  --- 
    icount = np.array([0.,0.,0.])
    jj=-1
    for vertex in face[0]:
        jj+=1
        fx.append(vx[vertex])
        fy.append(vy[vertex])
        fz.append(vz[vertex])
        icount[0]=icount[0] + abs(fx[jj]-fx[max(0,jj-1)])
        icount[1]=icount[1] + abs(fy[jj]-fy[max(0,jj-1)])
        #print abs(fy[jj]-fy[max(0,jj-1)])
        icount[2]+=abs(fz[jj]-fz[max(0,jj-1)])
        #print abs(fz[jj]-fz[max(0,jj-1)])
    #fx.append(vx[face[0]]+1e-30)
    #fy.append(vy[face[0]]+1e-30)
    #fz.append(vz[face[0]]+1e-30)
    print icount
    if icount[0]<=1e-30 and ((setface is "FaceX") or (setface is "allFaces")):
        print "FaceX"
        Y,Z = np.meshgrid(fy,fz)
        Z1,X = np.meshgrid(fz,fx)
        s1=ax.plot_surface(X,Y,Z,rstride=1, cstride=1,color="r",alpha=0.1)

    elif icount[1]<=1e-30 and ((setface is "FaceY") or (setface is "allFaces")):
        print "FaceY"
        X,Z = np.meshgrid(fx,fz)
        Z1,Y = np.meshgrid(fz,fy)
        s2=ax.plot_surface(X,Y,Z,rstride=1, cstride=1,color="r",alpha=0.1)
        
    elif icount[2]<=1e-30 and ((setface is "FaceZ") or (setface is "allFaces")):
        print "FaceZ"
        X,Y = np.meshgrid(fx,fy)
        Y1,Z = np.meshgrid(fy,fz)
        s3=ax.plot_surface(X,Y,Z,rstride=1, cstride=1,color="r",alpha=0.1)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()

