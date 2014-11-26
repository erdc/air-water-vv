###############INPUT
L = (1.6 ,0.61, 0.75)
obst_portions = (0.12,0.12,0.75) #(width_x,width_y,height)
obst_x_start = 0.9 # start x coordinate of the obstacle; caution to be in the domain's range
obst_x_end = obst_x_start + obst_portions[0] # end x coordinate of the obstacle; caution to be in the domain's range
obst_y_start =0.25  # start y coordinate of the obstacle; caution to be in the domain's range
obst_y_end=obst_y_start+obst_portions[1] # end y coordinate of the obstacle; caution to be in the domain's range
obst = (obst_x_start,obst_x_end,obst_y_start,obst_y_end,obst_portions[2]) #coordinates of the obstacle to be used to define the boundary

vertices=[[0.0,0.0,0.0], #0
                  [obst[0],0.0,0.0], #1
                  [obst[0],0.0,L[2]], #2
                  [0.0,0.0,L[2]], #3
                  #plane 2
                  [0.0,obst[2],0.0], #4
                  [obst[0],obst[2],0.0], #5
                  [obst[0],obst[2],L[2]], #6
                  [0.0,obst[2],L[2]], #7
                  #plane 3                  
                  [0.0,obst[3],0.0], #8
                  [obst[0],obst[3],0.0], #9
                  [obst[0],obst[3],L[2]], #10
                  [0.0,obst[3],L[2]], #11
                  #plane 4
                  [0.0,L[1],0.0], #12
                  [obst[0],L[1],0.0], #13
                  [obst[0],L[1],L[2]], #14
                  [0.0,L[1],L[2]], #15
                  #plane 1 extension 1
                  [obst[1],0.0,0.0], #16
                  [obst[1],0.0,L[2]], #17
                  #plane 2 extension 1
                  [obst[1],obst[2],0.0], #18
                  [obst[1],obst[2],L[2]], #19
                  #plane 3 extension 1
                  [obst[1],obst[3],0.0], #20
                  [obst[1],obst[3],L[2]], #21
                  #plane 4 extension 1                 
                  [obst[1],L[1],0.0], #22
                  [obst[1],L[1],L[2]], #23
                  #plane 1 extension 2
                  [L[0],0.0,0.0], #24
                  [L[0],0.0,L[2]], #25
                  #plane 2 extension 2
                  [L[0],obst[2],0.0], #26
                  [L[0],obst[2],L[2]], #27
                  #plane 3 extension 2
                  [L[0],obst[3],0.0], #28
                  [L[0],obst[3],L[2]], #29
                  #plane 4 extension 2
                  [L[0],L[1],0.0], #30
                  [L[0],L[1],L[2]] #31
                  ]                


     #    for v,vF in zip(vertices,vertexFlags):
             
     #        vertices.append([v[0],L[1],v[2]])
     #        vertexFlags.append(vF)
     

segments=[[0,1],
                  [1,16],
                  [16,24],
                  [24,25],
                  [25,17],
                  [17,2],
                  [2,3],
                  [3,0],
                  [12,13],
                  [13,22],
                  [22,30],
                  [30,31],
                  [31,23],
                  [23,14],
                  [14,15],
                  [15,12],
                  [0,4],
                  [4,8],
                  [8,12],
                  [24,26],
                  [26,28],
                  [28,30],
                  [25,27],
                  [27,29],
                  [29,31],
                  [3,7],
                  [7,11],
                  [11,15],
                  [5,18],
                  [18,19],
                  [19,6],
                  [6,5],
                  [9,20],
                  [20,21],
                  [21,10],
                  [10,9],
                  [5,9],
                  [18,20],
                  [19,21],
                  [6,10]]
  


        


facets=[[[0,4,8,12,15,11,7,3]],#left
                [[24,26,28,30,31,29,27,25]],#right
                [[5,18,19,6]],  #front facet of obstacle
                [[9,20,21,10]],#back facet of obstacle
                [[5,9,10,6]], #left facet of obstacle
                [[18,20,21,19]], #right facet of obstacle
                [[0,1,5,4]],#bottom
                [[4,5,9,8]],#bottom
                [[8,9,13,12]],#bottom
                [[1,16,18,5]],#bottom
                [[18,16,24,26]],#bottom
                [[18,26,28,20]],#bottom
                [[20,28,30,22]],#bottom
                [[20,22,13,9]],#bottom
                [[3,2,6,7]],  #top
                [[7,6,10,11]],#top
                [[11,10,14,15]], #top
                [[2,17,19,6]],#top
                [[19,17,25,27]],#top
                [[19,27,29,21]],#top
                [[21,29,31,23]],#top
                [[10,21,23,14]],#top
                [[0,1,16,24,25,17,2,3]], #front
                [[12,13,22,30,31,23,14,15]]    #back        
                ]
        

#=================================
# We can make it a function and directly load these from proteus by including the following line, or alternatively, by importing the namecase.py file
#def 3d_geom(vertices,segments,facets):
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

# Plotting facets, currently only valid for facets aligned with the main axes
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

