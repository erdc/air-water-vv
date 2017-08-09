from proteus import Domain
from proteus.mprans import SpatialTools as st
import numpy as np
from nbtools import plot_domain, plot_js


# plot_domain(domain3D)

import pandas as pd
import matplotlib.pyplot as plt

scale = 1.0 / 21.0

caisson_width = 10.5625*25*0.0254  # 1:25 (in) -> 1:1 (in) -> 1:1 (mts)



domain3D = Domain.PiecewiseLinearComplexDomain()
# myTank=st.Tank3D(domain3D,dim=(0.9144,3.0,1.0,))
he_min = 0.02


#caisson_coord=[0.9144/2, 1.5, .5]
caisson_coord=[0.5,0.5,0.5]


xl = pd.ExcelFile("Ribbon Bridge Section Line Segment.xlsx")
df = xl.parse(0)
y = np.asarray(df['x'].tolist())
z = np.asarray(df['y'].tolist())
# Offset to (0,0)
y -= 0.5 * (max(y) + min(y))
z -= min(z)

# Real dimension
prescale = 25 * 0.0254  # Scale of 25 * inches to meters
y = y * prescale
z = z * prescale
# Scale to experimental model
yp = y * scale
zp = z * scale
caisson_width *= scale

dim = [caisson_width, max(yp) - min(yp), max(zp) - min(zp)]

# Interpolate points to he_caisson
yd = np.diff(yp)
zd = np.diff(zp)
dist = np.sqrt(yd ** 2 + zd ** 2)
u = np.cumsum(dist)
u = np.hstack([[0], u])

t = np.linspace(0, u.max(), int(u.max() / he_min))
yBase = np.interp(t, u, yp)
zBase = np.interp(t, u, zp)

yTop = np.linspace(zp.min() + he_min, zp.max(), int((zp.max() - zp.min()) / he_min), endpoint=False)
zTop = np.full(len(yTop), zp[-1])
y1 = np.hstack([yBase, yTop])
z1 = np.hstack([zBase, zTop])

n1 = len(y1)

# Extrude for 3d object

x1 = np.empty(n1)
x1.fill(caisson_width / 2)

x2 = -x1
y2 = y1
z2 = z1

x = np.hstack([x1, x2])
y = np.hstack([y1, y2])
z = np.hstack([z1, z2])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# ax.scatter(x, y, z, c='b', marker='o')
#
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
#
# plt.show()

barycenter = (0, 0, 0.)

# Calculate Parameters
VCG=None
if VCG is None:
    VCG = dim[1]/2.
#I = ((dim[0]**2)+(dim[1]**2))*mass/12.
barycenter = (0, 0, 0.)


#Right Center
vertices1 = []
vertexFlags1 = []
segments1 = []
segmentFlags1 = []

v_start = 0
for i in range(len(x1)):
    v_start = len(vertices1)
    v = [[x1[i],y1[i],z1[i]]]
    if v_start >= 1:
        s = [[v_start-1, v_start]]
    else:
        s = []
    vertices1 += v
    vertexFlags1 += [1]*len(v)
    segments1 += s
    segmentFlags1 += [1]*len(s)
segments1[-1][1] = 0  # last segment links to vertex 0

facet1 = []
for i, vert in enumerate(vertices1):
    facet1 += [i]

#Left Side

vertices2 = []
vertexFlags2 = []
segments2 = []
segmentFlags2 = []
n1=len(vertices1)-1

v_start = 0
for i in range(len(x2)):
    v_start = len(vertices2)
    v = [[x2[i],y2[i],z2[i]]]
    if v_start >= 1:
        s = [[n1+v_start-1, n1+v_start]]
    else:
        s = []
    vertices2 += v
    vertexFlags2 += [1]*len(v)
    segments2 += s
    segmentFlags2 += [1]*len(s)

segments2[-1][1] = n1   # last segment links to vertex 0

facet2 = []
for i, vert in enumerate(vertices2):
    facet2 += [n1+i+1]

#Intermediate Segments

segments3=[]
segmentFlags3=[]
for i in range(n1+1):
    s = [[i , n1+i+1]]
    segments3 += s
    segmentFlags3 += [1]*len(s)


# Add Intermediate facets

facet3=[]
for i in range(n1):
    facet3 += [[[i , i+1, n1+i+2, n1+i+1]]]
facet3 += [[[i+1 , 0, n1+1, n1+i+2]]]



vertices=vertices1+vertices2
segments=segments1+segments2+segments3
segmentFlags=segmentFlags1+segmentFlags2
vertexFlags=vertexFlags1+vertexFlags2
facets=[[facet1],[facet2]]+facet3
facetFlags = np.array([1]*len(facets))
regionFlags = np.array([1])

boundaryTags = {'caisson': 1}

caisson = st.CustomShape(domain3D, barycenter=barycenter,
                         vertices=vertices,
                         vertexFlags=vertexFlags,
                         segments=segments,
                         segmentFlags=segmentFlags,
                         facets=facets,
                         facetFlags=facetFlags,
                         boundaryTags=boundaryTags
                         )

caisson.setHoles([[0.,0., dim[2] / 4]])
caisson.holes_ind = np.array([0])
caisson.translate(caisson_coord)

# myTank.setChildShape(caisson, 0)

plot_domain(domain3D)
plt.show()



