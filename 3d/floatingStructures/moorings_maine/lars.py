from pythreejs import *
from IPython.display import display
import sys, os, random, time
import numpy as np
from math import *


x = 0
y = 0
width = 640
height = 480

spheres =[]

from proteus.mbd import ChRigidBody as crb

g = np.array([0., 0., -9.807])
system = crb.ProtChSystem(g)
system.setTimeStep(0.001)
timestepper = "Euler"
if timestepper == "HHT":
    system.setTimestepperType(timestepper, verbose=False)
mesh = crb.Mesh(system)

length = np.array([5.672,0.126,4.,0.259])
nb_elems = np.array([40,1,40,2])
d = np.array([0.0049,0.006, 0.009, 0.008])
A = d**2/4.*np.pi
rho = np.array([0.402,1.558,0.00425,1.529])/A
E = np.array([2.0505e6,3.63e6,10.873e3,6.464e6])/A

length = np.array([7.00])
nb_elems = np.array([100])
d = np.array([2.5e-3])
A = d**2/4.*np.pi
rho = np.array([1.036/9.807])/A
w = np.array([1.036])
EA = np.array([560e3])
E = np.array([560e3])/A


anchor1 = [0.,0., 0.]
fairlead1 = [6.339, 0., 2.65]

import catenary as cat
l1 = cat.MooringLine(L=length, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1, tol=1e-8)
l1.setVariables()

func = lambda s: (s, s, s)

body = crb.ProtChBody(system)
body.setConstraints(free_x=np.array([0.,0.,0.]), free_r=np.array([0.,0.,0.]))
body.ChBody.SetMass(50.)
#body.setPosition(pos[-1])

moorings = crb.ProtChMoorings(system, mesh, length, nb_elems, d, rho, E, "CableANCF")
external_forces_bool = True
moorings.external_forces_manual = external_forces_bool
Cd = False
if external_forces_bool is True:
    moorings.setFluidDensityAtNodes(np.array([1000 for i in range(nb_elems+1)]))
    if Cd is True:
        moorings.setDragCoefficients(0.5, 2.5, 0)
        moorings.setAddedMassCoefficients(0., 3.8, 0)
moorings.setNodesPositionFunction(l1.s)
moorings.setNodesPosition()
moorings.buildNodes()

pos = moorings.getNodesPosition()



from proteus.mbd import pyChronoCore as cc
material = cc.ChMaterialSurfaceSMC()
material.SetYoungModulus(2e4)
material.SetFriction(0.3)
material.SetRestitution(0.2)
material.SetAdhesion(0)

moorings.attachBackNodeToBody(body)
moorings.fixFrontNode(True)
#moorings.fixBackNode(True)


moorings.setContactMaterial(material)


box_pos = np.array([0.,0.,-0.11])
box_dim = np.array([20.,1.,0.2])
vec = cc.ChVector(0., 0., -0.11)
#box = cc.ChBodyEasyBox(system, box_dim[0], box_dim[1], box_dim[2], 1000, True)
#box.SetPos(vec)
#box.SetMaterialSurface(material)
#box.SetBodyFixed(True)

# Some variables used inside the simulation loop
fps = 100
dt = 1.0/fps
running = True
state = 0
counter = 0
objcount = 0
lasttime = time.time()

f = 0.79/(2*np.pi)
amplitude = 0.1


moorings.setName('test')
system.calculate_init()

body.SetBodyFixed(True)
#body.setPrescribedMotionSine(a=amplitude, f=f)

T = 2*np.pi/w
t0_end = T
dt0 = 0.01
system.setTimeStep(1e-3)
print("starting_calculation")
#for i in range(int(t0_end/dt0)):
#    system.calculate(0.01)
#    Tf = moorings.getTensionBack()
#    print('T0', system.GetChTime(), Tf, np.linalg.norm(Tf))

Tf = moorings.getTensionBack()
print(system.GetChTime(), Tf, np.linalg.norm(Tf))

#dt = 0.01
#for i in range(20):
#    system.setTimeStep(1e-5)
#    system.calculate(dt)

system.setTimeStep(1e-3)
for i in range(int(10*T/dt)):
    system.calculate(dt)
    Tf = moorings.getTensionBack()
    print(system.GetChTime(), Tf, np.linalg.norm(Tf))
    

