"""
danbreak 2-D
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context,MeshTools)
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges
import math

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",0.5,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.25,"Desired CFL restriction"),
    ("he",0.01,"he relative to Length of domain in x"),
    ("refinement",3,"level of refinement"),
    ("genMesh",False,"generate mesh? "),
    ("usePUMI",False,"usePUMI?"),
    ("adapt", 0, "don't adapt by default"),
    ("mesh",'Reconstructed.smb',"mesh name?"),
    ])

# ****************** #
# ***** GAUGES ***** #
# ****************** #
height_gauges1 = LineGauges(gauges=((("phi",),
                                        (((2.724, 0.0, 0.0),
                                          (2.724, 1.8, 0.0)), # We consider this one in our paper
                                         ((2.228, 0.0, 0.0),
                                          (2.228, 1.8, 0.0)), # We consider this one in our paper
                                         ((1.732, 0.0, 0.0),
                                          (1.732, 1.8, 0.0)),
                                         ((0.582, 0.0, 0.0),
	                                  (0.582, 1.8, 0.0)))),),
                                        fileName="height1.csv")

height_gauges2 = LineGauges(gauges=((("phi",),
                                     (((0.0, 0.0, 0.0),
                                       (0.0, 0.0, -0.01)),
                                      ((0.0, 0.0, 0.0),
                                        (3.22, 0.0, 0.0)))),),
                            fileName="height2.csv")

pressure_gauges = PointGauges(gauges=((('p',),
                                      ((3.22, 0.16, 0.0), #P1
                                       (3.22, 0.584, 0.0), #P3
                                       (3.22, 0.12, 0.0))),), # This is the one considered in our paper
                                       fileName="pressure.csv")


# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
he = 0.015#0.03
tank_dim = (3.22,1.8)
refinement = opts.refinement
structured=False
if structured:
    nny = 5*(2**refinement)+1
    nnx = 2*(nnx-1)+1
    domain = Domain.RectangularDomain(tank_dim)
    boundaryTags = domain.boundaryTags
    triangleFlag=1

elif opts.usePUMI and not opts.genMesh:
    from proteus.MeshAdaptPUMI import MeshAdaptPUMI
    domain = Domain.PUMIDomain(dim=2) #initialize the domain
    #boundaryTags=baseDomain.boundaryFlags
    adaptMeshFlag = opts.adapt#1
    adaptMesh_nSteps = 50
    adaptMesh_numIter = 10
    hmax = he*2.0;
    hmin = he/2.0;
    hPhi = he/2.0;#/4.0
    domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI(hmax=hmax, hmin=hmin, hPhi = hPhi, adaptMesh=adaptMeshFlag, numIter=adaptMesh_numIter, numAdaptSteps=adaptMesh_nSteps,  sfConfig="combined",logType="off",reconstructedFlag=2,gradingFact=1.2,targetError=4.0)
    #read the geometry and mesh
    parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
    domain.MeshOptions.setParallelPartitioningType('element')
    domain.PUMIMesh.loadModelAndMesh("Reconstructed.dmg", opts.mesh)
    nnx = nny = None

else:
    nnx = nny = None
    domain = Domain.PlanarStraightLineGraphDomain()

if opts.genMesh and opts.usePUMI:
  from proteus.MeshAdaptPUMI import MeshAdaptPUMI
  domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI()

# ----- TANK ----- #
tank = Tank2D(domain, tank_dim)

# ----- EXTRA BOUNDARY CONDITIONS ----- #
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

#he = 0.06#tank_dim[0]*opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
if opts.genMesh and not opts.usePUMI:
  domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.0

waterLine_y = 0.6
waterLine_x = 1.2
class vof_init_cond(object):
    def uOfXT(self,x,t):
        if x[0] < waterLine_x and x[1] < waterLine_y:
            return -1.0
        elif x[0] > waterLine_x or x[1] > waterLine_y:
            return 1.0
        else:
            return 0.0

ecH = 3.0
def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_y = x[1] - waterLine_y
    if phi_x < 0.0:
        if phi_y < 0.0:
            return max(phi_x, phi_y)
        else:
            return phi_y
    else:
        if phi_y < 0.0:
            return phi_x
        else:
            return math.sqrt(phi_x ** 2 + phi_y ** 2)

class PerturbedSurface_phi:
    def uOfXT(self,x,t):
        return signedDistance(x)

from proteus.ctransportCoefficients import smoothedHeaviside
class PerturbedSurface_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(ecH*he,signedDistance(x))

initialConditions  = {0:PerturbedSurface_H()}

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
outputStepping.systemStepExact = True
initialConditions = {'pressure': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vof': PerturbedSurface_H(),
                     'ncls': PerturbedSurface_phi(),
                     'rdls': PerturbedSurface_phi()}
boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
    'vel_u_DBC': lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
    'vel_v_DBC': lambda x, flag: domain.bc[flag].v_dirichlet.init_cython(),
    'vel_w_DBC': lambda x, flag: domain.bc[flag].w_dirichlet.init_cython(),
    'vof_DBC': lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython(),
    'ncls_DBC': lambda x, flag: domain.bc[flag].ncls_dirichlet.init_cython(),
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': lambda x, flag: domain.bc[flag].p_advective.init_cython(),
    'vel_u_AFBC': lambda x, flag: domain.bc[flag].u_advective.init_cython(),
    'vel_v_AFBC': lambda x, flag: domain.bc[flag].v_advective.init_cython(),
    'vel_w_AFBC': lambda x, flag: domain.bc[flag].w_advective.init_cython(),
    'vof_AFBC': lambda x, flag: domain.bc[flag].vof_advective.init_cython(),
    'ncls_AFBC': lambda x, flag: domain.bc[flag].ncls_advective.init_cython(),
    # DIFFUSIVE FLUX BCs #
    'vel_u_DFBC': lambda x, flag: domain.bc[flag].u_diffusive.init_cython(),
    'vel_v_DFBC': lambda x, flag: domain.bc[flag].v_diffusive.init_cython(),
    'vel_w_DFBC': lambda x, flag: domain.bc[flag].w_diffusive.init_cython(),
    'vof_DFBC': lambda x, flag: None,
    'ncls_DFBC': lambda x, flag: None}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=structured,
                                             he=he,
                                             nnx=nnx,
                                             nny=nny,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=False)
myTpFlowProblem.Parameters.physical['gravity'] = np.array([0.0,-9.8,0.0])
myTpFlowProblem.useBoundaryConditionsModule = False
m = myTpFlowProblem.Parameters.Models
m.clsvof.p.CoefficientsOptions.disc_ICs = True
m.clsvof.auxiliaryVariables = [height_gauges1, height_gauges2]
m.pressure.auxiliaryVariables = [pressure_gauges]
m.rans2p.n.ShockCapturingOptions.shockCapturingFactor = 0.5

myTpFlowProblem.Parameters.mesh.he = he
myTpFlowProblem.Parameters.mesh.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)
