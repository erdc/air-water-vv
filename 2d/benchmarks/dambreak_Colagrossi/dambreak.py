"""
Dambreak flow - Collagrosi and Landrini (2003)
"""
import numpy as np
from math import sqrt
from proteus import (Domain, Context,
                     FemTools as ft,
                     #SpatialTools as st,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow


# predefined options
opts=Context.Options([
    # water column 
    ("mwl", 0.6, "Height of water column in m"),
    ("water_width", 1.2, "Width of  water column in m"),
    ("rho_0",998.2,"water density"),
    ("nu_0",1.004e-6,"water viscosity"),
    ("rho_1",1.205,"air density"),
    ("nu_1",1.500e-5,"air viscosity"),
    ("sigma",0.0, "surface tension"),
    # tank
    ("tank_dim", (3.22, 1.8), "Dimensions of the tank  in m"),
    #gravity 
    ("g",(0,-9.81,0), "Gravity vector in m/s^2"),
    # gauges
    ("gauge_output", True, "Produce gauge data"),
    ("gauge_location_p", (3.22, 0.12, 0), "Pressure gauge location in m"),
    # mesh refinement and timestep
    ("refinement", 32 ,"Refinement level, he = L/(4*refinement - 1), where L is the horizontal dimension"), 
    ("cfl", 0.33 ,"Target cfl"),
    ("ecH", 3,"Smoothing Coefficient"),
    # run time options
    ("T", 0.09 ,"Simulation time in s"),
    ("dt_fixed", 0.01, "Fixed time step in s"),
    ("dt_init", 0.001 ,"Maximum initial time step in s"),
    ("useHex", False, "Use a hexahedral structured mesh"),
    ("structured", False, "Use a structured triangular mesh"),
    ("gen_mesh", True ,"Generate new mesh"),
    ])



# CONTEXT

# water
waterLine_z = opts.mwl
waterLine_x = opts.water_width

# tank
tank_dim = opts.tank_dim

#     Discretization Input Options       


#[temp] temporary location
backgroundDiffusionFactor = 0.01

movingDomain = False
checkMass = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = 1
useHex = opts.useHex
structured = opts.structured
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
"""

# ----- INPUT CHECKS ----- #
if spaceOrder not in [1,2]:
    raise ValueError("INVALID: spaceOrder(" + str(spaceOrder) + ")")

if useRBLES not in [0.0, 1.0]:
    raise ValueError("INVALID: useRBLES(" + str(useRBLES) + ")")

if useMetrics not in [0.0, 1.0]:
    raise ValueError("INVALID: useMetrics(" + str(useMetrics) + ")")

# ----- DISCRETIZATION ----- #
nd = 2
if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = ft.C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = ft.C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, 4)

"""

# Numerical Options and Other Parameters             
weak_bc_penalty_constant = 100.0
nLevels = 1


# ----- TIME STEPPING & VELOCITY----- #

dt_init = min(0.1 * opts.dt_fixed, opts.dt_init)

nDTout = int(round(opts.T / opts.dt_fixed))

# ----- DOMAIN ----- #

# domain = Domain.PlanarStraightLineGraphDomain()
# domain replacement
if useHex:
    nnx = 4 * opts.refinement + 1
    nny = 2 * opts.refinement + 1
    hex = True
    domain = Domain.RectangularDomain(tank_dim)
elif structured:
    nnx = 4 * opts.refinement
    nny = 2 * opts.refinement
    domain = Domain.RectangularDomain(tank_dim)
    boundaryTags = domain.boundaryTags
else:
    domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #

tank = Tank2D(domain, tank_dim)

# ----- GAUGES ----- #

if opts.gauge_output:
    tank.attachPointGauges(
        'twp',
        gauges = ((('p',), (opts.gauge_location_p,)),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pressureGauge.csv'
    )


# ----- EXTRA BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

# ----- MESH CONSTRUCTION ----- #

he = tank_dim[0] / float(4 * opts.refinement - 1)
domain.MeshOptions.he = he
st.assembleDomain(domain)

# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# Numerics
def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x, phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x ** 2 + phi_z ** 2)


class P_IC:
    def __init__(self):
        self.mwl=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1] - (opts.mwl - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1]
class AtRest:
    def uOfXT(self, x, t):
        return 0.0

class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)
Duration= 10.
dt_output = 0.05
initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest()}

initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()


outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                       dt_init=opts.dt_init,
                                       # cfl=cfl,
                                       dt_output=dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None, # set with SpatialTools,
                                             )

params = myTpFlowProblem.Parameters

myTpFlowProblem.useSuperLu=False#True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma

# index in order of
m = params.Models
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4
# ----- TURBULENCE MODELS ----- #
useRANS=0

ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4

##########################################
#            Initial conditions for free-surface                     #
##########################################



#Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
