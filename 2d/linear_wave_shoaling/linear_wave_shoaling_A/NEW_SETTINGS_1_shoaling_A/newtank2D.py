from proteus import Domain
from proteus import SpatialTools as st
import ode
from math import *
import numpy as np


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank_dim = [30.0, 1.50]
slope=1.0/10.0
x_slope_start=5.0
slope_length=6.0
x_slope_end=x_slope_start + slope_length
y_slope_end=slope*slope_length
leftSponge = wavelength
rightSponge= wavelength*2.0
b_or = {'bottom': [0., -1.],
        'right': [1., 0.],
        'top': [0., 1.],
        'left': [-1., 0.],}
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,}
vertices=[[0.0,0.0],#0
          [x_slope_start,0.0],#1
          [x_slope_end,y_slope_end],#2
          [tank_dim[0],y_slope_end],#3
          [tank_dim[0],tank_dim[1]],#4
          [0.0,tank_dim[1]]]#5
vertexFlags=np.array([1, 1, 1, 1, 3, 3])
segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,0]]
segmentFlags=np.array([1, 1, 1, 2, 3, 4])
regions=[ [ 0.1*tank_dim[0] , 0.1*tank_dim[1] ],
          [0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]
regionFlags=np.array([1,2])
tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=b_or)


from proteus import BC
from proteus import WaveTools
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

# ----- BOUNDARY CONDITIONS ----- #

#They are bit different from the original tank with old BC
#Here I suppose there is not the smoothing function for velocity at boundary
H=smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
U=Vcomp*(1.0-H)
tank.BC.top.setOpenAir()
tank.BC.left.setTwoPhasesVelocityInlet(U=U, waterLevel=waterLevel, vert_axis=-1, air=1., water=0.)
tank.BC.bottom.setNoSlip()
tank.BC.right.reset()


# ----- WAVE CONDITIONS ----- #
period=1.94
waterLevel = 1
wavelength = 5.0
waveHeight = 0.025
waveinput=MonochromaticWaves(period=period,waterLevel=depth,g=g,waveType="Linear",wavelenght=wavelenght,waveHeight=waveHeight)

left = right = False
if leftSponge is not None: left = True
if rightSponge is not None: right = True
tank.setAbsorptionZones(left=left, right=left)

from proteus import Gauges

# ----- Output Gauges ----- #
line_output=LineGauges(gauges=((('u', 'v'), (((2.50, 1.50), (2.50, 0.)),
                                               ((5.50, 1.50), (5.50, slope*5.50)),
                                               ((6.50, 1.50), (6.50, y_slop_end)))),
                                 (('p'), (((2.50, 1.50), (2.50, 0.)),
                                               ((5.50, 1.50), (5.50, slope*5.50)),
                                               ((6.50, 1.50), (6.50, y_slop_end)))),
                                 ),
                         activeTime=(0, T),
                         sampleRate=1/dt_fixed,
                         filename='line_output_u,v_p.csv')
                         
                         
integral_output=LineIntegralGauges(gauges=((('vof'), (((2.50, 1.50), (2.50, 0.)),
                                               ((5.50, 1.50), (5.50, slope*5.50)),
                                               ((6.50, 1.50), (6.50, y_slop_end)))),
                                           ),
                         activeTime=(0, T),
                         sampleRate=1/dt_fixed,
                         filename='integral_output_vof.csv')



##########################################
# Numerical Options and other parameters #
##########################################

nd = 2

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,-9.8]

# ------------------------------------------------ #


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

he = tank_dim[0]/500
domain.writePoly("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

quad_order = 3

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = False
openSides = False
openEnd = False
smoothBottom = False
smoothObstacle = False
movingDomain=False
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=False

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 100
dt_fixed = period/21.0
dt_init = min(0.1*dt_fixed,0.001)
T = 50*period
nDTout= int(round(T/dt_fixed))
runCFL = 0.9

#----------------------------------------------------
water_depth  = waterLevel
depth = waterLevel

#  Discretization -- input options
genMesh=True
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.25
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.25
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.25
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-10,0.001*he**2)
vof_nl_atol_res = max(1.0e-10,0.001*he**2)
ls_nl_atol_res = max(1.0e-10,0.001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4

# Initial condition
waterLine_x = 2*tank_dim[0]
waterLine_z = waterLevel

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[1]-waterLine_z 
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)


#solution variables

def wavePhi(x,t):
    return x[1] - aa(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

