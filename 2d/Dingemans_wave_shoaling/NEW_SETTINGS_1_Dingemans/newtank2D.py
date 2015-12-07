from proteus import Domain
from proteus import SpatialTools as st
import ode
from math import *
import numpy as np


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank_dim = [58.0, 1.26]
y1=0.06
y2=0.66
L_leftSpo = 1.5
L_rightSpo = 1.5
b_or = {'bottom': [0., -1.],
        'right': [1., 0.],
        'top': [0., 1.],
        'left': [-1., 0.],
        'sponge': None,
        }
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,
                'sponge': 5,
               }
vertices=[[0.0,0.0],#0
          [9.22,0.0],#1
          [10.42,y1],#2
          [15.01,y1],#3
          [27.04,y2],#4
          [31.04,y2],#5
          [37.07,y1],#6
          [44.61,y1],#7
          [45.81,0.0],#8
          [tank_dim[0],0.0],#9
          [tank_dim[0],tank_dim[1]],#10
          [0.0,tank_dim[1]],#11
          [L_leftSpo,0.0],#12
          [L_leftSpo,y2],#13
          [tank_dim[0]-L_rightSpo,0.0],#14
          [tank_dim[0]-L_rightSpo,y2],#15
         ]

vertexFlags=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 1, 3, 1, 3])

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,8],
          [8,9],
          [9,10],
          [10,11],
          [11,0],
          [12,13],
          [14,15],
         ]

segmentFlags=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 5])
regions=[ [ 0.1*tank_dim[0] , 0.1*tank_dim[1] ],
          [ 0.5*tank_dim[0] , y2+0.1 ],
          [0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]
regionFlags=np.array([1,2,3])

left = right = False
if L_leftSpo is not None: left = True
if L_rightSpo is not None: right = True


tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=b_or)

abszone=tank.setAbsorptionZones(indice=[0, 2], epsFact_solid=L_leftSpo)



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
g = np.array([0.0,-9.8])



# ----- WAVE CONDITIONS ----- #
period=1.43
waterLevel = 0.86
Ycoeff=[0.08448147, 0.00451131, 0.00032646, 0.00002816,
         0.00000267, 0.00000027, 0.00000003, 0.00000001]
Bcoeff=[0.08677594, 0.00070042, -0.00001289,
            0.00000006, 0.00000001]
waveDir=np.array([1, 0])
mwl=0.
waveHeight=np.array(Ycoeff)
omega = float(2.0*pi/period)
wavelength = 8.12
inflowHeightMean=waterLevel
inflowVelocityMean = (0.0,0.0)
windVelocity = (0.0,0.0)
k = float(2.0*pi/wavelength)



#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 100
dt_fixed = period/21.0
dt_init = min(0.1*dt_fixed,0.001)
T = 50*period
nDTout= int(round(T/dt_fixed))
runCFL = 0.9



# ----- WAVE input ----- #

from proteus import WaveTools as wt

waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=wavelength,
                                  waveType="Linear",
                                  Ycoeff = Ycoeff,
                                  Bcoeff = Bcoeff,
                                  )




# ----- Output Gauges ----- #

from proteus import Gauges as ga

line_output=ga.LineGauges(gauges=((('u', 'v'), (((0., 1.26), (0., 0.)),
                                                ((20.04, 1.26), (20.04, 0.66)),
                                                ((30.04, 1.26), (30.04, 0.66)),
                                                ((36.04, 1.26), (36.04, 0.66))),
                                   ),
                                  (('p'), (((0., 1.26), (0., 0.)),
                                           ((20.04, 1.26), (20.04, 0.66)),
                                           ((30.04, 1.26), (30.04, 0.66)),
                                           ((36.04, 1.26), (36.04, 0.66))),
                                   ),
                                  ),
                          activeTime=(0, 71.5),
                          sampleRate=1/dt_fixed,
                          fileName='line_gauges.csv')
                         
                         
integral_output=ga.LineIntegralGauges(gauges=((('vof'), (((20.04, 1.26), (20.04, 0.66)),
                                                         ((30.04, 1.26), (30.04, 0.66)),
                                                         ((36.04, 1.26), (36.04, 0.66)))),
                                             ),
                                      activeTime=(0, 71.5),
                                      sampleRate=1/dt_fixed,
                                      fileName='line_integral_gauges.csv')




# ----- Mesh ----- #

from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

he = tank_dim[0]/2900
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
openEnd = True
smoothBottom = False
smoothObstacle = False
movingDomain=False
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True


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
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
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
waterLine_z = inflowHeightMean



def waveHeight(x,t):
   Y =  [0.08448147,    #Surface elevation Fourier coefficients for non-dimensionalised solution
         0.00451131,
         0.00032646,
         0.00002816,
         0.00000267,
         0.00000027,
         0.00000003,    
         0.00000001] 
   waterDepth = inflowHeightMean 
   for i in range(0,int(len(Y))):  
       waterDepth += Y[i]*cos((i+1)*theta(x,t))/k
   return waterDepth


def wavePhi(x,t):
    return x[1] - waveHeight(x,t)


# ----- BOUNDARY CONDITIONS ----- #

tank.BC.top.setOpenAir()
tank.BC.left.MonochromaticWaveBoundary(waveinput)
tank.BC.bottom.setNoSlip()
tank.BC.right.reset()
