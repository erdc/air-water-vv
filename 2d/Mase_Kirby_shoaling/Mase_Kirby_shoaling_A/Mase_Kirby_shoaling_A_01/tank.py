from math import *
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
from proteus.WaveTools import RandomWaves
from proteus.mprans import SpatialTools as st
from proteus.mprans import BoundaryConditions as BC
from proteus import MeshTools, AuxiliaryVariables
import proteus.MeshTools
import ode
import numpy as np


# Wave generator
windVelocity = [0.,0.,0.]
inflowHeightMean = 0.47
period = 1.0
waveheight = 0.062
amplitude = waveheight/2.0
wavelength = 1.56
waveDir = np.array([1,0,0])
waterLevel = inflowHeightMean
g = np.array([0,-9.81,0])
Tp = period
Hs = waveheight
depth = waterLevel
mwl = inflowHeightMean
N = 101
bandFactor = 2.0
spectName = "PM_mod"
spectral_params = None
phi = None
waves = RandomWaves(Tp,
                    Hs,
                    mwl,
                    depth,
                    waveDir,
                    g,
                    N,
                    bandFactor,
                    spectName,
                    spectral_params,
                    phi
                    )
             

# Discretization -- input options  
genMesh = True
movingDomain = False
applyRedistancing = True
checkMass = False
freezeLevelSet = False
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'
spaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
     

       
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
    


# Domain
domain = Domain.PlanarStraightLineGraphDomain()



# Shape
L = [24.0, 1.0]
he = wavelength/100 #try this first
domain.MeshOptions.elementSize(he)
GenerationZoneLength = 0.50*wavelength
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0
cot_slope = 20.0
x1 = 10.0
y1 = 14.0*float(1.0/cot_slope)


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


vertices = [[0.0,0.0],       #0
            [x1,0.0],        #1
            [L[0],y1],       #2
            [L[0],L[1]],     #3
            [0.0,L[1]],      #4
            [xSponge,0.0],   #5
            [xSponge,L[1]],  #6
            ]


vertexFlags = np.array([1, 1, 1, 3, 3, 1, 3])


segments = [[0,5],
            [5,1],
            [1,2],
            [2,3],
            [3,6],
            [6,4],
            [4,0],
            [5,6],
            ]


segmentFlags = np.array([1, 1, 1, 2, 3, 3, 4, 5])


regions = [ [ 0.1*L[0] , 0.1*L[1] ],
            [ 0.95*L[0] , 0.95*L[1] ] ]


regionFlags = np.array([1,2])


tank = st.CustomShape(domain,
                      vertices=vertices,
                      vertexFlags=vertexFlags,
                      segments=segments,
                      segmentFlags=segmentFlags,
                      regions=regions,
                      regionFlags=regionFlags,
                      boundaryTags=boundaryTags,
                      boundaryOrientations=b_or)


tank.setGenerationZones(flags=1,
                        epsFact_solid=xSponge/2.0,
                        center=[xSponge, 0],
                        orientation=[1, 0, 0],
                        waves=waves,
                        windSpeed=windVelocity,
                        dragAlphaTypes=0.5/1.005e-6,
                        dragBetaTypes=0.,
                        porosityTypes=1.)


# Physical parameters
rho_0 = 998.2
nu_0 = 1.004e-6

rho_1 = 1.205
nu_1 = 1.500e-5 

sigma_01 = 0.0


weak_bc_penalty_constant = 10.0/nu_0
nLevels = 1
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

st.assembleDomain(domain)
domain.writePoly("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes = False
quad_order = 3


# Time stepping
T = 20*period
dt_fixed = 0.05 #T  
dt_init = 0.001
runCFL = 0.9
nDTout = int(round(T/dt_fixed))


# GAUGES  ------------------------------------------------------
gauge_dx=2 #0.25
PGL=[]
LGL=[]
for i in range(0, int(L[0]/gauge_dx)): #+1 only if gauge_dx is an exact 
  if gauge_dx*i<x1:
      a = 0.0
  elif gauge_dx*i>=x1: 
      a = (1/cot_slope)*(gauge_dx*i-x1) 
 
  LGL.append([[gauge_dx*i,a,0],[gauge_dx*i,L[1],0]])
  PGL.append([gauge_dx*i,a,0])                          


gaugeLocations=tuple(map(tuple,PGL))

columnLines=tuple(map(tuple,LGL))


pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
                                (('p',), gaugeLocations)),
                  activeTime = (0, 20.0),
                  sampleRate = 0.,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt')

fields = (('vof',))

columnGauge = LineIntegralGauges(gauges=((fields, columnLines),),
                                 activeTime = (0, 20.0),
                                 sampleRate = 0.,
                                 fileName='column_gauge.csv')


domain.auxiliaryVariables += [pointGauges] 


# Numerical parameters
ns_forceStrongDirichlet = False    #True
backgroundDiffusionFactor = 0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor = 0.5
    rd_lag_shockCapturing = False
    epsFact_density = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0


ns_nl_atol_res = max(1.0e-10,0.001*he**2)
vof_nl_atol_res = max(1.0e-10,0.001*he**2)
ls_nl_atol_res = max(1.0e-10,0.001*he**2)
rd_nl_atol_res = max(1.0e-10,0.01*he)
mcorr_nl_atol_res = max(1.0e-10,0.0001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)


# Turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4
   
 
# Initial condition
waterLine_x = 2*L[0]
waterLine_z = inflowHeightMean


# Boundary conditions
tank.BC.top.setOpenAir()
tank.BC.bottom.setFreeSlip()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waves, vert_axis=1, windSpeed=windVelocity, air=1., water=0.)
tank.BC.right.setFreeSlip()
tank.BC.sponge.setNonMaterial()


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


def theta(x,t):
    return k*x[0] - omega*t


def z(x):
    return x[1] - inflowHeightMean


def wavePhi(x,t):
    return x[1] - inflowHeightMean - waves.eta(x,t)


