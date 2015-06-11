from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
import WaveTools


#wave generator
windVelocity = (0.0,0.0)
inflowHeightMean = 0.86
inflowVelocityMean = (0.0,0.0)
period = 1.43
omega = 2.0*math.pi/period
#waveheight = 
#amplitude = waveheight/ 2.0
wavelength =  3.037
k = 2.0*math.pi/wavelength

#  Discretization -- input options  
genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=False
timeDiscretization='be'#'be','vbdf','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
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
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
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
    
# Domain and mesh

L = (54.25,1.26)
he = float(wavelength)/65
x1=9.22 
x2=10.42 
x3=15.01 
x4=27.04 
x5=31.04 
x6=37.07 
x7=44.61 
x8=45.81
y1=0.06
y2=0.66

GenerationZoneLength = 4.0
AbsorptionZoneLength= 6.0
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0
#zone 2
xSponge_2 = L[0]-AbsorptionZoneLength
xRelaxCenter_2 = 0.5*(xSponge_2+L[0])
epsFact_solid_2 = AbsorptionZoneLength/2.0


weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False


gaugeLocations=((0.0,0.36,0.0),
                (xRelaxCenter,0.36,0.0),
                (xSponge,0.36,0.0),
                (24.04,0.5295,0.0),
                (30.04,0.66,0.0),
                (34.04,0.36,0.0),
                (xSponge_2,0.36,0.0),
                (xRelaxCenter_2,0.36,0.0),
                (L[0],0.36,0.0)) 


columnLines=( ((0.0,0.0,0.0),(0.0,L[1],0.0)),
              ((xRelaxCenter,0.0,0.0),(xRelaxCenter,L[1],0.0)),
              ((xSponge,0.0,0.0),(xSponge,L[1],0.0)),
              ((24.04,0.5295,0.0),(24.04,L[1],0.0)),
              ((30.04,0.66,0.0),(30.04,L[1],0.0)),
              ((34.04,0.362,0.0),(34.04,L[1],0.0)),
              ((xSponge_2,0.0,0.0),(xSponge_2,L[1],0.0)),
              ((xRelaxCenter_2,0.0,0.0),(xRelaxCenter_2,L[1],0.0)),
              ((L[0],0.0,0.0),(L[0],L[1],0.0)))


pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
                                (('p',),    gaugeLocations)),
                  activeTime = (0, 1000.0),
                  sampleRate = 0,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt')

fields = ('vof',)

columnGauge = LineIntegralGauges(gauges=((fields, columnLines),),
                                 fileName='column_gauge.csv')

#lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_y=0':((0.0,0.0,0.0),(L[0],0.0,0.0))},linePoints=24)

#lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)


if useHex:   
    nnx=ceil(L[0]/he)+1
    nny=ceil(L[1]/he)+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=ceil(L[0]/he)+1
        nny=ceil(L[1]/he)+1
    elif spongeLayer:
        vertices=[[0.0,0.0],#0
                  [x1,0.0],#1
                  [x2,y1],#2
                  [x3,y1],#3
                  [x4,y2],#4
                  [x5,y2],#5
                  [x6,y1],#6
                  [x7,y1],#7
                  [x8,0.0],#8
                  [L[0],0.0],#9
                  [L[0],L[1]],#10
                  [0.0,L[1]],#11
                  [xSponge,0.0],#12
                  [xSponge,L[1]],#13
                  [xSponge_2,0.0],#14
                  [xSponge_2,L[1]]]#15

        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],                    
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],  
                     boundaryTags['top'],  
                     boundaryTags['bottom'],  
                     boundaryTags['top'],
                     boundaryTags['bottom'],  
                     boundaryTags['top']]

        segments=[[0,12],
                  [12,1],
                  [1,2],
                  [2,3],
                  [3,4],
                  [4,5],
                  [5,6],  
                  [6,7],
                  [7,8],
                  [8,14],
                  [14,9],
                  [9,10],
                  [10,15],
                  [15,13],
                  [13,11],
                  [11,0],
                  [12,13],
                  [14,15]]

        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],                    
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['left'],
                      0,
                      0]
    
        regions=[[xRelaxCenter, 0.5*L[1]],
                 [xRelaxCenter_2, 0.5*L[1]],
                 [0.5*L[0],0.5*L[1]]]
        regionFlags=[1,2,3]
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)

        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
        porosityTypes      = numpy.array([1.0,
                                          1.0,
                                          1.0,
                                          1.0])
        dragAlphaTypes = numpy.array([0.0,
                                      0.5/1.004e-6,
                                      0.5/1.004e-6,
                                      0.0])
        dragBetaTypes = numpy.array([0.0,0.0,0.0,0.0])
        
        epsFact_solidTypes = np.array([0.0,epsFact_solid,epsFact_solid_2,0.0])

    else:
        vertices=[[0.0,0.0],#0
                  [L[0],0.0],#1
                  [L[0],L[1]],#2
                  [0.0,L[1]]]#3

        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],  
                     boundaryTags['top']]
        segments=[[0,1],
                  [1,2],
                  [2,3],  
                  [3,0]
                  ]
        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['left']]

        regions=[ [ 0.1*L[0] , 0.1*L[1] ],
                  [0.95*L[0] , 0.95*L[1] ] ]
        regionFlags=[1,2]
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)

        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
# Time stepping
T=60*period
dt_fixed = T
dt_init = min(0.1*dt_fixed,0.1)
runCFL=0.9
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.0
if useMetrics:
    ns_shockCapturingFactor  = 0.25
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.35
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.35
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.75
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 1.5
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
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

ns_nl_atol_res = max(1.0e-10,0.00001*he**2)
vof_nl_atol_res = max(1.0e-10,0.00001*he**2)
ls_nl_atol_res = max(1.0e-10,0.0001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.0001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
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

# Initial condition
waterLine_x = 2*L[0]
waterLine_z = inflowHeightMean


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
    return k*x[0] - omega*t + pi/2.0

def z(x):
    return x[1] - inflowHeightMean

sigma = omega - k*inflowVelocityMean[0]
h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0

Y =  [0.08448147, #Surface elevation Fourier coefficients for non-dimensionalised solution
         0.00451131,
         0.00032646,
         0.00002816,
         0.00000267,
         0.00000027,
         0.00000003,    
         0.00000001] 

def waveHeight(x,t):
  waterDepth = inflowHeightMean 
  for i in range(0,int(len(Y))):  
       waterDepth += Y[i]*cos((i+1)*theta(x,t))/k
  return waterDepth
 
B = [0.08677594,
     0.00070042,
    -0.00001289,
     0.00000006,     
     0.00000001]  
 
def waveVelocity_u(x,t):
   wu=0
   for i in range(0,int(len(B))): 
     wu += sqrt(abs(g[1])/k)*(i+1)*B[i]*cosh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*cos((i+1)*theta(x,t))
    
   return wu

def waveVelocity_v(x,t):
   wv=0
   for i in range(0,int(len(B))): 
     wv += sqrt(abs(g[1])/k)*(i+1)*B[i]*sinh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*sin((i+1)*theta(x,t)) 

   return wv

 
#solution variables

def wavePhi(x,t):
    return x[1] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = waveVelocity_v(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

outflowHeight=inflowHeightMean

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)

def outflowPhi(x,t):
    return x[1] - outflowHeight

def outflowPressure(x,t):
  if x[1]>inflowHeightMean:
    return (L[1]-x[1])*rho_1*abs(g[1])
  else:
    return (L[1]-inflowHeightMean)*rho_1*abs(g[1])+(inflowHeightMean-x[1])*rho_0*abs(g[1])


    #p_L = L[1]*rho_1*g[1]
    #phi_L = L[1] - outflowHeight
    #phi = x[1] - outflowHeight
    #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
    #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def twpflowVelocity_w(x,t):
    return 0.0

def zeroVel(x,t):
    return 0.0

from collections import  namedtuple

RelaxationZone = namedtuple("RelaxationZone","center_x sign u v w")

class RelaxationZoneWaveGenerator(AV_base):
    """ Prescribe a velocity penalty scaling in a material zone via a Darcy-Forchheimer penalty
    
    :param zones: A dictionary mapping integer material types to Zones, where a Zone is a named tuple
    specifying the x coordinate of the zone center and the velocity components
    """
    def __init__(self,zones):
        assert isinstance(zones,dict)
        self.zones = zones
    def calculate(self):
        for l,m in enumerate(self.model.levelModelList):
            for eN in range(m.coefficients.q_phi.shape[0]):
                mType = m.mesh.elementMaterialTypes[eN]
                if self.zones.has_key(mType):
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN,k]
                        m.coefficients.q_phi_solid[eN,k] = self.zones[mType].sign*(self.zones[mType].center_x - x[0])
                        m.coefficients.q_velocity_solid[eN,k,0] = self.zones[mType].u(x,t)
                        m.coefficients.q_velocity_solid[eN,k,1] = self.zones[mType].v(x,t)
                        #m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
        m.q['phi_solid'] = m.coefficients.q_phi_solid
        m.q['velocity_solid'] = m.coefficients.q_velocity_solid

rzWaveGenerator = RelaxationZoneWaveGenerator(zones={1:RelaxationZone(xRelaxCenter,
                                                                      -1.0,
                                                                      twpflowVelocity_u,
                                                                      twpflowVelocity_v,
                                                                      twpflowVelocity_w),
                                                    2:RelaxationZone(xRelaxCenter_2,
                                                                     1.0,
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel)})


 

