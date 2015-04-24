from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
#import WaveTools

airVent=True
#wave generator
windVelocity = (0.0,0.0)
inflowHeightMean = 1.5
outflowHeightMean = 0.50
inflow_velocity = 1.345
outflow_velocity=4.035


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
L = (13.5 , 2.1)

#Obstacle or weir dimensions/position
obst_portions = (3.5,0.5) #(width,height)

obst_x_start = 6.0 # start x coordinate of the obstacle; caution to be in the domain's range
cot_upstream_slope=2.0
obst_high_point_x= obst_x_start + obst_portions[1]*cot_upstream_slope #x coordinate of the highest point of the weir
obst_x_end = obst_x_start + obst_portions[0] # end x coordinate of the obstacle; caution to be in the domain's range
obst = (obst_x_start,obst_portions[1],obst_x_end,obst_high_point_x) #coordinates of the obstacle to be used to define the boundary


#Background refinement
Refinement = 68 #0.05
he = L[0]/float(4*Refinement-1)

# Refinement parameters
#x_refine = (1.0 , 1.5 , 3.0) #end of zone 1, end of zone 2, end of zone 3 (zone 4 is up to the right wall)
#refinementLevel = (4 , 2) #refinemnt level for zone 1 and zone 2 and 4 respectively zone 3 has the basic refinement level

#Left boundary imposed velocity 
GenerationZoneLength = 2.0
AbsorptionZoneLength= 2.0
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0
#zone 2
xSponge_2 = L[0]-AbsorptionZoneLength
xRelaxCenter_2 = 0.5*(xSponge_2+L[0])
epsFact_solid_2 = AbsorptionZoneLength/2

weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False

gauge_dx=0.25
PGL=[]
LGL=[]
for i in range(0,int(L[0]/gauge_dx+1)): #+1 only if gauge_dx is an exact 
  PGL.append([gauge_dx*i,0.5,0])
  if gauge_dx*i<obst_x_start or gauge_dx*i>obst_x_end:
      a = 0.0
  elif gauge_dx*i>=obst_x_start and gauge_dx*i<obst_high_point_x:
      a= (obst[1]/(obst_high_point_x-obst_x_start)*(gauge_dx*i-obst_x_start))
  else: 
      a=obst[1]+(obst[1]/(obst_x_end-obst_high_point_x)*(gauge_dx*i-obst_high_point_x))
 # print a , gauge_dx*i
  
  LGL.append([(gauge_dx*i,a,0),(gauge_dx*i,L[1],0)])
 

gaugeLocations=tuple(map(tuple,PGL)) 
columnLines=tuple(map(tuple,LGL)) 


pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
                                (('p',),    gaugeLocations)),
                  activeTime = (0, 1000.0),
                  sampleRate = 0,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt')

#print gaugeLocations
#print columnLines

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
    boundaries=['left','right','bottom','top','front','back','airvent']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=ceil(L[0]/he)+1
        nny=ceil(L[1]/he)+1
    elif spongeLayer:
        vertices=[[0.0,0.0],#0
                  [obst[0],0.0], #1
                  [obst[3],obst[1]], #2   
                  [obst[2],0.0],#3 
                  [L[0],0.0],#4
                  [L[0],L[1]],#5
                  [0.0,L[1]], #6                 
                  [xSponge,0.0],#7
                  [xSponge,L[1]],#8
                  [xSponge_2,0.0],#9
                  [xSponge_2,L[1]]]#10

        vertexFlags=[boundaryTags['bottom'],
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

        segments=[[0,7],
                  [7,1],
                  [1,2],                
                  [2,3],
                  [3,9],
                  [9,4],                  
                  [4,5],
                  [5,10],
                  [10,8],            
                  [8,6],                  
                  [6,0],
                  [7,8],
                  [9,10]]

        segmentFlags=[boundaryTags['bottom'],
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
T=30.0
dt_fixed = 0.25
dt_init = min(0.1*dt_fixed,0.1)
runCFL=0.75
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = True
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.75
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.75
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.75
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.75
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
waterLine_x = obst_x_start
waterLine_x2 = obst_x_end 
waterLine_z = inflowHeightMean
waterLine_z2 = outflowHeightMean

def signedDistance(x):
   if x[0]<waterLine_x:
      phi_z = x[1]-waterLine_z
   elif x[0]<waterLine_x2:
      phi_z = x[1]-waterLine_z+((waterLine_z-waterLine_z2)/(waterLine_x2-waterLine_x)*(x[0]-waterLine_x))
   else: 
      phi_z=x[1]-waterLine_z2
 
   return phi_z
 
#solution variables

def wavePhi(x,t):
    return x[1] - inflowHeightMean

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = inflow_velocity
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = 0.0
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowVelocity_w(x,t):
    return 0.0

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowPressure(x,t):
 # if x[1]>inflowHeightMean:
    return (L[1]-x[1])*rho_1*abs(g[1])
 # else:
 #   return (L[1]-inflowHeightMean)*rho_1*abs(g[1])+(inflowHeightMean-x[1])*rho_0*abs(g[1])


    #p_L = L[1]*rho_1*g[1]
    #phi_L = L[1] - outflowHeight
    #phi = x[1] - outflowHeight
    #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
    #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
def outflowPhi(x,t):
   return x[1] - outflowHeightMean

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,outflowPhi(x,t))

def outflowVel(x,t):
    waterspeed = outflow_velocity
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,outflowPhi(x,t)-epsFact_consrv_heaviside*he)
    u = (1.0-H)*waterspeed
    return u
 
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
                                                                      -1.0, #change after compliling to 1.0
                                                                      twpflowVelocity_u,
                                                                      twpflowVelocity_v,
                                                                      twpflowVelocity_w),
                                                    2:RelaxationZone(xRelaxCenter_2,
                                                                     1.0, #change after compliling to -1.0
                                                                     outflowVel,
                                                                     zeroVel,
                                                                     zeroVel)})



