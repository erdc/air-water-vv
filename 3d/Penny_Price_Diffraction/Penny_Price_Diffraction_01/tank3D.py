from __future__ import print_function
from __future__ import division
from builtins import map
from builtins import zip
from builtins import range
from past.utils import old_div
from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges

#wave generator
windVelocity = (0.0,0.0,0.0)
inflowHeightMean = 1.0
inflowVelocityMean = (0.0,0.0,0.0)
period = 1.94
omega = 2.0*math.pi/period
#waveheight = 0.10
#amplitude = waveheight/ 2.0
wavelength =  5.0
k = 2.0*math.pi/wavelength
 

#  Discretization -- input options  

genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=False
timeDiscretization='be'#'vbdf'#'be','flcbdf'
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
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()    
    
if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES) 
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
    sys.exit()
    
#  Discretization   
nd = 3
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
L = (float(6.0*wavelength), 30.0, 1.50)
x1=2.0*wavelength
x2=x1+0.01
y1=2.0*wavelength

he = old_div(wavelength,20)

GenerationZoneLength = wavelength*1.0
AbsorptionZoneLength= wavelength*2.0
spongeLayer = True #False  
levee=spongeLayer
slopingSpongeLayer=spongeLayer
xSponge = GenerationZoneLength
xRelaxCenter = old_div(xSponge,2.0)
epsFact_solid = old_div(xSponge,2.0)
#zone 2
xSponge_2 = L[0]-AbsorptionZoneLength
ySponge_2= L[1]-AbsorptionZoneLength
xRelaxCenter_2 = 0.5*(xSponge_2+L[0])
yRelaxCenter_2 = 0.5*(ySponge_2+L[1])
epsFact_solid_2 = old_div(AbsorptionZoneLength,2.0)

nLevels = 1
weak_bc_penalty_constant = 100.0
quasi2D=False
if quasi2D:#make tank one element wide
    L = (L[0],he,L[2])

#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

structured=False  

gauge_dx=0.5
PGL=[]
LGL=[]
for i in range(0,int(old_div(L[0],gauge_dx)+1)): #+1 only if gauge_dx is an exact 
  PGL.append([gauge_dx*i,old_div(L[1],2.0),0.5])
  LGL.append([(gauge_dx*i,old_div(L[1],2.0),0),(gauge_dx*i,old_div(L[1],2.0),L[2])])
 

gaugeLocations=tuple(map(tuple,PGL)) 
columnLines=tuple(map(tuple,LGL)) 


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
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back','obst']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
        domain = Domain.RectangularDomain(L)
    elif spongeLayer:
        vertices=[[0.0,0.0,0.0],#0
                  [xSponge,0.0,0.0],#1  
                  [xSponge_2,0.0,0.0],#2 
                  [L[0],0.0,0.0],#3
                  [L[0],y1,0.0], #4
                  [L[0],ySponge_2,0.0],#5
                  [L[0],L[1],0.0],#6
                  [xSponge_2,L[1],0.0],#7
                  [x2,L[1],0.0],#8
                  [x2,ySponge_2,0.0],#9
                  [x2,y1,0.0],#10
                  [x1,y1,0.0],#11
                  [x1,L[1],0.0],#12
                  [xSponge,L[1],0.0],#13
                  [0.0,L[1],0.0], #14
                  [0.0,y1,0.0], #15
                  [xSponge,y1,0.0], #16
                  [xSponge_2,y1,0.0], #17
                  [xSponge_2,ySponge_2,0.0]] #18
             
        vertexFlags=[]
        for i in range(0,int(len(vertices))):
            vertexFlags.append(boundaryTags['bottom'])
        
        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1],L[2]])
            vertexFlags.append(boundaryTags['top'])
       
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
                  [11,12],
                  [12,13],
                  [13,14],
                  [14,15],
                  [15,0],
                  [1,16],
                  [16,13],
                  [15,16],
                  [16,11],
                  [2,17],
                  [17,18],
                  [18,7],                   
                  [10,17],
                  [17,4],
                  [9,18],
                  [18,5]]
                 
        segmentFlags=[boundaryTags['front'],
                     boundaryTags['front'],
                     boundaryTags['front'],                 
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['back'],
                     boundaryTags['back'],
                     boundaryTags['obst'],
                     boundaryTags['obst'],
                     boundaryTags['obst'],
                     boundaryTags['obst'],
                     boundaryTags['back'],
                     boundaryTags['back'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty']]
        

        facets=[]
        facetFlags=[]
        print(int(len(vertices)))
        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+int(old_div(len(vertices),2)),s[0]+int(old_div(len(vertices),2))]])
            facetFlags.append(sF)

        bf=[[0,1,16,15],[1,2,17,10,11,16],[2,3,4,17],[15,16,13,14],[16,11,12,13],[10,17,18,9],[9,18,7,8],[17,4,5,18],[18,5,6,7]]
        tf=[]
 
        for i in range(0,int(len(bf))):
         facets.append([bf[i]])
         tf=[ss + int(old_div(len(vertices),2)) for ss in bf[i]]
         facets.append([tf])

        for i in range(0,int(len(bf))):
         facetFlags.append(boundaryTags['bottom'])
         facetFlags.append(boundaryTags['top'])

        print(facets)
        print(facetFlags)

        regions=[[xRelaxCenter, 0.01,0.0],#1
                 [xRelaxCenter_2, 0.01, 0.0],#2
                 [xRelaxCenter, 0.99*L[1],0.0],#3
                 [xRelaxCenter_2, y1+0.1,0.0],#4
                 [xRelaxCenter_2, 0.99*L[1],0.0],#5
                 [x2+0.01,yRelaxCenter_2,0.0],#6
                 [x1,0.01, 0.0]]#7

        regionFlags=[1,2,3,4,5,6,7]

        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags)

        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % (old_div((he**3),6.0),)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

        porosityTypes      = numpy.array([1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0])

        dragAlphaTypes = numpy.array([0.0, #1
                                      old_div(0.5,1.004e-6),#1
                                      old_div(0.5,1.004e-6), #2
                                      0.0, #2
                                      0.0, #3
                                      old_div(0.5,1.004e-6),#3
                                      old_div(0.5,1.004e-6),#4
                                      0.0,#4
                                      old_div(0.5,1.004e-6),#5
                                      0.0,#5
                                      old_div(0.5,1.004e-6),#6
                                      0.0])#6

        dragBetaTypes = numpy.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        
        epsFact_solidTypes = np.array([0.0,epsFact_solid,epsFact_solid_2,0.0,0.0,epsFact_solid,epsFact_solid_2,0.0,epsFact_solid_2,0.0,epsFact_solid_2,0.0])

    else:             
        vertices=[[0.0,0.0,0.0],#0
                  [L[0],0.0,0.0],#1
                  [L[0],L[1],0.0],#2       
                  [0.0,L[1],0.0]]#3
        
               
        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom']]


        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1],L[2]])
            vertexFlags.append(boundaryTags['top'])

        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,0]]
                 
        segmentFlags=[boundaryTags['front'],                   
                     boundaryTags['right'],
                     boundaryTags['back'],
                     boundaryTags['left']]

        facets=[]
        facetFlags=[]

        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+4,s[0]+4]])
            facetFlags.append(sF)

        bf=[[0,1,2,3]]
        tf=[]
        for i in range(0,1):
         facets.append([bf[i]])
         tf=[ss + 4 for ss in bf[i]]
         facets.append([tf])

        for i in range(0,1):
         facetFlags.append(boundaryTags['bottom'])
         facetFlags.append(boundaryTags['top'])

        for s,sF in zip(segments,segmentFlags):
            segments.append([s[1]+4,s[0]+4])
            segmentFlags.append(sF)
        

        regions=[[0.5*L[0],0.5*L[1], 0.0]]
        regionFlags=[1]

        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % (old_div((he**3),6.0),)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

# Time stepping
T=20.0*period
dt_fixed = T
dt_init = min(0.1*dt_fixed,0.1*he)
runCFL=0.90
nDTout = int(round(old_div(T,dt_fixed)))

# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.25
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
    epsFact_consrv_diffusion = 0.1
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
    epsFact_consrv_diffusion = 1.0
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
g = [0.0,0.0,-9.8]

# Initial condition
waterLine_x =  2*L[0]
waterLine_z =  inflowHeightMean
waterLine_y =  2*L[1]

def signedDistance(x):
    phi_z = x[2]-waterLine_z 
    return phi_z

def theta(x,t):
    return k*x[0] - omega*t + old_div(math.pi,2.0)

def z(x):
    return x[2] - inflowHeightMean

def ramp(t):
  t0=10 #ramptime
  if t<t0:
    return 1
  else:
    return 1 

h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0
sigma = omega - k*inflowVelocityMean[0]

Y =  [ 0.06239652,     #Surface elevation Fourier coefficients for non-dimensionalised solution
          0.00364527,     
          0.00024734,    
          0.00001974,  
          0.00000174,
          0.00000016,    
          0.00000002] 

def waveHeight(x,t):
   waterDepth = inflowHeightMean 
   for i in range(0,int(len(Y))):  
       waterDepth += Y[i]*cos((i+1)*theta(x,t))/k
   return waterDepth*ramp(t)
 
B = [0.06758337,     
     0.00125642,    
    -0.00000264,    
    -0.00000062]    

def waveVelocity_u(x,t):
   wu=0
   for i in range(0,int(len(B))): 
     wu += sqrt(old_div(abs(g[2]),k))*(i+1)*B[i]*cosh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*cos((i+1)*theta(x,t))
    
   return wu*ramp(t)

def waveVelocity_v(x,t):
   wv=0
   for i in range(0,int(len(B))): 
     wv += sqrt(old_div(abs(g[2]),k))*(i+1)*B[i]*sinh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*sin((i+1)*theta(x,t)) 

   return wv*ramp(t)
  
#solution variables

def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

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

def twpflowVelocity_w(x,t):
    return 0.0

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - inflowHeightMean)

def outflowPressure(x,t):
  if x[2]>inflowHeightMean:
    return (L[2]-x[2])*rho_1*abs(g[2])
  else:
    return (L[2]-inflowHeightMean)*rho_1*abs(g[2])+(inflowHeightMean-x[2])*rho_0*abs(g[2])

def waterVelocity(x,t):
   if x[2]>inflowHeightMean:
     return 0.0
   else: 
     ic=inflowVelocityMean[0]
     return ic

def zeroVel(x,t):
    return 0.0

def phi_solid(x,t):
  if x[0]<=x2:
    return xRelaxCenter-x[0] #region 1 and 3
  elif x[1]<=ySponge_2:   
    return xRelaxCenter_2-x[0] #region 2 and 4
  elif x[0]<=xSponge_2:
    return yRelaxCenter_2-x[1] #region 6
  else:
    return epsFact_solid_2 -sqrt((x[0]-xSponge_2)**2 + (x[1]-ySponge_2)**2) #region 5


from collections import  namedtuple

RelaxationZone = namedtuple("RelaxationZone","phi_solid sign u v w")

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
                if mType in self.zones:
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN,k]
                        m.coefficients.q_phi_solid[eN,k] = self.zones[mType].sign*self.zones[mType].phi_solid(x,t)
                        m.coefficients.q_velocity_solid[eN,k,0] = self.zones[mType].u(x,t)
                        m.coefficients.q_velocity_solid[eN,k,1] = self.zones[mType].v(x,t)
                        m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
        m.q['phi_solid'] = m.coefficients.q_phi_solid
        m.q['velocity_solid'] = m.coefficients.q_velocity_solid

rzWaveGenerator = RelaxationZoneWaveGenerator(zones={1:RelaxationZone(phi_solid,
                                                                      1.0,
                                                                      twpflowVelocity_u,
                                                                      twpflowVelocity_v,
                                                                      twpflowVelocity_w),
                                                    2:RelaxationZone(phi_solid,
                                                                     -1.0,
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel),
                                                    3:RelaxationZone(phi_solid,
                                                                      1.0,
                                                                      twpflowVelocity_u,
                                                                      twpflowVelocity_v,
                                                                      twpflowVelocity_w),
                                                    4:RelaxationZone(phi_solid,
                                                                     -1.0,
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel),
                                                    5:RelaxationZone(phi_solid,
                                                                     -1.0,
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel),
                                                    6:RelaxationZone(phi_solid,
                                                                     -1.0,
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel)})

