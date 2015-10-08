from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
from proteus import WaveTools as WT

from proteus import Context
opts=Context.Options([
    ("wave_type", 'linear', "type of waves generated: 'linear', 'Nonlinear', 'single-peaked', 'double-peaked'"),
    ("depth", 0.457, "water depth at leading edge of vegetation (not at wave generator)[m]"),
    ("wave_height", 0.192, "wave height at leading edget of vegetation [m]"),
    ("peak_period", 2.0, "Peak period [s]"),
    ("peak_period2", 1.5, "Second peak period (only used in double-peaked case)[s]"),
    ("peak_wavelength",3.91,"Peak wavelength in [m]"),
    ("parallel", False, "Run in parallel")])

#wave generator
windVelocity = (0.0,0.0,0.0)
#veg_platform_height = 17.2/44.0 + 6.1/20.0
depth = opts.depth
inflowHeightMean = depth
inflowVelocityMean = (0.0,0.0,0,0)
period = opts.peak_period
omega = 2.0*math.pi/period
waveheight = opts.wave_height
amplitude = waveheight/ 2.0
wavelength = opts.peak_wavelength
k = 2.0*math.pi/wavelength
waveDir = numpy.array([1,0,0])
g = numpy.array([0.0,0.0,-9.81])
if opts.wave_type == 'linear':
    waves = WT.MonochromaticWaves(period = period, # Peak period
                                  waveHeight = waveheight, # Height
                                  depth = depth, # Depth
                                  mwl = inflowHeightMean, # Sea water level
                                  waveDir = waveDir, # waveDirection
                                  g = g, # Gravity vector, defines the vertical
                                  waveType="Linear")
elif opts.wave_type == 'Nonlinear':
    waves = WT.MonochromaticWaves(period = period, # Peak period
                                  waveHeight = waveheight, # Height
                                  wavelength = wavelength,
                                  depth = depth, # Depth
                                  mwl = inflowHeightMean, # Sea water level
                                  waveDir = waveDir, # waveDirection
                                  g = g, # Gravity vector, defines the vertical
                                  waveType="Fenton",
                                  Ycoeff = [0.04160592, #Surface elevation Fourier coefficients for non-dimensionalised solution
                                       0.00555874,
                                       0.00065892,
                                       0.00008144,
                                       0.00001078,
                                       0.00000151,
                                       0.00000023,
                                       0.00000007],
                                  Bcoeff = [0.05395079,
                                       0.00357780,
                                       0.00020506,
                                       0.00000719,
                                       -0.00000016,
                                       -0.00000005,
                                       0.00000000,
                                       0.00000000])
elif opts.wave_type == 'single-peaked':
    waves = WT.RandomWaves( Tp = period, # Peak period
                            Hs = waveheight, # Height
                            d = depth, # Depth
                            fp = 1./period, #peak Frequency
                            bandFactor = 2.0, #fmin=fp/Bandfactor, fmax = Bandfactor * fp
                            N = 101, #No of frequencies for signal reconstruction
                            mwl = inflowHeightMean, # Sea water level
                            waveDir = waveDir, # waveDirection
                            g = g, # Gravity vector, defines the vertical
                            gamma=3.3,
                            spec_fun = WT.JONSWAP)
# elif opts.wave_type == 'double-peaked':
#     waves = WT.DoublePeakedRandomWaves( Tp = period, # Peak period
#                                         Hs = waveheight, # Height
#                                         d = depth, # Depth
#                                         fp = 1./period, #peak Frequency
#                                         bandFactor = 2.0, #fmin=fp/Bandfactor, fmax = Bandfactor * fp
#                                         N = 101, #No of frequencies for signal reconstruction
#                                         Hs = 2.0,         #m significant wave height
#                                         d = 2.0,           #m depth
#                                         fp = 1.0/5.0,      #peak  frequency
#                                         bandFactor = 2.0, #controls width of band  around fp
#                                         N = 101,          #number of frequency bins
#                                         mwl = 0.0,        #mean water level
#                                         waveDir = np.array([1,0,0]),
#                                         g = np.array([0, -9.81, 0]),         #accelerationof gravity
#                                         spec_fun = JONSWAP,
#                                         gamma=3.3,                                        mwl = inflowHeightMean, # Sea water level
#                                         waveDir = waveDir, # waveDirection
#                                         g = g, # Gravity vector, defines the vertical
#                                         gamma=10.0,
#                                         spec_fun = WT.JONSWAP,
#                                         Tp_2 = opts.peak_period2)

#  Discretization -- input options
genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=not opts.parallel
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

#for debugging, make the tank short

he = float(wavelength)/30.0#0.0#100
L = (15.0, 1.5, 1.0)

GenerationZoneLength = 1.2
AbsorptionZoneLength= 2.8
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0
#zone 2
xSponge_2 = L[0] - AbsorptionZoneLength
xRelaxCenter_2 = 0.5*(xSponge_2+L[0])
epsFact_solid_2 = AbsorptionZoneLength/2.0

weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False

if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back']
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
                  [L[0],L[1],0.0],#4
                  [xSponge_2,L[1],0.0],#5
                  [xSponge,L[1],0.0],#6
                  [0.0,L[1],0.0]]#7
        
               
        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],            
                     boundaryTags['bottom'],       
                     boundaryTags['bottom']]


        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1],L[2]])
            vertexFlags.append(boundaryTags['top'])

        print vertices
        print vertexFlags

        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,4],
                  [4,5],
                  [5,6],
                  [6,7],
                  [7,0],
                  [1,6],
                  [2,5]]
                 
        segmentFlags=[boundaryTags['front'],
                     boundaryTags['front'],
                     boundaryTags['front'],                   
                     boundaryTags['right'],
                     boundaryTags['back'],
                     boundaryTags['back'],
                     boundaryTags['back'],
                     boundaryTags['left'],
                     boundaryTags['empty'],
                     boundaryTags['empty'] ]
        

        facets=[]
        facetFlags=[]

        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+8,s[0]+8]])
            facetFlags.append(sF)

        bf=[[0,1,6,7],[1,2,5,6],[2,3,4,5]]
        tf=[]
        for i in range(0,3):
         facets.append([bf[i]])
         tf=[ss + 8 for ss in bf[i]]
         facets.append([tf])

        for i in range(0,3):
         facetFlags.append(boundaryTags['bottom'])
         facetFlags.append(boundaryTags['top'])

        print facets
        print facetFlags

        regions=[[xRelaxCenter, 0.5*L[1],0.0],
                 [xRelaxCenter_2, 0.5*L[1], 0.0],
                 [0.5*L[0],0.1*L[1], 0.0]]
        regionFlags=[1,2,3]

        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags,
                                                     )
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

        porosityTypes      = numpy.array([1.0,
                                          1.0,
                                          1.0,
                                          1.0])
        dragAlphaTypes = numpy.array([0.0,
                                      0.0,
                                      0.0,
                                      0.5/1.004e-6])

        dragBetaTypes = numpy.array([0.0,0.0,0.0,0.0])
        
        epsFact_solidTypes = np.array([0.0,0.0,0.0,0.0])

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
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))



# Time stepping
T=40*period
dt_fixed = period/11.0#2.0*0.5/20.0#T/2.0#period/21.0
dt_init = min(0.001*dt_fixed,0.001)
runCFL=0.90
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.0
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
    epsFact_redistance = 1.5
    epsFact_consrv_diffusion = 10.0
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

ns_nl_atol_res = max(1.0e-10,0.0001*he**2)
vof_nl_atol_res = max(1.0e-10,0.0001*he**2)
ls_nl_atol_res = max(1.0e-10,0.0001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.0001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
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
waterLine_x = 2*L[0]
waterLine_z = inflowHeightMean


def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[2]-waterLine_z
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
    return x[2] - inflowHeightMean

#sigma = omega - k*inflowVelocityMean[0]
h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0

#waveData

def waveHeight(x,t):
    return inflowHeightMean + waves.eta(x[0],x[1],x[2],t)
def waveVelocity_u(x,t):
    return waves.u(x[0],x[1],x[2],t,"x")
def waveVelocity_v(x,t):
    return waves.u(x[0],x[1],x[2],t,"y")
def waveVelocity_w(x,t):
    return waves.u(x[0],x[1],x[2],t,"z")

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

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

outflowHeight=inflowHeightMean

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - outflowHeight)

def outflowPhi(x,t):
    return x[2] - outflowHeight

def outflowPressure(x,t):
  if x[2]>inflowHeightMean:
    return (L[2]-x[2])*rho_1*abs(g[2])
  else:
    return (L[2]-inflowHeightMean)*rho_1*abs(g[2])+(inflowHeightMean-x[2])*rho_0*abs(g[2])


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

rzWaveGenerator = RelaxationZoneWaveGenerator(zones={
                                                    # 1:RelaxationZone(xRelaxCenter,
                                                    #                  1.0,
                                                    #                  twpflowVelocity_u,
                                                    #                  twpflowVelocity_v,
                                                    #                  twpflowVelocity_w),
                                                    1:RelaxationZone(xRelaxCenter_2,
                                                                     -1.0, #currently Hs=1-exp_function
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel)})

beam_quadOrder=3
beam_useSparse=False
beamFilename="wavetankBeams"
#nBeamElements=max(nBeamElements,3)

#beam info
beamLocation=[]
beamLength=[]
beamRadius=[]
EI=[]
GJ=[]
lam = 0.25 #0.05 #3.0*2.54/100.0 #57.4e-3
lamx = 3.0**0.5*lam
xs = 1.2
ys = 0.0
xList=[]
yList = []
while xs <= 11.0:
    xList.append(xs)
    xs += lam
while ys<= L[1]:
    yList.append(ys)
    ys+=lamx
for i in xList:
    for j in yList:
        beamLocation.append((i,j))
        beamLength.append(0.415)
        beamRadius.append(0.0032)
        EI.append(3.0e-4) # needs to be fixed
        GJ.append(1.5e-4) # needs to be fixed

xs = 1.2+0.5*lam
ys = 0.5*lamx
xList=[]
yList = []
while xs <= 11.0:
    xList.append(xs)
    xs += lam

while ys<= L[1]:
    yList.append(ys)
    ys+=lamx

for i in xList:
    for j in yList:
        beamLocation.append((i,j))
        beamLength.append(0.415)
        beamRadius.append(0.0032)
        EI.append(3.0e-4) # needs to be fixed
        GJ.append(1.5e-4) # needs to be fixed
nBeamElements = int(beamLength[0]/he*0.5)
nBeamElements=max(nBeamElements,3)
print nBeamElements
