from math import *
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
from proteus.WaveTools import MonochromaticWaves
from proteus.mprans import SpatialTools as st
from proteus import MeshTools, AuxiliaryVariables
import proteus.MeshTools
import ode
import numpy as np


# Wave generator
windVelocity = np.array([0.0, 0.0, 0.0])
inflowHeightMean = 1.0
period = 1.94
omega = 2.0*math.pi/period
wavelength = 7.42 #=5.0 without the current
k = 2.0*math.pi/wavelength
rampTime = 2.0*period
meanFrameVelocity = 2.825 #calculated from FFT
outflowHeightMean = inflowHeightMean
netcurrentVelocity = wavelength/period-meanFrameVelocity

Y = [0.00668842,  # Surface elevation Fourier coefficients for non-dimensionalised solution, calculated from FFT
     0.00008619, 
     0.00000106,
     0.00000001]

B = [0.00805507,   # Velocities Fourier coefficients for non-dimensionalised solution, calculated from FFT
     0.00004774, 
     0.00000019,
     0.00000000]

waveheight = np.array(Y)
waterLevel = inflowHeightMean
waveDir = np.array([1,0,0])
g = np.array([0, -9.81, 0])
mwl = inflowHeightMean

waves = MonochromaticWaves(period = period,
                           waveHeight = waveheight,
                           mwl = mwl,
                           depth = waterLevel,
                           g = g,
                           waveDir = waveDir,
                           wavelength = wavelength,
                           waveType = "Fenton",
                           Ycoeff = Y,
                           Bcoeff = B,
                           meanVelocity = ([0.,0.,0.]),
                           phi0 = 0.)



#  Discretization -- input options  
genMesh = True
movingDomain = False
applyRedistancing = True
checkMass=False
freezeLevelSet=False
useOldPETSc = False
useSuperlu = False
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
L = [float(7.0*wavelength),1.50]
he = wavelength/200 # Background refinement
domain.MeshOptions.elementSize(he)

# Refinement parameters
#x_refine = (1.0 , 1.5 , 3.0) #end of zone 1, end of zone 2, end of zone 3 (zone 4 is up to the right wall)
#refinementLevel = (4 , 2) #refinemnt level for zone 1 and zone 2 and 4 respectively zone 3 has the basic refinement level

#Left boundary imposed velocity
leftSponge = 7.42    # generationzonelength=wavelength
rightSponge = 14.84  # absorptionzonelength=wavelength*2.0

"""GenerationZoneLength = wavelength
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0"""


#zone 2
#xSponge_2 = L[0]-AbsorptionZoneLength
#xRelaxCenter_2 = 0.5*(xSponge_2+L[0])
#epsFact_solid_2 = AbsorptionZoneLength/2.0


tank = st.Tank2D(domain, L, leftSponge=leftSponge, rightSponge=rightSponge)
left = right = False
if leftSponge is not None: left = True
if rightSponge is not None: right = True


tank.setAbsorptionZones(left=left, right=right,
                        dragAlphaTypes=0.5/1.005e-6,
                        dragBetaTypes=0.,
                        porosityTypes=1.)

# I can remove completely the arguments if I am using the deault values

"""genzone = tank.setGenerationZones(indice=[1],
                                  epsFact_solid=xSponge/2.0,
                                  waves=waves,
                                  wind=0.,
                                  dragAlphaTypes=0.5/1.005e-6,
                                  dragBetaTypes=0.,
                                  porosityTypes=1.)"""


weak_bc_penalty_constant = 100.0
nLevels = 1
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

st.assembleDomain(domain)
domain.writePoly("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes = False

quad_order = 3


# Boundary conditions and other flags
openTop=False
openSides=False
openEnd=True
smoothBottom=False
smoothObstacle=False



# Gauges
gauge_dx=0.37
PGL=[]
LGL=[]
for i in range(0,int(L[0]/gauge_dx)): #+1 only if gauge_dx is an exact 
  PGL.append([gauge_dx*i,0.5,0])
  LGL.append([(gauge_dx*i,0.0,0),(gauge_dx*i,L[1],0)])
 

gaugeLocations=tuple(map(tuple,PGL)) 
columnLines=tuple(map(tuple,LGL)) 


pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
                                (('p',),    gaugeLocations)),
                  activeTime = (0, 77.60),
                  sampleRate = 0.097,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt')


fields = ('vof',)

columnGauge = LineIntegralGauges(gauges=((fields, columnLines),),
                                 fileName='column_gauge.csv')


#lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_y=0':((0.0,0.0,0.0),(L[0],0.0,0.0))},linePoints=24)

#lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)



#domain.auxiliaryVariables += [pointGauges, columnGauge] 

    
        
        
"""porosityTypes      = numpy.array([1.0,
                                  1.0,
                                  1.0])
dragAlphaTypes = numpy.array([0.0,
                              0.5/1.004e-6,                                     
                              0.0])
dragBetaTypes = numpy.array([0.0,0.0,0.0])

epsFact_solidTypes = np.array([0.0,epsFact_solid_2,0.0])"""


        
# Time stepping
T=40.0 * period
dt_fixed = T
dt_init = min(0.1*dt_fixed,0.1)
runCFL=0.9
nDTout = int(round(T/dt_fixed))


# Numerical parameters
ns_forceStrongDirichlet = False #True
backgroundDiffusionFactor = 0.0
if useMetrics:
    ns_shockCapturingFactor  = 0.25
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.35
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.35
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.75
    rd_lag_shockCapturing = False
    epsFact_density = 3.0
    epsFact_viscosity = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
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
    ns_shockCapturingFactor = 0.9
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


ns_nl_atol_res = max(1.0e-10,0.00001*he**2)
vof_nl_atol_res = max(1.0e-10,0.00001*he**2)
ls_nl_atol_res = max(1.0e-10,0.0001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.0001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)


# Turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4

    
# Physical parameters                                                     
rho_0 = 998.2
nu_0  = 1.004e-6

rho_1 = 1.205
nu_1  = 1.500e-5 

sigma_01 = 0.0
                                                      


# Initial condition
waterLine_x = 2*L[0]
waterLine_z = inflowHeightMean



# Boundary Conditions
tank.BC.top.setOpenAir()
tank.BC.bottom.setFreeSlip()

tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waves, vert_axis=1, windSpeed=windVelocity, air=1., water=0., smooth=False)

#tank.BC.right.hydrostaticPressureOutletWithDepth(seaLevel=outflowHeightMean, rhoUp=rho_1, rhoDown=rho_0, g=g, refLevel=, pRef=0.0, vert_axis=1, air=1.0, water=0.0)

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

        
# Solution variables
def ramp(t):
 if t<rampTime:
   return 1#/rampTime*t
 else:
   return 1

def theta(x,t):
    return k*x[0] - omega*t + pi/2.0

def z(x):
    return x[1] - inflowHeightMean

h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0
    

def waveHeight(x,t):
   waterDepth = inflowHeightMean 
   for i in range(0,int(len(Y))):  
       waterDepth += Y[i]*cos((i+1)*theta(x,t))/k
   return waterDepth*ramp(t)
  
       
def waveVelocity_u(x,t):
   wu = wavelength/period-meanFrameVelocity
   for i in range(0,int(len(B))): 
     wu += sqrt(abs(g[1])/k)*(i+1)*B[i]*cosh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*cos((i+1)*theta(x,t))
    
   return wu*ramp(t)


def waveVelocity_v(x,t):
   wv=0
   for i in range(0,int(len(B))): 
     wv += sqrt(abs(g[1])/k)*(i+1)*B[i]*sinh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*sin((i+1)*theta(x,t)) 

   return wv*ramp(t)


def wavePhi(x,t):
    return x[1] - inflowHeightMean


def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def twpflowVelocity_u(x,t):
    waterspeed =  waveVelocity_u(x,t)
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


def outflowPressure(x,t):
  if x[1]>outflowHeightMean:
    return (L[1]-x[1])*rho_1*abs(g[1])
  else:
    return (L[1]-outflowHeightMean)*rho_1*abs(g[1])+(outflowHeightMean-x[1])*rho_0*abs(g[1])


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
    waterspeed = netcurrentVelocity
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,outflowPhi(x,t)-epsFact_consrv_heaviside*he)
    u = (1.0-H)*waterspeed
    return u
 

def zeroVel(x,t):
    return 0.0
 




