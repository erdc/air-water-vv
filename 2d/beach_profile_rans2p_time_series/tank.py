from proteus import Domain, Context
#from proteus.mprans
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import proteus.MeshTools
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 1.425, "Height of free surface above bottom"), # Choose waterLevels=[0.425, 0.463]
    # waves
    ('waveType', 'TimeSeries', 'Wavetype for regular waves,  Linear or Fenton for regular waves, TimeSeries for importing field data'),
    ("wave_period", 10, "Period of the waves"), # Choose periods=[0.5, 0.8, 1.1, 1.5, 2.8, 3.9, 4.0]
    ("wave_height", 0.234, "Height of the waves"), # # Choose for d=0.425-->[0.025, 0.075, 0.125, 0.234]. Choose for d=0.463-->[0.025, 0.075, 0.125, 0.254].
    ('wavelength', 9, 'Wavelength only if Fenton is activated'), # Choose for d=0.425-->[0.4, 1.0, 1.8, 2.9]. Choose for d=0.463-->[0.4, 1.0, 1.8, 2.9, 3.0, 5.7, 5.9, 8.8, 9.4].
    # Geometry of the tank - left lower boundary at (0.,0.,0.)
    ("Ls",   2.0, "Distance of the front toe of the structure end from generation zone in wavelengths"),
    ("Lend", 2.0, "Distance of the back toe of the structure end from absorption zone in wavelengths"),
    ("Lgen", 1., "Length of generation zone in wavelegths"),
    ("Labs", 2., "Length of absorption zone in wavelegths"),
    ("h", 1.0, "Height of domain in meters"),
    ("use_gmsh",False, "use_gmsh"),
    ("he", 0.01, "he"),# this is not the he used
    
    # numerical options
    ("refinement_level", 10. ,"he=walength/refinement_level"),
    ("cfl", 0.9 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing"),
    ('movingDomain', not True, "Moving domain and mesh option"),
    ('conservativeFlux', not True,'Fix post-processing velocity bug for porous interface'),
    ])


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
period=opts.wave_period

waterLevel=opts.water_level

waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface

waveHeight=opts.wave_height

inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])


# ----- Phisical constants ----- #

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))




# ----- WAVE input ----- #

if opts.waveType=='Fenton' or opts.waveType=='Linear' :
    waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  Nf=8,
                                  Ycoeff=np.zeros(8),
                                  Bcoeff=np.zeros(8),
                                  waveDir=waveDir,
                                  wavelength=None, # if wave is linear I can use None
                                  waveType=opts.waveType)


if opts.waveType=='TimeSeries':
    waveinput = wt.TimeSeries(timeSeriesFile = "TimeSeries.txt",
                              skiprows = 0,
                              timeSeriesPosition = np.array([0,0,0]),
                              depth = waterLevel,
                              N = 32,
                              mwl = waterLevel,
                              waveDir = np.array([1,0,0]),
                              g = np.array([0,-9.81,0]),
                              Lgen = np.array([10,0,0]),
                              )
                              


#---------Domain Dimension
nd = 2
wl = waveinput.wavelength

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = opts.Lgen*wl
#L_rightSpo = opts.Labs*wl





boundaryOrientations = {'bottom': np.array([0., -1.,0.]),
                        'right': np.array([1., 0.,0.]),
                        'top': np.array([0., 1.,0.]),
                        'left': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'beach_profile': None,
                       }
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,
                'sponge': 5,
                'beach_profile': 6,
               }



##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################


data = np.loadtxt("profile_data.txt")
data = data.transpose()


vertices = []

for i in range(0, len(data)):
    vertices.append((data[i,0]+L_leftSpo,data[i,1]))
    
    

vertices.append((max(data[:,0])+L_leftSpo,min(data[:,1])))

vertices.append((min(data[:,0])+L_leftSpo,max(data[:,1])))

vertices.append((min(data[:,0]),max(data[:,1])))

vertices.append((min(data[:,0]),min(data[:,1])))


vertexFlags = np.array(np.ones(len(data)))
vertexFlags = vertexFlags *6

vertexFlags = np.append(vertexFlags, [1, 3, 3, 1])

segments = []

for i in range(0, len(data)-1):
    segments.append((i,i+1))


segments.append((1295,1296))
segments.append((1296,0))
segments.append((1295,1297))
segments.append((1297,0))
segments.append((1297,1298))
segments.append((1298,1299))
segments.append((1299,0))


segmentFlags = np.array(np.ones(len(data)-1))
segmentFlags = segmentFlags *6

segmentFlags = np.append(segmentFlags, [2, 1, 3, 5, 3, 4, 1])

regions = [ [ 0.1 , 0.10 ],
            [ 10 , 1 ],
            [ 130 , 0.1 ]
            ]

regionFlags=np.array([1, 2, 3])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)





#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################



tank.BC['top'].setAtmosphere()
tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, smoothing=3.0*opts.he)
#tank.BC.bottom.setFreeSlip()
tank.BC['bottom'].setNoSlip()
#tank.BC.right.setFreeSlip()
tank.BC['right'].setNoSlip()
tank.BC['sponge'].setNonMaterial()
#tank.BC.moving_sponge.setNonMaterial()
tank.BC['beach_profile'].setFreeSlip()



########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################

omega = 2*np.pi/opts.wave_period

tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, smoothing= 3*opts.he, dragAlpha=10.*omega/nu_0)
                        


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = 7.*period

PG=[]
LG=[]
LG1=[]
LG2=[]
LG3=[]
LG4=[]



######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

#he = waveinput.wavelength/opts.refinement_level
he = 0.1
domain.MeshOptions.he = he 
domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.use_gmsh = opts.use_gmsh



from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

st.assembleDomain(domain)

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = 0.1
dt_init = min(0.1*dt_fixed,0.001)
T = T
nDTout= int(round(T/dt_fixed))
runCFL = opts.cfl

#----------------------------------------------------
#  Discretization -- input options
#----------------------------------------------------

checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=opts.freezeLevelSet
useOnlyVF = False # if TRUE  proteus uses only these modules --> twp_navier_stokes_p + twp_navier_stokes_n
                  #                                              vof_p + vof_n
movingDomain=opts.movingDomain
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988

genMesh= True

# By DEFAULT on the other files.py -->  fullNewtonFlag = True 
#                                       multilevelNonlinearSolver & levelNonlinearSolver == NonlinearSolvers.Newton

useOldPETSc=False # if TRUE  --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.PETSc
                  # if FALSE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.KSP_petsc4py

useSuperlu = False #if TRUE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.LU

spaceOrder = 1
useHex     = False # used for discretization, if 1.0 --> CubeGaussQuadrature
                   #                          ELSE   --> SimplexGaussQuadrature

useRBLES   = 0.0 # multiplied with subGridError
useMetrics = 1.0 # if 1.0 --> use of user's parameters as (ns_shockCapturingFactor, ns_lag_shockCapturing, ecc ...)
useVF = opts.useVF # used in the smoothing functions as (1.0-useVF)*smoothedHeaviside(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf))


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
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5 # magnifies numerical viscosity in NS (smoothening velocity fields)
    ns_lag_shockCapturing = True # lagging numerical viscosity speedsup Newton but destabilzes the solution
    ns_lag_subgridError = True # less nonlinear but less stable
    ls_shockCapturingFactor  = 0.5 # numerical diffusion of level set (smoothening phi)
    ls_lag_shockCapturing = True # less nonlinear but less stable
    ls_sc_uref  = 1.0 # reference gradient in numerical solution (higher=more diffusion)
    ls_sc_beta  = 1.5 # 1 is fully nonlinear, 2 is linear
    vof_shockCapturingFactor = 0.5 # numerical diffusion of level set (smoothening volume of fraction)
    vof_lag_shockCapturing = True # less nonlinear but less stable
    vof_sc_uref = 1.0 
    vof_sc_beta = 1.5 
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0 # control width of water/air transition zone
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0 # affects smoothing diffusion in mass conservation
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True # False
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

ns_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = max(1.0e-12,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

# Initial condition
waterLine_x = 2*max(data[:,0])
waterLine_z = waterLevel

tank_dim = [max(data[:,0]), max(data[:,1])]


def waveHeight(x,t):
    waterDepth = waveinput.eta(x, t) + waveinput.mwl
    return waterDepth


def wavePhi(x,t):
    [nd-1]- waveHeight(x,t)
    

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[nd-1]-waterLine_z

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
