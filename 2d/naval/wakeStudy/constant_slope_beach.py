from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd

opts=Context.Options([
    # predefined test cases
    ("water_level", 5., "Height of free surface above bottom"),
    # Geometry
    ("tank_dim", (180., 7.5,), "Dimensions of the tank"),
    ("tank_sponge", (1., 0.), "Length of relaxation zones (front/back, left/right)"),
    ("tank_BC", 'FreeSlip', "tank boundary conditions: NoSlip or FreeSlip"),
    # waves
    ('wave', True, 'Enable  generation'),
    ("period", 3.5 , "wave period"),
    ("wave_height", 0.5, "wave height"), 
    ("wave_dir", np.array([1.,0.,0.]),"Direction of the waves"),
    # numerical options
    ("refinement_level", 0.0,"he=walength/refinement_level"),
    ("he", 0.2,"he=walength/refinement_level"),
    ("cfl", 0.4,"Target cfl"),
    ("duration", 60., "Durarion of the simulation"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', not True, "Moving domain and mesh option"),
    ('conservativeFlux', True,'Fix post-processing velocity bug for porous interface'),
    # obstacle dimensions (tweak to set up the geometry of the sloping beach)
    ("slope_start", 100, "x coordinate where the beach slope is starting within the tank"), # approx. 3 wave lengths away form the generation
    ("beach_slope", 10, "1/slope"),
    ("beach_crest", 7, "elevation of the highest point of the beach")
    ])


# --- DOMAIN
domain = Domain.PlanarStraightLineGraphDomain()


# --- Phisical constants
rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g =np.array([0.,-9.81,0.])
gAbs=sqrt(sum(g**2))
waterLevel = opts.water_level

# --- WAVE input
period = opts.period
height = opts.wave_height
mwl = opts.water_level
depth = opts.water_level
direction = opts.wave_dir
N = 32
Nwaves = 5
overl = 0.7
cutoff = 0.1
#wave = wt.MonochromaticWaves(period, height, mwl, depth, g, direction)

#Class loading
wave1 = wt.TimeSeries(
    timeSeriesFile="Time_series_short.txt",
    skiprows=0,
    timeSeriesPosition=np.array([-1.,0.,0.]),
    depth = opts.water_level,
    N = N,
    mwl=opts.water_level,
    waveDir = opts.wave_dir,
    g =g,
    rec_direct=False,
    window_params = {"Nwaves":Nwaves,"Tm":opts.period,"Window":"costap","Overlap":overl,"Cutoff":cutoff},
    Lgen = np.array([1,0,0])
    )

    
wave2 = wt.TimeSeries(
    timeSeriesFile="Time_series_long.txt",
    skiprows=0,
    timeSeriesPosition=np.array([-1.,0.,0.]),
    depth = opts.water_level,
    N = N,
    mwl=opts.water_level,
    waveDir = opts.wave_dir,
    g =g,
    rec_direct=True,
    Lgen = np.array([1,0,0])
    )
   
wave =  wt.CombineWaves([wave1,wave2])

time = np.loadtxt("Time_series_short.txt")[:,0]
time = time 
data = time.copy()
x0 = np.array([-1.,0.,0.])
for i,t in enumerate(time):
    data[i] = wave.eta(x0,t)

np.savetxt("out_short.txt",zip(time,data))

   



#######################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = opts.tank_sponge[0]
L_rightSpo = opts.tank_sponge[1]


tank_dim = opts.tank_dim


boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'porousLayer': None,
                        'moving_porousLayer': None,
                       }

boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'sponge' : 5,
                'porousLayer' : 6,
                'moving_porousLayer' : 7,
               }


##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################
obstacle = [
           [ [opts.slope_start ,0.],
            [opts.slope_start + (opts.beach_crest*opts.beach_slope),opts.beach_crest],
            [opts.slope_start + (opts.beach_crest*opts.beach_slope),0.]
            ]
            ]
                     

tank = st.TankWithObstacles2D(domain = domain, dim = tank_dim, obstacles = obstacle, hole = True)
#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################


# --- Tank
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

#tank.BC['x-'].setFixedNodes()
#tank.BC['x+'].setFixedNodes()
#tank.BC['y+'].setTank()  # sliding mesh nodes
#tank.BC['y-'].setTank()  #sliding mesh nodes
tank.BC['sponge'].setNonMaterial()
tank.BC['sponge'].setFixedNodes()

########################################################################################################################################################################################################################################################################################################################################################
# -----  ABSORPTION ZONE BEHIND PADDLE  ----- #
########################################################################################################################################################################################################################################################################################################################################################
tank_sponge = opts.tank_sponge
dragAlpha = 5*(2*np.pi/opts.period)/1e-6
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
left = True
smoothing = opts.he*3.
tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=dragAlpha)


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
"""
T = opts.duration

gauge_dx=0.25
tank_dim_x=int(tank_dim[0])
nprobes=int(tank_dim_x/gauge_dx)+1
probes=np.linspace(0., tank_dim_x, nprobes)
PG=[]
if opts.paddle2D:
    zProbes=yc1*0.5
else:
    zProbes=opts.water_level*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

if opts.paddle2D:
    gauge_dy=0.01
    tol=np.array([1*(10**-5),1*(10**-5),0.])
    i_point_f=np.array([paddle.vertices[0][0],paddle.vertices[0][1],0.])
    i_point_f += -tol #to avoid floating point error
    i_point_b=np.array([paddle.vertices[1][0],paddle.vertices[1][1],0.])
    i_point_b += tol #to avoid floating point error
    yProbes = np.linspace(i_point_f[1],i_point_f[1]+dimy, int(dimy/gauge_dy)+1)
    LG1=[]
    LG2=[]
    for j in yProbes:
        LG1.append((i_point_f[0],j,0.),)
        LG2.append((i_point_b[0],j,0.),)

tank.attachPointGauges(
        'ls',
        gauges=((('phi',),PG),),
        activeTime = (0., T),
        sampleRate=0.,
        fileName='levelset_gauges.csv')

"""
######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

he = opts.he
domain.MeshOptions.he = he


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
T = opts.duration
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

genMesh=True

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
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = ecH = epsFact_density
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
waterLine_x = 2*tank_dim[0]
waterLine_z = opts.water_level




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

