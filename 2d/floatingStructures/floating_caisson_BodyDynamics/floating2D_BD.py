from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd

opts=Context.Options([
    # predefined test cases
    ("water_level", 0.900, "Height of free surface above bottom"),
    # Geometry
    ("tank_dim", (1., 1.2,), "Dimensions of the tank"),
    ("tank_sponge", (2., 2.), "Length of relaxation zones (front/back, left/right)"),
    ("tank_BC", 'FreeSlip', "tank boundary conditions: NoSlip or FreeSlip"),
    # waves
    ('wave', False, 'Enable  generation'),
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 1.30, "Period of the waves"),
    ("wave_height", 0.167, "Height of the waves"),
    ('wavelength', 2.121, 'Wavelength only if Fenton is activated'),
    ('Ycoeff', [0.21107604, 0.07318902, 0.02782228, 0.01234846, 0.00618291, 0.00346483, 0.00227917, 0.00194241], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [0.23112932, 0.03504843, 0.00431442, 0.00036993, 0.00004245, 0.00001877, 0.00000776, 0.00000196], 'Bcoeff only if Fenton is activated'),
    ('Nf', 8 ,'Number of frequency components for fenton waves'),
    ('meanVelocity', [ 0., 0., 0.],'Velocity used for currents'),
    ('phi0', 0.0 ,'Initial phase for waves'),
    ('Uwind', [0.0, 0.0, 0.0], 'Set air velocity'),
    ('fast', True ,'Switches ON fast cosh approximation'),
   # caisson
    ("caisson2D", True, "Switch on/off caisson2D"),
    ("center", (0.5, 0.9),"Coord of the caisson center"),
    ('dimx', 0.3, 'X-dimension of the caisson2D'),
    ('dimy', 0.1, 'Y-dimension of the caisson2D'),
    ('width', 0.9, 'Z-dimension of the caisson2D'),
    ('mass', 125., 'Mass of the caisson2D [kg]'),
    ('caisson_BC', 'FreeSlip', 'caisson2D boundaries: NoSlip or FreeSlip'),
    ("free_x", (0., 0., 0.), "Translational DOFs"),
    ("free_r", (0., 0., 1.0), "Rotational DOFs"),
    ("caisson_inertia", 0.236, "Inertia of the caisson [kg m2]"),
    ("rotation_angle", 15., "Initial rotation angle (in degrees)"),
    ("Tn", 0.93, "Roll natural period"),
    ("overturning", True, "Switch on/off overturning module"),
    ('scheme', 'Forward_Euler', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Forward_Euler)'),
    # numerical options
    ("refinement_level", 0.0,"he=walength/refinement_level"),
    ("he", 0.04,"he=walength/refinement_level"),
    ("cfl", 0.4,"Target cfl"),
    ("duration", 10., "Durarion of the simulation"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', True, "Moving domain and mesh option"),
    ('conservativeFlux', True,'Fix post-processing velocity bug for porous interface'),
    ])


# --- DOMAIN
domain = Domain.PlanarStraightLineGraphDomain()



# --- WAVE CONDITIONS
period=opts.wave_period

waterLevel=opts.water_level

waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface

waveHeight=opts.wave_height

inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])


# --- Phisical constants
rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))


# --- WAVE input
if opts.wave == True:
    waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=opts.wavelength,       # used by fenton waves
                                  waveType=opts.waveType, 
                                  Ycoeff=np.array(opts.Ycoeff),     # used by fenton waves
                                  Bcoeff=np.array(opts.Bcoeff),     # used by fenton waves
                                  Nf=opts.Nf,                       # used by fenton waves
                                  meanVelocity = np.array(opts.meanVelocity),
                                  phi0 = opts.phi0,
                                  fast = opts.fast,)

# --- Domain Dimension
nd = 2
if opts.wave == True:
    wl = waveinput.wavelength
else:
    wl = 1.0

# --- MESH SIZE
if opts.he == 0.0:
    he = wl/opts.refinement_level
else:
    he = opts.he 


#######################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = opts.tank_sponge[0]*wl
L_rightSpo = opts.tank_sponge[1]*wl

# --- caisson2D
dimx = opts.dimx
dimy = opts.dimy
ctr = opts.center

xc1 = ctr[0]-0.5*dimx
yc1 = ctr[1]-0.5*dimy
xc2 = xc1+dimx
yc2 = yc1
xc3 = xc2
yc3 = yc2+dimy
xc4 = xc1
yc4 = yc3

tank_dim = opts.tank_dim
xEnd = tank_dim[0]+L_rightSpo

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
# caisson2D
############################################################################################################################################################################################################

if opts.caisson2D:
    dimx = dimx
    dimy = dimy
    dim = (dimx,dimy)
    coords = ctr
    VCG = dim[1]/2.
    width = opts.width
    mass = opts.mass
    volume = float(dimx*dimy*width)
    density = float(mass/volume)
    inertia = opts.caisson_inertia/mass/width

# --- Shape properties setup
    caisson = st.Rectangle(domain, dim=dim, coords=coords)
    xc1, yc1 = caisson.vertices[0][0], caisson.vertices[0][1]
    xc2, yc2 = caisson.vertices[1][0], caisson.vertices[1][1]
   

# --- Body properties setup
    caisson2D = bd.RigidBody(shape=caisson, substeps=20)
    free_x = np.array(opts.free_x)
    free_r = np.array(opts.free_r)
    caisson2D.setMass(mass)
    caisson2D.setConstraints(free_x=free_x, free_r=free_r)
    rotation = np.radians(opts.rotation_angle)
    caisson.rotate(rotation)
    caisson2D.It = inertia
    caisson2D.setNumericalScheme(scheme=opts.scheme)
    caisson2D.setRecordValues(filename='caisson2D', all_values=True)

##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

tank = st.Tank2D(domain, tank_dim)

#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

# --- Caisson2D
for bc in caisson.BC_list:
    if opts.caisson_BC == 'NoSlip':
        bc.setNoSlip()
    if opts.caisson_BC == 'FreeSlip':
        bc.setFreeSlip()
# --- Tank
tank.BC['y+'].setAtmosphere()
if opts.tank_BC == 'NoSlip':
    tank.BC['y-'].setNoSlip()
if opts.tank_BC == 'FreeSlip':
    tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes

########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################

tank_sponge = opts.tank_sponge
if opts.wave is True:
    dragAlpha = 5*(2*np.pi/period)
else:
    dragAlpha = 5*(2*np.pi/opts.Tn)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if left:
    if opts.wave is True:
        smoothing = opts.he*3.
        tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=dragAlpha)
        tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
    else:
        tank.setAbsorptionZones(x_n=left, dragAlpha=dragAlpha)
        tank.BC['x-'].setNoSlip()
else:
    tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right, dragAlpha=dragAlpha)
if opts.caisson2D:
    tank.setChildShape(caisson, 0)

############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = opts.duration

gauge_dx=0.25
tank_dim_x=int(tank_dim[0])
nprobes=int(tank_dim_x/gauge_dx)+1
probes=np.linspace(0., tank_dim_x, nprobes)
PG=[]
if opts.caisson2D:
    zProbes=yc1*0.5
else:
    zProbes=opts.water_level*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

if opts.caisson2D:
    gauge_dy=0.01
    tol=np.array([1*(10**-5),1*(10**-5),0.])
    i_point_f=np.array([caisson.vertices[0][0],caisson.vertices[0][1],0.])
    i_point_f += -tol #to avoid floating point error
    i_point_b=np.array([caisson.vertices[1][0],caisson.vertices[1][1],0.])
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


######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

he = he
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
waterLine_z = waterLevel


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

