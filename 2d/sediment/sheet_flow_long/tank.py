from proteus import StepControl
from math import *
import proteus.MeshTools
from proteus import Domain, Context
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from proteus.mprans.SedClosure import  HsuSedStress
from proteus.mprans import BoundaryConditions as bc

from proteus.mprans import BodyDynamics as bd

opts=Context.Options([
	# TANK PARAMETERS
    ("Lx", 100.0, "Length of x dim (cm)"),
    ("Ly", 6.5, "Length of y dim (cm)"),
    ("dtout", 0.05, "Time interval for output"),
    #fluid parameters
    ("rho_0", 1.06, "'water' density (g/cm^3)"),
    ("rho_1", 1.06, "'air' density (g/cm^3)"),
    ("nu_0", 3.019, "water kin viscosity (cm^2/s)"),
    ("nu_1", 3.019, "air kin viscosity (cm^2/s)"),
    ('g',np.array([0.0, -980.0, 0.0]),'Gravitational acceleration (cm/s^2)'),
    ('waterLevel', 3.0, 'Water level (cm)'),
    # current
#    ("current",not True, "yes or no"),
    ("inflow_vel", 0.25, "inflow velocity (cm/s)"),
    ("GenZone", not True, "on/off"),
    ("AbsZone", not True, "on/off"),
    # sediment parameters
    ('cSed', 0.6,'Sediment concentration'),
    ('hf',0.5871,'Initial distance between height of sediment and top of tank (cm)'),
    ('rho_s',2.230 ,'sediment material density (g/cm^3)'),
    ('alphaSed', 150.0,'laminar drag coefficient'),
    ('betaSed', 0.0,'turbulent drag coefficient'),
    ('grain',0.11, 'Grain size (cm)'),
    ('packFraction',0.2,'threshold for loose / packed sediment'),
    ('packMargin',0.01,'transition margin for loose / packed sediment'),
    ('maxFraction',0.635,'fraction at max  sediment packing'),
    ('frFraction',0.57,'fraction where contact stresses kick in'),
    ('sigmaC',1.1,'Schmidt coefficient for turbulent diffusion'),
    ('C3e',1.2,'Dissipation coefficient '),
    ('C4e',1.0,'Dissipation coefficient'),
    ('eR', 0.8, 'Collision stress coefficient (module not functional)'),
    ('fContact', 0.05,'Contact stress coefficient'),
    ('mContact', 3.0,'Contact stress coefficient'),
    ('nContact', 5.0,'Contact stress coefficient'),
    ('angFriction', pi/6., 'Angle of friction'),
    ('vos_limiter', 0.62, 'Weak limiter for vos'),
    ('mu_fr_limiter', 1e-3,'Hard limiter for contact stress friction coeff'),
    # numerical options
    ("he", 0.5,"he"),
    ("sedimentDynamics", True, "Enable sediment dynamics module"),
    ("openTop",  True, "Enable open atmosphere for air phase on the top"),
    ("cfl", 0.25 ,"Target cfl"),
    ("duration", 5.0 ,"Duration of the simulation"),
    ("PSTAB", 1.0, "Affects subgrid error"),
    ("res", 1.0e-10, "Residual tolerance"),
    ("epsFact_density", 3.0, "Control width of water/air transition zone"),
    ("epsFact_consrv_diffusion", 1.0, "Affects smoothing diffusion in mass conservation"),
    ("vos_SC",0.9,"vos shock capturing"),
    ("useRANS", 0, "Switch ON turbulence models: 0-None, 1-K-Epsilon, 2-K-Omega1998, 3-K-Omega1988"), # ns_closure: 1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("sigma_e", 1.0, "sigma_e coefficient for the turbulence model"),
    ("K", 0.41, "von Karman coefficient for the turbulence model"),
    ("B", 5.57, "Wall coefficient for the turbulence model"),
    ("Cmu", 0.09, "Cmu coefficient for the turbulence model"),
    ])

# SO Models

VOS_model = 0
VOF_model = 1
LS_model = 2
RD_model = 3
MCORR_model =4
SED_model =5
V_model =6
DP_model = 7
P_model = 8

if opts.useRANS:
    K_model = 9
    EPS_model = 10
    PI_model = 11
else:
    K_model = None
    EPS_model = None   
    PI_model = 9
    
# ----- Sediment stress ----- #

sedClosure = HsuSedStress(opts.alphaSed,
                          opts.betaSed,
                          opts.grain,
                          opts.packFraction,
                          opts.packMargin,
                          opts.maxFraction,
                          opts.frFraction,
                          opts.sigmaC,
                          opts.C3e,
                          opts.C4e,
                          opts.eR,
                          opts.fContact,
                          opts.mContact,
                          opts.nContact,
                          opts.angFriction,
                          opts.vos_limiter,
                          opts.mu_fr_limiter,
                          )

# ----- Phisical constants ----- #
  
# Water
rho_0 = opts.rho_0
nu_0 = opts.nu_0

# Air
rho_1 = opts.rho_1
nu_1 = opts.nu_1

# Sediment
rho_s = opts.rho_s
nu_s = 1000000.0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Gravity
g = opts.g
gamma_0 = abs(g[1])*rho_0

# Initial condition
waterLevel = opts.waterLevel  
waterLine_z= opts.waterLevel 
sediment_level = opts.Ly-opts.hf 

################################################
# Domain and mesh
################################################

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
nd = 2

he = opts.he
L = (opts.Lx, opts.Ly)
dim = dimx, dimy = L
coords = [ dimx/2., dimy/2. ]

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                           }

boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                       }

tank = st.Rectangle(domain, dim=dim, coords=coords)

###################################
# ----- BOUNDARY CONDITIONS ----- #
###################################

tank.BC['y-'].setFreeSlip()
tank.BC['y+'].setNoSlip()

steady_current = wt.SteadyCurrent(U=np.array([opts.inflow_vel,0.0,0.0]),
                                 mwl=opts.waterLevel,
                                 rampTime=0.8)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=steady_current,
                                               smoothing = 3.0*he,
                                               vert_axis=1)
#tank.BC['x-'].pInit_advective.setConstantBC(0.0)
#tank.BC['x-'].pInit_diffusive.setConstantBC(0.0)
#tank.BC['x-'].pInc_diffusive.setConstantBC(0.0)
#tank.BC['x-'].pInc_advective.uOfXT = lambda x,t: -opts.inflow_vel
#tank.BC['x-'].p_advective.setConstantBC(0.0)

tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=opts.waterLevel,
                                                    rhoUp=rho_1,
                                                    rhoDown = rho_0,
                                                    g=g,
                                                    refLevel= opts.Ly,
                                                    smoothing = 3.0*he,
                                                    orientation=np.array([1.0,0.0,0.0]))
#tank.BC['x+'].u_dirichlet.uOfXT = None
#tank.BC['x+'].v_dirichlet.uOfXT = None
#tank.BC['x+'].u_advective.setConstantBC(0.0)
#tank.BC['x+'].v_advective.setConstantBC(0.0)
#tank.BC['x+'].u_diffusive.setConstantBC(0.0)
#tank.BC['x+'].v_diffusive.setConstantBC(0.0)
#tank.BC['x+'].pInc_dirichlet.setConstantBC(0.0)

###############################################
# Turbulence
###############################################

if opts.useRANS:
    kInflow = 1e-6
    dissipationInflow = 1e-6
    tank.BC['x-'].setTurbulentZeroGradient()
    tank.BC['x+'].setTurbulentZeroGradient()
    tank.BC['y-'].setTurbulentZeroGradient()
    tank.BC['y+'].setTurbulentZeroGradient()

###############################################
# Gauges and probes #
###############################################

T=opts.duration
PG = []
gauge_dy=0.01
tank_dim_y=opts.Ly#tank_dim_y
nprobes=int(tank_dim_y/gauge_dy)+1
probes=np.linspace(0., tank_dim_y, nprobes)

for i in probes:
    PG.append((opts.Lx/2., i, 0.),)
v_output = ga.PointGauges(gauges=((('u',), PG),
                                  (('v',), PG),),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='fluidVelocityGauges.csv')
vs_output = ga.PointGauges(gauges=((('us',),PG), 
                                   (('vs',),PG),),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='solidVelocityGauges.csv')
vos_output = ga.PointGauges(gauges=((('vos',),PG),),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='solidFractionGauges.csv')

##########################################
# Numerical Options and other parameters #
##########################################
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

T=opts.duration
weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = opts.dtout
dt_init = min(0.1*dt_fixed,0.001)
nDTout= int(round(T/dt_fixed))
runCFL = opts.cfl

sedimentDynamics=opts.sedimentDynamics

#----------------------------------------------------
#  Discretization -- input options
#----------------------------------------------------

genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'#'vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = 1
pspaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
useRANS = opts.useRANS  # 0 -- None
                        # 1 -- K-Epsilon
                        # 2 -- K-Omega


KILL_PRESSURE_TERM = False
fixNullSpace_PresInc = True
INTEGRATE_BY_PARTS_DIV_U_PresInc = True
CORRECT_VELOCITY = True
STABILIZATION_TYPE = 0 #0: SUPG, 1: EV via weak residual, 2: EV via strong residual

# Input checks
if spaceOrder not in [1, 2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES)
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
    sys.exit()

#  Discretization
nd = tank.nd

if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 4)

if pspaceOrder == 1:
    if useHex:
        pbasis = C0_AffineLinearOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineLinearOnSimplexWithNodalBasis
elif pspaceOrder == 2:
    if useHex:
        pbasis = C0_AffineLagrangeOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineQuadraticOnSimplexWithNodalBasis


###############################################
# Numerical parameters
###############################################

ns_forceStrongDirichlet = False
ns_sed_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01

if useMetrics:
    ns_shockCapturingFactor = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.5
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor =  opts.vos_SC # <------------------------------------- 
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.5
    rd_lag_shockCapturing = False
    epsFact_vos =opts.epsFact_density
    epsFact_density = opts.epsFact_density # 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = opts.epsFact_consrv_diffusion # 0.1
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.9
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor = .9
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = opts.epsFact_density # 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_vos = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = opts.epsFact_consrv_diffusion # 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5

ns_nl_atol_res = max(opts.res, 0.001 * he ** 2)
ns_sed_nl_atol_res = max(opts.res, 0.001 * he ** 2)
vof_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
vos_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
ls_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
rd_nl_atol_res =  max(opts.res, 0.005 * he)
mcorr_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)
phi_nl_atol_res = max(opts.res, 0.001 * he ** 2)
pressure_nl_atol_res = max(opts.res, 0.001 * he ** 2)


########################
# Turbulence
###########################

ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4


#####################################################################
# Functions for model variables - Initial conditions
#####################################################################
hf=opts.hf
x0=6.0
l=opts.Lx/2-x0
x1=x0+2*l
vcut = opts.Ly-hf

def signedDistance(x):
    phi_z = x[1] - 6.0
    return phi_z

def vos_signedDistance(x):
    if abs(x[0]-50)<l and x[1]<vcut:
        phi_z = -min(x[0]-x0,x1-x[0],vcut-x[1])
    elif abs(x[0]-50)<l:
        phi_z = x[1]-vcut
    elif x[1] < vcut:
        if x[0]<x0:
            phi_z = x0-x[0]
        else:
            phi_z = x[0]-x1
    elif x[0]<x0:
        phi_z = sqrt((x[0]-x0)**2 + (x[1]-vcut)**2)
    else:
        phi_z = sqrt((x[0]-x1)**2 + (x[1]-vcut)**2)
    return phi_z

class Suspension_class:
    def __init__(self):
        pass
    def uOfXT(self, x, t=0):
        phi = vos_signedDistance(x)
        smoothing = (epsFact_consrv_heaviside)*he/2.
        Heav = smoothedHeaviside(smoothing, phi)      
        if phi <= -smoothing:
            return opts.cSed
        elif -smoothing < phi < smoothing:
            return opts.cSed * (1.-Heav)            
        else:
            return 1e-10    

Suspension = Suspension_class()

vos_function = Suspension.uOfXT
