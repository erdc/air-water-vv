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
from proteus.mbd import ChRigidBody as crb
from proteus.mbd import pyChronoCore as pych
from proteus.ctransportCoefficients import smoothedHeaviside 
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.mprans import BodyDynamics as bd
from proteus import WaveTools as wt



opts=Context.Options([
    ("dtout", 0.05, "Time interval for output"),
    ("Refiment", 4, "refinement"),
    ("tank_dim_x", 1.6, "x_dim"),
    ("tank_dim_y", 0.6, "y_dim"),
    ("hole_tank", True, "hole"),
    ("waterLevel" , 0.35+0.16, "waterLevel"),
    # current
    ("inflow_vel", 1e-10, "inflow velocity"),
    ("GenZone", not True, "on/off"),
    ("AbsZone", not True, "on/off"),
    # cylinder
    ("cylinder_radius", 0.05, "radius of cylinder"),
    ("cylinder_pos_x", 0.8, "x position of cylinder"),
    ("cylinder_pos_y", 0.085+0.05+0.16, "y position of cylinder"),
    ("circle2D", False, "switch on/off cylinder"),
    ("circleBC", 'NoSlip','circle BC'),
    # sediment parameters
    ('cSed', 0.62,'Sediment concentration'),
    # numerical options
    ("he", 0.04,"he"),
    ("sedimentDynamics", True, "Enable sediment dynamics module"),
    ("openTop",  True, "Enable open atmosphere for air phase on the top"),
    ("cfl", 0.5 ,"Target cfl"),
    ("duration", 1.0 ,"Duration of the simulation"),
    ("PSTAB", 1.0, "Affects subgrid error"),
    ("res", 1.0e-10, "Residual tolerance"),
    ("epsFact_density", 3.0, "Control width of water/air transition zone"),
    ("epsFact_consrv_diffusion", 1.0, "Affects smoothing diffusion in mass conservation"),
    ("vos_SC",0.9,"vos shock capturing"),
    ("useRANS", 1, "Switch ON turbulence models: 0-None, 1-K-Epsilon, 2-K-Omega1998, 3-K-Omega1988"), # ns_closure: 1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("sigma_e", 1.0, "sigma_e coefficient for the turbulence model"),
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

# Domain dimensions

nd = 2





steady_current = wt.SteadyCurrent(U=np.array([opts.inflow_vel,0.,0.]),mwl=opts.waterLevel,rampTime=0.1)
I = 0.03
kInflow = 1.5*abs(opts.inflow_vel*I)
dissipationInflow = 0.09*(kInflow**1.5)*opts.waterLevel


# ----- Sediment stress ----- #

sedClosure = HsuSedStress(aDarcy =  150.0,
                          betaForch =  1.75,
                          grain =  0.0002,
                          packFraction =  0.2,
                          packMargin =  0.01,
                          maxFraction =  0.635,
                          frFraction =  0.57,
                          sigmaC =  1.1,
                          C3e =  1.2,
                          C4e =  1.0,
                          eR =  0.8,
                          fContact =  0.05,
                          mContact =  3.0,
                          nContact =  5.0,
                          angFriction =  pi/6.,
                          vos_limiter = 0.62,
                          mu_fr_limiter = 1e-1,
                          )

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()







# ----- Physical constants ----- #
  
# Water
rho_0 = 998.2
nu_0 = 1.004e-6

# Air
rho_1 =  1.205 #
nu_1 =  1.500e-5 # 

# Sediment

rho_s = 2600.0 # rho_0
nu_s = 1000000.0 # 0.0 # nu_0 # 
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Gravity
g = np.array([0.0, -9.8, 0.0])
gamma_0 = abs(g[1])*rho_0

# Initial condition
#waterLine_x = opts.waterLine_x
#waterLine_z = opts.waterLine_z
waterLevel = opts.waterLevel 

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Domain and mesh
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#L = (opts.Lx, opts.Ly)
#he = L[0]/opts.refinement
he = opts.he
#dim = dimx, dimy = L
#coords = [ dimx/2., dimy/2. ]

############ TANK ###################



tank_dim = (opts.tank_dim_x, opts.tank_dim_y)

x1=0.2




### CYLINDER #####
cylinder_pos = np.array([opts.cylinder_pos_x, opts.cylinder_pos_y, 0.])
cylinder_radius = opts.cylinder_radius





if opts.circle2D:
    circle = st.Circle(domain=domain, radius=cylinder_radius, coords=(cylinder_pos[0],cylinder_pos[1]), barycenter=(cylinder_pos[0],cylinder_pos[1]), nPoints=28)

    circle2D = bd.RigidBody(shape=circle)
    free_x = (0.0,0.0,0.0)
    free_r = (0.0,0.0,0.0)
    circle2D.setConstraints(free_x=free_x,free_r=free_r)
    circle2D.setNumericalScheme(None)
    circle2D.setRecordValues(filename='circle2D',all_values=True)





boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1, 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'circle': None,
                        'hole_x-': np.array([-1., 0.,0.]),
                        'hole_x+': np.array([+1., 0.,0.]),
                        'hole_y-': np.array([0., -1.,0.]),


                           }
boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                    'circle':6,
                    'hole_x-':7,
                    'hole_x+':8,
                    'hole_y-':9,
                       }




if opts.hole_tank:

    vertices=[[0.0, 0.24],#0
              [x1,  0.24],#1
              [0.3,  0.24],#2
              [0.3,  0.0],#3
              [1.3,  0.0],#4
              [1.3,  0.24],#5
              [1.4,  0.24],#6
              [tank_dim[0],0.24],#7
              [tank_dim[0],tank_dim[1]+0.16], #8
              [1.4, tank_dim[1]+0.16], #9
              [0.2, tank_dim[1]+0.16], #10
              [0.0, tank_dim[1]+0.16], #11
              ]

    vertexFlags=np.array([1, 1, 1,
                          9, 9,
                          1, 1, 1,
                          3, 3, 3, 3,
                          ])
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

              [11,0],

              [1,10],
              [6,9],
             ]

    segmentFlags=np.array([ 1, 1, 
                            7, 9, 8,
                            1, 1,
                            2,
                            3, 3, 3, 
                            4,
                            5, 5,
                         ])


regions = [ [ 0.90*x1 , 0.50*tank_dim[1] ],
            [ 0.7 , 0.50*tank_dim[1] ],
            [ 0.95*tank_dim[0] , 0.50*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)





#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

if opts.circle2D:
    
    for bc in circle.BC_list:
        if opts.circleBC == 'FreeSlip':
            bc.setFreeSlip()
        if opts.circleBC == 'NoSlip':
            bc.setNoSlip()



tank.BC['y-'].setNoSlip()

################################## y+ ######################
tank.BC['y+'].setFreeSlip()
if opts.openTop:
    tank.BC['y+'].setAtmosphere(orientation=np.array([0,1,0]))

####################################### x- ########################



tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=steady_current, smoothing = 3*he, vert_axis=1, kInflow = kInflow, dInflow = dissipationInflow)
    

    


#################################### x+ ###############################
#tank.BC['x+'].setFreeSlip()


    
tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=opts.waterLevel, rhoUp=rho_1, rhoDown = rho_0, g=g, refLevel= tank_dim[1], smoothing = 3*he, kInflow = kInflow, kInflowAir=kInflow)
    

tank.BC['hole_x+'].setNoSlip()


tank.BC['hole_x-'].setNoSlip()


tank.BC['hole_y-'].setNoSlip()

###########################Wall functions########################

if opts.useRANS:
    # Skin friction and friction velocity for defining initial shear stress at the wall
    from proteus.mprans import BoundaryConditions as bc
    d = opts.he
    LL = opts.tank_dim_y
    U0 = opts.inflow_vel
    Cmu = 0.09
    Re0 = U0*d/nu_0
    cf = 0.045*(Re0**(-1./4.))
    Ut = U0*sqrt(cf/2.)
    kappaP = (Ut**2)/sqrt(Cmu)
    Y_ = opts.he 
    Yplus = Y_*Ut/nu_0
    dissipationP = (Ut**3)/(0.41*Y_)
    
    # Pipeline initial conditions, Schlichting
    
    kWallBottom = bc.kWall(Y=Y_, Yplus=Yplus, b_or=boundaryOrientations['y-'], nu=nu_0)
    kWallhole1 = bc.kWall(Y=Y_, Yplus=Yplus, b_or=boundaryOrientations['hole_x-'], nu=nu_0)
    kWallhole2 = bc.kWall(Y=Y_, Yplus=Yplus, b_or=boundaryOrientations['hole_x+'], nu=nu_0)
    kWallhole3 = bc.kWall(Y=Y_, Yplus=Yplus, b_or=boundaryOrientations['hole_y-'], nu=nu_0)
    kWalls = [kWallBottom, kWallhole1,kWallhole2,kWallhole3]
    
    # Attached to 'twp' in auxiliary variables
    U0 = U0 * np.array([1,0,0])
    wallBottom = bc.WallFunctions(turbModel='ke', kWall=kWallBottom, b_or=boundaryOrientations['y-'], Y=Y_, Yplus=Yplus, U0=U0, nu=nu_0)
    wallhole1 = bc.WallFunctions(turbModel='ke', kWall=kWallhole1, b_or=boundaryOrientations['hole_x-'], Y=Y_, Yplus=Yplus, U0=U0, nu=nu_0 )
    wallhole2 = bc.WallFunctions(turbModel='ke', kWall=kWallhole2, b_or=boundaryOrientations['hole_x+'], Y=Y_, Yplus=Yplus, U0=U0, nu=nu_0)
    wallhole3 = bc.WallFunctions(turbModel='ke', kWall=kWallhole3, b_or=boundaryOrientations['hole_y-'], Y=Y_, Yplus=Yplus, U0=U0, nu=nu_0)

    walls = [wallBottom, wallhole1,wallhole2,wallhole3]

    tank.setTurbulentWall(walls)
    tank.setTurbulentKWall(kWalls)
    tank.BC['y-'].setWallFunction(walls[0])    
    tank.BC['hole_x-'].setWallFunction(walls[1])
    tank.BC['hole_x+'].setWallFunction(walls[2])
    tank.BC['hole_y-'].setWallFunction(walls[3])



######################################################################################################################################################################################################################
# Gauges and probes #
######################################################################################################################################################################################################################

T=opts.duration
PG = []
gauge_dy=0.01
tank_dim_y=tank_dim[1]
nprobes=int(tank_dim_y/gauge_dy)+1
probes=np.linspace(0., tank_dim_y, nprobes)
for i in probes:
    PG.append((tank_dim[0]/2., i, 0.),)
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


######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################
 
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
dt_init = 1e-5
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
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
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


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Numerical parameters
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

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
    vos_shockCapturingFactor = opts.vos_SC # <------------------------------------- 
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
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
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
    vos_shockCapturingFactor = 2.
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
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

ns_nl_atol_res = max(opts.res, 0.001 * he ** 2)
ns_sed_nl_atol_res = max(opts.res, 0.001 * he ** 2)
vof_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
vos_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
ls_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
rd_nl_atol_res =  max(opts.res, 0.005 * he)
mcorr_nl_atol_res =  max(opts.res, 0.001 * he ** 2)
kappa_nl_atol_res =  max(opts.res,  0.005* he ** 2)
dissipation_nl_atol_res =  max(opts.res, 0.005 * he ** 2)
phi_nl_atol_res = max(opts.res, 0.001 * he ** 2)
pressure_nl_atol_res = max(opts.res, 0.001 * he ** 2)


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Turbulence
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Functions for model variables - Initial conditions
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

waterLine_z = opts.waterLevel

def signedDistance(x):
    phi_z = x[1] - waterLine_z
    return phi_z

def vos_signedDistance(x):
    phi_z = x[1] - 0.75*waterLine_z
    return phi_z

class Suspension_class:
    def __init__(self):
        pass
    def uOfXT(self, x, t=0):
        phi = signedDistance(x) + 0.29   # 0.6 is the distance between the free surface and the sediments
        smoothing = (epsFact_consrv_heaviside)*he/2.
        Heav = smoothedHeaviside(smoothing, phi)      
        if phi <= -smoothing:
            return opts.cSed
        elif -smoothing < phi < smoothing:
            return opts.cSed * (1.-Heav)            
        else:
            return 1e-10    

def vos_function(x, t=0):
    phi = signedDistance(x) + 0.29
    smoothing = (epsFact_consrv_heaviside)*he/2.
    Heav = smoothedHeaviside(smoothing, phi)      
    if phi <= -smoothing:
        return opts.cSed
    elif -smoothing < phi < smoothing:
        return opts.cSed * (1.-Heav)            
    else:
        return 1e-10    


Suspension = Suspension_class()

