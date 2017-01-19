from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.75, "Height of free surface above bottom"),
    # Geometry
    ('Lgen', 0.5, 'Genaration zone in terms of wave lengths'),
    ('Labs', 0.5, 'Absorption zone in terms of wave lengths'),
    ('Ls', 0.5, 'Length of domain from genZone to the front toe of rubble mound in terms of wave lengths'),
    ('Lend', 0.5, 'Length of domain from absZone to the back toe of rubble mound in terms of wave lengths'),
    ('th', 1.25, 'Total height of numerical tank'),
    # waves
    ('wave', True, 'Turn on wave generation'),
    ('waveType', 'Linear', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 2.0, "Period of the waves"),
    ("wave_height", 0.1, "Height of the waves"), 
    ('wavelength', 0.0, 'Wavelength only if Fenton is activated'), 
    ('Ycoeff', [ 0.0 ], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [ 0.0 ], 'Bcoeff only if Fenton is activated'),
    # circle2D
    ("circle2D", True, "Switch on/off circle2D"),
    ('radius', 0.075, 'Radius of the circle2D'),
    ('width', 1.0, 'Z-dimension of the circle2D'),
    ('mass', 50.00, 'Mass of the caisson2D [kg]'),
    ("rotation", False, "Initial position for free oscillation"),
    ('circleBC', 'FreeSlip', 'circle2D boundaries: NoSlip or FreeSlip'),
    ('InputMotion', True, 'If True, set a motion as input rather than calculate it'),
    ('scheme', 'Runge_Kutta', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Central_Difference)'),
    ("springs", True, "Switch on/off soil module"),
    ("Kx", 7.0*(10.**2), "Horizontal stiffness in N/m"),
    ("Ky", 3.5*(10.**3), "Vertical stiffness in N/m"),
    ("Krot", 0.0, "Rotational stiffness in N*m"),
    ("Cx", 0.0, "Damping factor in N/m*s "),
    ("Cy", 0.0, "Damping factor in N/m*s "),
    ("Crot", 0.0, "Rotational damping factor in N*m*s "),
    # numerical options
    ("refinement_level", 100.,"he=walength/refinement_level"),
    ("he", 0.00,"Mesh size"),
    ("cfl", 0.90 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', True, "Moving domain and mesh option"),
    ('conservativeFlux', False,'Fix post-processing velocity bug for porous interface'),
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

if opts.wave==True:
    if opts.waveType=='Linear':
        waveinput = wt.MonochromaticWaves(period=period,
                                      waveHeight=waveHeight,
                                      mwl=mwl,
                                      depth=waterLevel,
                                      g=g,
                                      waveDir=waveDir,
                                      wavelength=None, # if wave is linear I can use None
                                      waveType="Linear")

    if opts.waveType=='Fenton':
        waveinput = wt.MonochromaticWaves(period=period,
                                      waveHeight=waveHeight,
                                      mwl=mwl,
                                      depth=waterLevel,
                                      g=g,
                                      waveDir=waveDir,
                                      wavelength=opts.wavelength, # if wave is linear I can use None
                                      waveType="Fenton",
                                      Ycoeff=opts.Ycoeff,
                                      Bcoeff=opts.Bcoeff,
                                      )

#---------Domain Dimension

nd = 2.
if opts.wave==True:
    wl = waveinput.wavelength
else:
    wl=5.0

#---------MESH SIZE
if opts.he == 0.0:
    he = wl/opts.refinement_level
else:
    he = opts.he 



####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


if opts.circle2D:
    L_leftSpo  = opts.Lgen*wl
    L_rightSpo = opts.Labs*wl

#-circle
    radius=opts.radius
    b=2.0*radius
    h=2.0*radius

#-Tank
    x1=0.0+L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+b
    x4=x3+opts.Lend*wl
    x5=x4+L_rightSpo

    xc1=(x2+x3)*0.5
    yc1=0.5

    tank_dim = [x5, opts.th]

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([0., -1.,0.]),
                        'y+': np.array([0., -1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'circle': None,
                           }
boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                    'circle':6,
                       }


##############################################################################################################################################################################################################
# circle2D
############################################################################################################################################################################################################
if opts.circle2D:
    width=opts.width           # The 3rd dimension
    mass=opts.mass             # kg
    I=100000000000.0

# --- Shape properties setup
    circle = st.Circle(domain=domain, radius=opts.radius, coords=(xc1, yc1), barycenter=(xc1, yc1), nPoints=32)

# --- Body properties setup
    circle2D = bd.circularBody(shape=circle, substeps=20)
    free_x=(0.0, 0.0, 0.0) # Translational DOFs
    free_r=(0.0, 0.0, 0.0) # Rotational DOFs
    if opts.movingDomain==True:
        free_x=(1.0, 1.0, 0.0) # Translational DOFs
        free_r=(0.0, 0.0, 0.0) # Rotational DOFs
    circle2D.setMass(mass)
    circle2D.It= I/circle2D.mass/width
    circle2D.setConstraints(free_x=free_x, free_r=free_r)
    circle2D.setNumericalScheme(scheme=opts.scheme)
    circle2D.setRecordValues(filename='circle2D', all_values=True)

#----- Spring setup
    springs=opts.springs
    Kx = opts.Kx
    Ky = opts.Ky
    Krot = opts.Krot
    Cx = opts.Cx
    Cy = opts.Cy
    Crot = opts.Crot
    circle2D.setSprings(springs, Kx, Ky, Krot, Cx, Cy, Crot)


##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

if opts.circle2D==True:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2,  0.0],#2
              [x3,  0.0],#3
              [x4,  0.0],#4
              [x5,  0.0],#5
              [x5,  tank_dim[1]],#6
              [x4,  tank_dim[1]],#7
              [x1,  tank_dim[1]],#8
              [0.0, tank_dim[1]],#9
              ]

    vertexFlags=np.array([1, 1, 1, 1, 1, 1,
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
              [9,0],

              [1,8],
              [4,7],
             ]

    segmentFlags=np.array([ 1, 1, 1, 1,   
                            1,                        
                            2, 3, 3, 3, 4,
                            5, 5,
                         ])


regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2])



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

if opts.wave==True:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1)
else:
    tank.BC['x-'].setFreeSlip()
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

if opts.movingDomain==True:
    for tb in [tank.BC['x+'], tank.BC['x-'], tank.BC['y+'], tank.BC['y-'], tank.BC['sponge']]:
        tb.setFixedNodes()




########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################

tank.setAbsorptionZones(flags=2, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(tank_dim[0]-L_rightSpo/2.), 0., 0.),
                        )


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = 30.*period

gauge_dx=0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)
PG=[]
zProbes=opts.water_level*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

levelset_output=ga.PointGauges(gauges=((('phi',),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='levelset_gauges.csv')


######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

T = 6.*opts.wave_period
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

