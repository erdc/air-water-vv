from proteus import Domain, Context, Comm
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.Profiling import logEvent
from proteus.mbd import ChRigidBody as crb
from math import *
import numpy as np

opts=Context.Options([
    # predefined test cases
    ("water_level", 0.9, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (1., 1.8,), "Dimensions of the tank"),
    ("tank_sponge", (2., 2.), "Length of relaxation zones (front/back, left/right)"),
    ("tank_BC", 'FreeSlip', "Tank boundary conditions: NoSlip or FreeSlip"),
    ("gauge_output", False, "Places Gauges in tank"),
    ("gauge_fixed", False, "Places Gauges in tank"),
    # caisson
    ("caisson", True, "Switch on/off caisson"),
    ("caisson_dim", (0.3, 0.1), "Dimensions of the caisson"),
    ("caisson_ycoord", 0.9, "y coord of the caisson center"),
    ("caisson_xcoord", 0.5, "x coord of the caisson center"),
    ("caisson_width", 0.9, "Width of the caisson"),
    ("caisson_BC", 'noslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (0., 0., 0.), "Translational DOFs"),
    ("free_r", (0., 0., 1.0), "Rotational DOFs"),
    ("VCG", 0.05, "vertical position of the barycenter of the caisson"),
    ("caisson_mass", 125., "Mass of the caisson"),
    ("caisson_inertia", 0.236, "Inertia of the caisson"),
    ("rotation_angle", 15., "Initial rotation angle (in degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # mesh refinement
    ("he", 0.04, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("sc", 0.25, "shockCapturing factor"),
    ("weak_factor", 10., "weak bc penalty factor"),
    ("strong_dir", False, "strong dirichlet (True/False)"),
    ("parallel", True ,"Run in parallel")])


# ----- CONTEXT ------ #

waterLevel = opts.water_level
rotation_angle = np.radians(opts.rotation_angle)

tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge


dim = opts.caisson_dim
VCG = opts.VCG
if VCG is None:
    VCG = dim[1]/2.
free_x = opts.free_x
free_r = opts.free_r
rotation = np.radians(opts.rotation_angle)
coords = [0,0]
coords[0] = opts.caisson_xcoord or tank_dim[0]/2.
coords[1] = opts.caisson_ycoord or tank_dim[1]/2.
barycenter = (0, -dim[1]/2.+VCG, 0.)
width = opts.caisson_width
inertia = opts.caisson_inertia/width


caisson_coords = coords
caisson_dim = opts.caisson_dim

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options
if opts.caisson is True:
    caisson = st.Rectangle(domain, dim=opts.caisson_dim)

    ang = rotation_angle
    caisson.setHoles([[0., 0.]])
    caisson.holes_ind = np.array([0])
    caisson.translate([caisson_coords[0], caisson_coords[1]])
    # system = crb.System(np.array([0., -9.81, 0.]))
    # rotation = np.array([1, 0., 0., 0.])
    rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)
    system = crb.ProtChSystem(np.array([0., -9.81, 0.]))
    system.setTimeStep(opts.chrono_dt)
    body = crb.ProtChBody(shape=caisson,
                        system=system)
    from proteus.mbd import pyChronoCore as pych
    x, y, z = caisson.barycenter
    pos = pych.ChVector(x, y, z)
    e0, e1, e2, e3 = rotation_init
    rot = pych.ChQuaternion(e0, e1, e2, e3)
    inertia = pych.ChVector(1., 1., inertia)
    # chrono functions
    body.ChBody.SetPos(pos)
    body.ChBody.SetRot(rot)
    body.ChBody.SetMass(opts.caisson_mass)
    body.ChBody.SetInertiaXX(inertia)
    # customised functions
    body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))
    body.setRecordValues(all_values=True)
    for bc in caisson.BC_list:
        if opts.caisson_BC == 'noslip':
            bc.setNoSlip()
        if opts.caisson_BC == 'freeslip':
            bc.setFreeSlip()


# ----- SHAPES ----- #
tank = st.Tank2D(domain, tank_dim)

# Generation / Absorption zones
dragAlpha = 10*np.pi/1.004e-6
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if left:
    tank.setAbsorptionZones(x_n=left, dragAlpha=dragAlpha)
    tank.BC['x-'].setNoSlip()
else:
    tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right, dragAlpha=dragAlpha)
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)


# ----- BOUNDARY CONDITIONS ----- #

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



domain.MeshOptions.he = opts.he
domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
mesh_fileprefix ='mesh0'
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)

# ASSEMBLE DOMAIN
st.assembleDomain(domain)

# MESH REFINEMENT

he = opts.he




##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = opts.weak_factor/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = opts.useRANS # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
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
if opts.sc == 0.25:
    sc = 0.25 # default: 0.5. Test: 0.25
    sc_beta = 1. # default: 1.5. Test: 1.
    epsFact_consrv_diffusion = 0.1 # default: 1.0. Test: 0.1. Safe: 10.
elif opts.sc == 0.5:
    sc = 0.5
    sc_beta = 1.5
    epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
else:
    import sys
    sys.quit()
ns_forceStrongDirichlet = opts.strong_dir
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
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
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

tolfac = 0.001
mesh_tol = 0.001
ns_nl_atol_res = max(1.0e-8,tolfac*he**2)
vof_nl_atol_res = max(1.0e-8,tolfac*he**2)
ls_nl_atol_res = max(1.0e-8,tolfac*he**2)
mcorr_nl_atol_res = max(1.0e-8,0.1*tolfac*he**2)
rd_nl_atol_res = max(1.0e-8,tolfac*he)
kappa_nl_atol_res = max(1.0e-8,tolfac*he**2)
dissipation_nl_atol_res = max(1.0e-8,tolfac*he**2)
mesh_nl_atol_res = max(1.0e-8,mesh_tol*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
