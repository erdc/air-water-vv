from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.mbd import ChRigidBody as crb
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 1., "Height of free surface above bottom"),
    # tank
    ("tank_dim", (3., 2.,), "Dimensions of the tank"),
    # caisson
    ("caisson_dim", (0.5, 0.2), "Dimensions of the caisson"),
    ("caisson_coords", (1.5, 1.), "Dimensions of the caisson"),
    ("caisson_width", 1., "Width of the caisson"),
    ("free_x", (1., 1., 1.), "Translational DOFs"),
    ("free_r", (0., 0., 0.), "Rotational DOFs"),
    ("VCG", None, "vertical position of the barycenter of the caisson"),
    ("caisson_mass", 50., "Mass of the caisson"),
    ("caisson_inertia", 0.1, "Inertia of the caisson"),
    ("rotation_angle", 0., "intial rotation angle of body (degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("he", 0.01, "mesh element size"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("parallel", True ,"Run in parallel")])


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]

# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

tank_dim = opts.tank_dim



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# caisson options
free_x = opts.free_x
free_r = opts.free_r

caisson_dim = opts.caisson_dim
caisson_coords = opts.caisson_coords
VCG = opts.VCG
if VCG is None:
    VCG = caisson_dim[1]/2.
barycenter = (0, -caisson_dim[1]/2.+VCG, 0.)
inertia = opts.caisson_inertia


caisson = st.Rectangle(domain, dim=opts.caisson_dim)
caisson.setHoles([[0., 0.]])
caisson.holes_ind = np.array([0])  # for gmsh (index of facet that is a hole)

caisson.translate([caisson_coords[0], caisson_coords[1]])

ang = np.radians(opts.rotation_angle)
rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
caisson.rotate(ang, pivot=caisson.barycenter)

# CHRONO
system = crb.System(np.array([0., -9.81, 0.]))
system.setTimeStep(opts.chrono_dt)
body = crb.RigidBody(shape=caisson,
                    system=system,
                    center=caisson.barycenter[:2],
                    rot=rotation_init,
                    mass = opts.caisson_mass,
                    inertia = np.array([0., 0., inertia]),
                    free_x = np.array(opts.free_x),
                    free_r = np.array(opts.free_r))
body.setRecordValues(all_values=True)



for bc in caisson.BC_list:
    bc.setNoSlip()


tank = st.Rectangle(domain, tank_dim)
# let gmsh know that the caisson is IN the tank
tank.setChildShape(caisson, 0)
tank.translate(np.array(tank_dim)/2.)


# ----- BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setNoSlip()

# moving mesh BC
tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes


domain.MeshOptions.he = opts.he
st.assembleDomain(domain)






##########################################
# Numerical Options and other parameters #
##########################################






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
weak_bc_penalty_constant = 10./nu_0#Re
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
epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01

ns_shockCapturingFactor  = 0.5
ns_lag_shockCapturing = True
ns_lag_subgridError = True
ls_shockCapturingFactor  = 0.5
ls_lag_shockCapturing = True
ls_sc_uref  = 1.0
ls_sc_beta  =1.5
vof_shockCapturingFactor = 0.5
vof_lag_shockCapturing = True
vof_sc_uref = 1.0
vof_sc_beta =1.5
rd_shockCapturingFactor  =0.5
rd_lag_shockCapturing = False
epsFact_density    = 3.
epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
epsFact_redistance = 0.33
epsFact_consrv_diffusion = epsFact_consrv_diffusion
redist_Newton = True#False
kappa_shockCapturingFactor = 0.5
kappa_lag_shockCapturing = False#True
kappa_sc_uref = 1.0
kappa_sc_beta =1.5
dissipation_shockCapturingFactor = 0.5
dissipation_lag_shockCapturing = False#True
dissipation_sc_uref = 1.0
dissipation_sc_beta =1.5

he = opts.he
ns_nl_atol_res = max(1.e-8,0.001*he**2)
vof_nl_atol_res = max(1.e-8,0.001*he**2)
ls_nl_atol_res = max(1.e-8,0.001*he**2)
mcorr_nl_atol_res = max(1.e-8,0.1*0.001*he**2)
rd_nl_atol_res = max(1.e-8,0.001*he)
kappa_nl_atol_res = max(1.e-8,0.001*he**2)
dissipation_nl_atol_res = max(1.e-8,0.001*he**2)
mesh_nl_atol_res = max(1.e-8,0.001*he**2)

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
