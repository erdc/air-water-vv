import sys
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("test_case", None, "case number"),
    ("test_type", 3, "type of test: 1=heave, 2=pitch, or 3=roll"),
    ("water_level", 0.6, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (5., 1.5, 1.5), "Dimensions of the tank"),
    ("tank_sponge", (2., 0.), "Length of absorption zones (front/back, left/right)"),
    #caisson
    ("caisson_dim", (0.594, 0.416, 0.640), "Dimensions of the bar"),
    ("free_x", (0.0, 0.0, 0.0), "Translational DOFs"),
    ("free_r", (1.0, 0.0, 0.0), "Rotational DOFs"),
    ("VCG", 0.175, "vertical position of the barycenter of the caisson"),
    ("draft", 0.425, "Draft of the caisson"),
    ("It", (4.956, 6.620, 5.515), "Inertia tensor: Ixx, Iyy, and Izz components"),
    ("rotation_angle", np.pi/12., "Initial rotation angle (in radians)"),
    ("rotation_axis", (1.,0.,0.), "Axis for initial rotation"),
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("refinement_level", 1 ,"Set maximum element diameter to he/2**refinement_level"),
    ("T", 10.0 ,"Simulation time"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.33 ,"Target cfl"),
    ("nsave",  20,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel")])


# ----- CONTEXT ----- #

# general options
waterLevel = opts.water_level

# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# caisson options
dim = opts.caisson_dim
VCG = opts.VCG
draft = opts.draft
Ixx, Iyy, Izz = opts.It
free_x = opts.free_x
free_r = opts.free_r
rotation_angle = opts.rotation_angle
rotation_axis = opts.rotation_axis


# PREDEFINED CASES

case = opts.test_case
test = opts.test_type
if case is not None:
    # dimension of caisson
    dim1 = (0.594, 0.416, 0.640)
    dim2 = (0.594, 0.608, 0.640)
    dim3 = (1.143, 0.900, 0.640)
    dims = [dim1,  dim1,  dim2,  dim2,  dim2,  dim2,  dim3, dim3,  dim3,  dim3]
    # vertical center of gravity
    VCGs = [0.191, 0.201, 0.188, 0.189, 0.189, 0.203, None,  0.178, 0.183, 0.200]
    # Draft
    drafts = [0.425, 0.500, 0.350, 0.351, 0.425, 0.500, 0.1942, 0.350, 0.425, 0.500]
    # components of mass matrix
    Ixx =  [4.956, 5.356, 9.089, 9.789, 9.789, 11.089, None, 28.944, 32.544, 37.044]
    Iyy =  [6.620, 7.220, 8.848, 9.448, 9.448, 10.648, None, 47.389, 54.689, 62.889]
    Izz =  [5.515, 6.215, 9.307, 10.807, 10.807, 12.307, None, 58.050, 68.550, 79.050]
    # rotation options
    free_x = opts.free_x
    free_r = opts.free_r
    # get right value from index of test case
    if 1 <= case <= 10:
        ind = case-1
        dim = dims[ind]
        VCG = VCGs[ind]
        draft = drafts[ind]
        Ixx = Ixx[ind]
        Iyy = Iyy[ind]
        Izz = Izz[ind]
        if test == 1:
            free_x = (0., 0., 1.)
            free_r = (0., 0., 0.)
        elif test == 2:
            free_x = (0., 0., 0.)
            free_r = (0., 1., 0.)
            rotation_axis = (0., 1., 0.)
        elif test == 3:
            free_x = (0., 0., 0.)
            free_r = (1., 0., 0.)
            rotation_axis = (1., 0., 0.)
    else:
        sys.exit()

coords = (tank_dim[0]/2., tank_dim[1]/2., waterLevel+dim[2]/2.-draft)
barycenter = (coords[0], coords[1], coords[2]-dim[2]/2.+VCG)




# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()


# ----- SHAPES ----- #

front = back = tank_sponge[0]
right = left = tank_sponge[1]
fbt = rlt = False
if front and back:
    fbt = True
if right and left:
    rlt = True

tank = st.Tank3D(domain, tank_dim)
tank.setSponge(front=front, back=back, right=right, left=left)
tank.setAbsorptionZones(front=fbt, back=fbt, right=rlt, left=rlt)

caisson3D = st.Cuboid(domain, dim=dim, coords=coords)
caisson3D.setRigidBody()
caisson3D.setBarycenter(barycenter)
caisson3D.setConstraints(free_x=free_x, free_r=free_r)
# mass is not real mass ---> inerta tensor provided and scaled by mass
mass = 15
caisson3D.setMass(mass)
caisson3D.It = np.array([[Ixx, 0., 0.],
                         [0., Iyy, 0.],
                         [0., 0., Izz]])/caisson3D.mass
caisson3D.rotate(rotation_angle, rotation_axis)
caisson3D.setRecordValues(pos=True, rot=True, F=True, M=True)

# ----- BOUNDARY CONDITIONS ----- #

for bc in caisson3D.BC_list:
    bc.setNoSlip()

tank.BC.bottom.setFreeSlip()
tank.BC.front.setFreeSlip()
tank.BC.left.setFreeSlip()
tank.BC.back.setFreeSlip()
tank.BC.right.setFreeSlip()
tank.BC.top.setOpenAir()
tank.BC.sponge.setNonMaterial()






##########################################
# Numerical Options and other parameters #
##########################################

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., 0., -9.81]
# ------------------



refinement_level = opts.refinement_level
he = (caisson3D.dim[2])/2.0*(0.5**refinement_level)
domain.MeshOptions.elementSize(he)
st.assembleDomain(domain)


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral



quad_order = 3

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = False
openSides = False
openEnd = True
smoothBottom = False
smoothObstacle = False
movingDomain=True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
dt_out =  (T-dt_init)/nDTout
runCFL = opts.cfl

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
useRANS = 1 # 0 -- None
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
nd = 3
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
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
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
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*he**2)
rd_nl_atol_res = max(1.0e-12,0.01*he)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x,t):
    p_L = 0.0
    phi_L = tank.dim[2] - waterLevel
    phi = x[2] - waterLevel
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))




