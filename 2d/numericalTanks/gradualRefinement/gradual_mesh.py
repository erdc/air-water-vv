"""
Linear Wave Theory
"""
import numpy as np
import math
from proteus import (Domain, Context,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

opts = Context.Options([
    # test options
    ("water_level", 1., "Height of free surface above bottom"),
    # tank
    ("tank_dim", (10., 2.,), "Dimensions of the tank"),
    #gravity 
    ("g", [0, -9.81, 0], "Gravity vector"),
    # refinement
    ("cfl", 0.33, "Target cfl"),
    # run time
    ("T", 30.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("use_gmsh", True, "use gmsh"),
    ("refinement_grading", 1.2, "refinement grading"),
    ("refinement_type", 'point', "refinement types: 'band' or 'point'"),
    ("nsave", 20., "Number of time steps to save per period"),
    ("parallel", True, "Run in parallel")])

# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

# tank options
tank_dim = opts.tank_dim

##########################################
#     Discretization Input Options       #
##########################################

# ----- From Context.Options ----- #
genMesh = opts.gen_mesh
useHex = False
structured = False

# ----- SpaceOrder & Tool Usage ----- #
spaceOrder = 1
useOldPETSc = False
useSuperlu = False
useRBLES = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega

# ----- BC & Other Flags ----- #
movingDomain = False
checkMass = False
applyCorrection = True
applyRedistancing = True
freezeLevelSet = True

# ----- INPUT CHECKS ----- #
if spaceOrder not in [1, 2]:
    raise ValueError("INVALID: spaceOrder(" + str(spaceOrder) + ")")

if useRBLES not in [0.0, 1.0]:
    raise ValueError("INVALID: useRBLES(" + str(useRBLES) + ")")

if useMetrics not in [0.0, 1.0]:
    raise ValueError("INVALID: useMetrics(" + str(useMetrics) + ")")

# ----- DISCRETIZATION ----- #

nd = 2
if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = ft.C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd - 1, 3) #[temp] 3? Others have 2.
    else:
        basis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = ft.C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, 4)

##########################################
#   Physical, Time, & Misc. Parameters   #
##########################################

# ----- PHYSICAL PROPERTIES ----- #

# Water
rho_0 = 998.2
nu_0 = 1.004e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

# Surface Tension
sigma_01 = 0.0

# Gravity
g = opts.g

# ----- TIME STEPPING & VELOCITY ----- #

runCFL = opts.cfl
T = opts.T
dt_init = opts.dt_init
nDTout = int(opts.T*opts.nsave)
dt_out= (T-dt_init)/nDTout

# ----- MISC ----- #

weak_bc_penalty_constant = 10 / nu_0
nLevels = 1
backgroundDiffusionFactor = 0.01

##########################################
#              Mesh & Domain             #
##########################################

# ----- DOMAIN ----- #

if useHex:
    nnx = refinement_x + 1
    nny = refinement_y + 1
    hex = True
    domain = Domain.RectangularDomain(tank_dim)
elif structured:
    nnx = refinement_x
    nny = refinement_y
    domain = Domain.RectangularDomain(tank_dim)
    boundaryTags = domain.boundaryTags
else:
    domain = Domain.PlanarStraightLineGraphDomain()

#refinement
he = 0.008
smoothing = he*3.

# ----- TANK ------ #

tank = st.Tank2D(domain, tank_dim)

tank.BC['y+'].setAtmosphere()
tank.BC['x-'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()


# ----- MESH CONSTRUCTION ----- #


domain.MeshOptions.he = he
st.assembleDomain(domain)


grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
geofile = 'mesh'+str(int(he*1000))
if opts.use_gmsh is True:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    he_max = 10.
    field_list = []
    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)

    def dist_plane(xn, xp, plane='x'):
        x_range = abs(xp-xn)
        dist = '0.5*(abs({plane}-({xn}))+abs({plane}-({xp}))-{x_range})'.format(xn=xn, xp=xp, x_range=x_range, plane=plane)
        return dist

    me1 = py2gmsh.Fields.MathEval(mesh=mesh)
    if opts.refinement_type == 'band':
        box = 0.  # width of band around free surface
        dist_x = dist_plane(xn=0, xp=tank_dim[0], plane='x')
        dist_y = dist_plane(xn=waterLevel-box, xp=waterLevel+box, plane='y')
    elif opts.refinement_type == 'point':
        dist_x = '({center_x}-x)'.format(center_x=tank_dim[0]/2.)
        dist_y = '({center_y}-y)'.format(center_y=tank_dim[1]/2.)
    dist = 'sqrt(({dist_x})^2+({dist_y})^2)'.format(dist_x=dist_x, dist_y=dist_y)
    me1.F = mesh_grading(start=dist, he=he, grading=grading)
    field_list += [me1]

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')

domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.gen_mesh
domain.geofile = geofile
domain.MeshOptions.use_gmsh = opts.use_gmsh

# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# ----- NUMERICAL PARAMETERS ----- #

if useMetrics:
    ns_shockCapturingFactor = 0.25
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
    epsFact_density = epsFact_viscosity = epsFact_curvature \
                    = epsFact_vof = ecH = epsFact_consrv_dirac \
                    = 3.0
    epsFact_redistance = 1.5
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
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
    epsFact_density = epsFact_viscosity = epsFact_curvature \
                    = epsFact_vof = ecH = epsFact_consrv_dirac \
                    = 1.5
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

# ----- NUMERICS: TOLERANCES ----- #

ns_nl_atol_res = max(1.0e-10, 0.001 * domain.MeshOptions.he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.001 * domain.MeshOptions.he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.001 * domain.MeshOptions.he ** 2)
mcorr_nl_atol_res = max(1.0e-10, 0.0001 * domain.MeshOptions.he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.01 * domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-10, 0.001 * domain.MeshOptions.he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.001 * domain.MeshOptions.he ** 2)

# ----- TURBULENCE MODELS ----- #
#1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure = 4
else:
    ns_closure = 0

##########################################
#            Boundary Edit               #
##########################################
def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd - 1] - waterLevel
    phi = x[nd - 1] - waterLevel
    return p_L - g[nd - 1] * (rho_0 * (phi_L - phi) + (rho_1 - rho_0) * (
    smoothedHeaviside_integral(ecH * he, phi_L)
    - smoothedHeaviside_integral(ecH * he, phi)))

tank.BC['y+'].p_dirichlet.uOfXT = twpflowPressure_init
