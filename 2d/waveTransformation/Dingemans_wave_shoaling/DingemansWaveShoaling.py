"""
Dingemans Wave Shoaling
"""
import numpy as np
from math import sqrt, cos, pi, cosh, sinh
from proteus import (Domain, Context,
                     FemTools as ft,
                     #SpatialTools as st,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent

opts = Context.Options([
    # test options
    ("water_level", 0.86, "Height of (mean) free surface above bottom"),
    # tank
    ("tank_dim", (58., 1.26), "Dimensions (x,y) of the tank"),
    ("tank_sponge", (5., 5.), "Length of generation/absorption zone"),
    # waves
    ("wave_period", 2.02, "Period of the waves"),
    ("wave_height", 0.02, "Height of the waves"),
    ("wave_depth", 1., "Wave depth"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 5.037, "Wavelength"),  # calculated by FFT
    ("y_coeff", np.array([0.01246994, 0.00018698, 0.00000300, 0.00000006, 0.00000000,
                          0.00000000, 0.00000000, 0.00000000]), "YCoeff array from Fenton calculation tool"),
    ("b_coeff", np.array([0.01402408, 0.00008097, 0.00000013, 0.00000000, 0.00000000,
                          0.00000000, 0.00000000, 0.00000000]), "Bcoeff array from calculation tool"),
    # gauges
    ("column_gauge_output", True, "Produce column gauge output"),
    ("gauge_dx", 0.25, "Horizontal spacing of point gauges/column gauges"),
    # refinement
    ("refLevel", 300, "Refinement level (w/respect to wavelength)"),
    ("cfl", 0.9, "Target cfl"),
    # run time
    ("T", 30.0, "Simulation time (in numbers of wave_period's)"),
    ("dt_init", 0.1, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")])

# ----- CONTEXT ------ #

# tank
tank_dim = opts.tank_dim

# water
waterLine_z = opts.water_level
waterLine_x = 2 * tank_dim[0]

##########################################
#     Discretization Input Options       #
##########################################

# ----- From Context.Options ----- #
refinement_level = opts.refLevel
genMesh = opts.gen_mesh

# ----- Structured Meshes ----- #
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

# ----- Parallel Options ----- #
parallelPartitioningType = mt.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

# ----- BC & Other Flags ----- #
timeDiscretization='be'  #'vbdf'#'be','flcbdf'
movingDomain = False
checkMass = False
applyCorrection = True
applyRedistancing = True

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
        elementQuadrature = ft.CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd - 1, 2)
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

weak_bc_penalty_constant = 100.0
nLevels = 1
backgroundDiffusionFactor = 0.0

# ----- PHYSICAL PROPERTIES ----- #

# Water
rho_0 = 998.2
nu_0 = 1.004e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

# Surface Tension
sigma_01 = 0.0

# gravity
g = np.array([0., -9.81, 0.])

# wind
windVelocity = (0.0, 0.0, 0.0)

# ----- WAVES ----- #
period = opts.wave_period
waveheight = opts.wave_height
depth = opts.wave_depth
waveDir = np.array(opts.wave_dir)
wavelength = opts.wavelength

waves = wt.MonochromaticWaves(period=period,
                              waveHeight=waveheight,
                              mwl=waterLine_z,
                              depth=depth,
                              g=g,
                              waveDir=waveDir,
                              wavelength=wavelength,
                              waveType='Fenton',
                              Ycoeff=opts.y_coeff,
                              Bcoeff=opts.b_coeff,
                              Nf=len(opts.y_coeff),
                              fast=True)

# ----- TIME STEPPING & VELOCITY----- #

T = opts.T * period
dt_fixed = T
dt_init = min(0.1 * dt_fixed, opts.dt_init)
runCFL = opts.cfl
nDTout = int(round(T / dt_fixed))

##########################################
#              Mesh & Domain             #
##########################################

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #

sloped_shore = [[[9.22, 0.],
                 [9.64, 0.06],
                 [15.01, 0.06],
                 [27.04, 0.66],
                 [31.04, 0.66],
                 [37.07, 0.06],
                 [45.39, 0.06],
                 [45.81, 0.]],]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=sloped_shore)

tank_sponge = opts.tank_sponge
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
tank.setGenerationZones(x_n=True, waves=waves)
tank.setAbsorptionZones(x_p=True)

# ----- GAUGES ----- #

gauge_x = [6.26, 10.26, 12.66, 23.26, 27.26, 29.26, 31.26, 33.66, 36.86, 40.26, 44.26]
gauge_y = []
column_gauge_locations = []

for i in range(len(gauge_x)):
    
    if 9.22 < gauge_x[i] < 9.64:
        gauge_y.append( (gauge_x[i]-9.22)*0.06/(9.64-9.22) )
    elif 9.64 <= gauge_x[i] <= 15.01:
        gauge_y.append(0.06)
    elif 15.01 < gauge_x[i] < 27.04:
        gauge_y.append( 0.06+(gauge_x[i]-15.01)*(0.66-0.06)/(27.04-15.01) )
    elif 27.04 <= gauge_x[i] <= 31.04:
        gauge_y.append(0.66)
    elif 31.04 < gauge_x[i] < 37.07:
        gauge_y.append( 0.66+(gauge_x[i]-31.04)*(0.06-0.66)/(37.07-31.04) )
    elif 37.07 <= gauge_x[i] <= 45.39:
        gauge_y.append(0.06)
    elif 45.39 < gauge_x[i] < 45.81:
        gauge_y.append( 0.06+(gauge_x[i]-45.39)*(0.-0.06)/(45.81-45.39) )
    else:
        gauge_y.append(0.)
        
    column_gauge_locations.append(((gauge_x[i], gauge_y[i], 0.), (gauge_x[i], tank_dim[1], 0.)))

tank.attachLineIntegralGauges('vof', gauges=((('vof',),column_gauge_locations),), fileName='column_gauges.csv')

# ----- EXTRA BOUNDARY CONDITIONS ----- #

he = wavelength / refinement_level
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waves, smoothing=he*3., vert_axis=1)
tank.BC['sponge'].setNonMaterial()

# ----- MESH CONSTRUCTION ----- #

domain.MeshOptions.he = he
st.assembleDomain(domain)

##########################################
# Numerical Options and Other Parameters #
##########################################

# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# ----- NUMERICAL PARAMETERS ----- #

if useMetrics:
    ns_shockCapturingFactor = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.25
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.25
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.25
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_viscosity = epsFact_curvature \
                    = epsFact_vof = ecH = epsFact_consrv_dirac \
                    = 3.0
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
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
    epsFact_consrv_diffusion = 10.0
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

ns_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.005 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)

# ----- TURBULENCE MODELS ----- #
#1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4
else:
    ns_closure = 2
    
##########################################
#       Signed Distance & Wave Phi       #
##########################################

# Wave parameters
omega = 2.0 * pi / period
k = 2.0 * pi / wavelength

def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x, phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x ** 2 + phi_z ** 2)

def theta(x, t):
    return k * x[0] - omega * t

def waveHeight(x, t):
    Y = opts.y_coeff
    waterDepth = waterLine_z
    for i in range(len(Y)):
        waterDepth += Y[i] * cos((i + 1) * theta(x, t)) / k
    return waterDepth

def wavePhi(x, t):
    return x[1] - waveHeight(x, t)
