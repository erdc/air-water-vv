"""
Wave Breaking (Ting & Kirby 1994)
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
    ("water_level", 0.40, "Height of (mean) free surface above bottom"),
    ("free_slip", True, "Free slip BC's enforced (otherwise, no slip)"),
    # tank
    ("tank_dim", (40.0, 1.0), "Dimensions (x,y) of the tank"),
    ("cot_slope", 35.0, "Cotangent slope"),
    ("slope_start_x", 17.0, "Start of constant sloped region."),
    ("slope_initial", (2.0, 0.02), "Dimensions (x,y) of the initial slope - "
                                   "caused by the physical joint between the"
                                   " tank and the sloping bottom."),
    # waves
    ("wave_period", 5., "Period of the waves"),
    ("wave_height", 0.065, "Height of the waves"),
    ("wave_depth", 1., "Wave depth"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 10.655, "Wavelength"),  # calculated by FFT
    ("y_coeff", (0.02107762,
                 0.01588268,
                 0.01059044,
                 0.00658221,
                 0.00394857,
                 0.00233566,
                 0.00138150,
                 0.00082688,
                 0.00050907,
                 0.00033254,
                 0.00024420,
                 0.00021730), "YCoeff array from Fenton calculation tool"),
    ("b_coeff", (0.04205140,
                 0.01597969,
                 0.00713356,
                 0.00328541,
                 0.00151320,
                 0.00068665,
                 0.00030328,
                 0.00012864,
                 0.00005141,
                 0.00001879,
                 0.00000623,
                 0.00000148), "Bcoeff array from calculation tool"),
    # gauges
    ("point_gauge_output", True, "Produce point gauge data"),
    ("column_gauge_output", True, "Produce column gauge data"),
    ("evenly_spaced_gauges", True, "Produce point/column gauge data at evenly"
                                   " spaced (by gauge_dx) intervals."),
    ("experimental_gauges", True, "Produce point/column gauge data at specified"
                                  " points (by gauge_depths)"),
    ("gauge_dx", 0.25, "x spacing between gauges/columns of gauges"),
    ("gauge_depths", [0.169, 0.156, 0.142, 0.128, 0.113, 0.096, 0.079],
     "Gauges placed at specific depths (based on experiment setup)."),
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
slope_x0 = opts.slope_start_x - opts.slope_initial[0]
slope_x1 = opts.slope_start_x
slope_y1 = opts.slope_initial[1]
slope_y2 = (tank_dim[0] - slope_x1) / opts.cot_slope

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
amplitude = waveheight / 2.0

waves = wt.MonochromaticWaves(period=period,
                              waveHeight=waveheight,
                              mwl=waterLine_z,
                              depth=depth,
                              g=g,
                              waveDir=waveDir,
                              wavelength=wavelength,
                              waveType='Fenton',
                              Ycoeff=opts.y_coeff,
                              Bcoeff=opts.b_coeff)

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

sloped_shore = [[[slope_x0, 0],
                 [slope_x1, slope_y1],
                 [tank_dim[0], slope_y2]],]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=sloped_shore)

# # Test Tank
# tank = st.Tank2D(domain=domain,
#                  dim=tank_dim)

# ----- GAUGES ----- #

# Evenly Spaced Gauges

if opts.evenly_spaced_gauges:
    point_gauge_locations = []
    column_gauge_locations = []

    if opts.point_gauge_output or opts.column_gauge_output:
        number_of_gauges = tank_dim[0] / opts.gauge_dx + 1
        for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):
            if gauge_x < slope_x0:
                gauge_y = 0.0
            elif gauge_x > slope_x1:
                gauge_y = (1 / opts.cot_slope) * (gauge_x - slope_x1) \
                          + slope_y1 + 0.001
            else:
                gauge_y = slope_y1 * (gauge_x - slope_x0) \
                          / (slope_x1 - slope_x0) + 0.001

            if waterLine_z - gauge_y > 0.2:
                point_gauge_locations.append((gauge_x, gauge_y, 0), )
            column_gauge_locations.append(((gauge_x, gauge_y, 0.),
                                           (gauge_x, tank_dim[1], 0.)))

        if opts.point_gauge_output:
            tank.attachPointGauges('twp',
                                   gauges=((('u','v'),
                                            point_gauge_locations),
                                           (('p',),
                                            point_gauge_locations),),
                                   fileName='combined_0_0.5_sample_all.txt')

        if opts.column_gauge_output:
            tank.attachLineIntegralGauges('vof',
                                          gauges=(
                                          (('vof',), column_gauge_locations),),
                                          fileName='column_gauges.csv')

if opts.experimental_gauges:
    point_gauge_locations_exp = []
    column_gauge_locations_exp = []

    if opts.point_gauge_output or opts.column_gauge_output:
        gauges_by_depth = opts.gauge_depths
        for depth in gauges_by_depth:
            gauge_y = max(waterLine_z - depth, 0.0)
            gauge_x = slope_x1 + (opts.cot_slope * (gauge_y - slope_y1))

            if waterLine_z - gauge_y > 0.095:
                point_gauge_locations_exp.append((gauge_x, gauge_y, 0), )
            column_gauge_locations_exp.append(((gauge_x, gauge_y, 0.),
                                           (gauge_x, tank_dim[1], 0.)))

        if opts.point_gauge_output:
            tank.attachPointGauges('twp',
                                   gauges=((('u', 'v'),
                                            point_gauge_locations_exp),
                                           (('p',),
                                            point_gauge_locations),),
                                   fileName='exp_gauges.txt')

        if opts.column_gauge_output:
            tank.attachLineIntegralGauges('vof',
                                          gauges=((('vof',),
                                                   column_gauge_locations_exp),
                                                  ),
                                          fileName='exp_column_gauge.csv')

# ----- EXTRA BOUNDARY CONDITIONS ----- #

# open top
tank.BC['y+'].setAtmosphere()

# free/no slip
if opts.free_slip:
    tank.BC['y-'].setFreeSlip()
    tank.BC['x+'].setFreeSlip()
else:  # no slip
    tank.BC['y-'].setNoSlip()
    tank.BC['x+'].setNoSlip()

# waves
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waves,
                                               wind_speed=windVelocity)


# ----- MESH CONSTRUCTION ----- #

he = wavelength / refinement_level
domain.MeshOptions.he = he
st.assembleDomain(domain)

##########################################
# Numerical Options and Other Parameters #
##########################################

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
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
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

ns_nl_atol_res = max(1.0e-10, 0.00001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.00001 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.0001 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.005 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.0001 * he ** 2)
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

#
# def z(x):
#     return x[1] - inflowHeightMean

#
# sigma = omega - k * inflowVelocityMean[0]
# h = waterLine_z  # - transect[0][1] if lower left hand corner is not at z=0


def waveHeight(x, t):
    Y = opts.y_coeff

    waterDepth = waterLine_z
    for i in range(len(Y)):
        waterDepth += Y[i] * cos((i + 1) * theta(x, t)) / k
    return waterDepth


# B = [0.04205140,
#      0.01597969,
#      0.00713356,
#      0.00328541,
#      0.00151320,
#      0.00068665,
#      0.00030328,
#      0.00012864,
#      0.00005141,
#      0.00001879,
#      0.00000623,
#      0.00000148]


# def waveVelocity_u(x, t):
#     wu = 0
#     for i in range(0, int(len(B))):
#         wu += sqrt(abs(g[1]) / k) * (i + 1) * B[i] * cosh(
#             (i + 1) * k * (z(x) + h)) / cosh((i + 1) * k * h) * cos(
#             (i + 1) * theta(x, t))
#
#     return wu
#
#
# def waveVelocity_v(x, t):
#     wv = 0
#     for i in range(0, int(len(B))):
#         wv += sqrt(abs(g[1]) / k) * (i + 1) * B[i] * sinh(
#             (i + 1) * k * (z(x) + h)) / cosh((i + 1) * k * h) * sin(
#             (i + 1) * theta(x, t))
#
#     return wv


#solution variables

def wavePhi(x, t):
    return x[1] - waveHeight(x, t)

# def waveVF(x, t):
#     return smoothedHeaviside(epsFact_consrv_heaviside * he, wavePhi(x, t))
#
#
# def twpflowVelocity_u(x, t):
#     waterspeed = waveVelocity_u(x, t)
#     H = smoothedHeaviside(epsFact_consrv_heaviside * he,
#                           wavePhi(x, t) - epsFact_consrv_heaviside * he)
#     u = H * windVelocity[0] + (1.0 - H) * waterspeed
#     return u
#
#
# def twpflowVelocity_v(x, t):
#     waterspeed = 0.0
#     H = smoothedHeaviside(epsFact_consrv_heaviside * he,
#                           wavePhi(x, t) - epsFact_consrv_heaviside * he)
#     return H * windVelocity[1] + (1.0 - H) * waterspeed
#
#
# def twpflowFlux(x, t):
#     return -twpflowVelocity_u(x, t)
#
#
# def outflowVF(x, t):
#     return smoothedHeaviside(epsFact_consrv_heaviside * he,
#                              x[1] - outflowHeight)
