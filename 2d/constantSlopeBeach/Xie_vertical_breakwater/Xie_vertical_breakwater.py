"""
Vertical Breakwater (Xie 1981)
"""
import numpy as np
from math import sqrt, cos, pi
from proteus import (Domain, Context,
                     FemTools as ft,
                     #SpatialTools as st,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent

opts = Context.Options([
    # test options
    ("water_level", 0.45, "Height of (mean) free surface above bottom"),
    ("free_slip", True, "Free slip BC's enforced (otherwise, no slip)"),
    # tank
    ("tank_dim", (14.5, 0.7), "Dimensions (x,y) of the tank"),
    ("slope_start_ratio", 0.27, "x coordinate of start of slope as a fraction"
                                " of total length (tank_dim[0] * [value] )"),
    ("slope_end_ratio", 0.586, "x coordinate of end of slope as a fraction"
                               " of total length (tank_dim[0] * [value] )"),
    ("slope_height_ratio", 0.214, "y coordinate of end of slope as a fraction"
                                  " of total height (tank_dim[1] * [value] )"),
    # waves
    ("wave_period", 1.53, "Period of the waves"),
    ("wave_height", 0.065, "Height of the waves"),
    ("wave_depth", 1., "Wave depth"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 2.4, "Wavelength"),
    # gauges
    ("point_gauge_output", False, "Produce point gauge data"),
    ("column_gauge_output", False, "Produce column gauge data"),
    ("gauge_dx", 0.25, "x spacing between gauges/columns of gauges"),
    # refinement
    ("refLevel", 100, "Refinement level (w/respect to wavelength)"),
    ("cfl", 0.9, "Target cfl"),
    # run time
    ("T", 10.0, "Simulation time (in numbers of wave_period's)"),
    ("dt_fixed_fraction", 21.0, "Time steps per period: dt_fixed = period/[value]"),
    ("dt_init", 0.001, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")])

# ----- CONTEXT ------ #

# tank
tank_dim = opts.tank_dim
slope_x0 = opts.slope_start_ratio * tank_dim[0]
slope_x1 = opts.slope_end_ratio * tank_dim[0]
slope_y1 = opts.slope_height_ratio * tank_dim[1]

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
timeDiscretization='be'#'vbdf'#'be','flcbdf'
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

waves = wt.RandomWaves(Tp = period,
                       Hs = waveheight,
                       mwl = waterLine_z,
                       depth = depth,
                       waveDir = waveDir,
                       g = g,
                       bandFactor = 2.0,
                       N = 101,
                       spectName = 'JONSWAP')
#[temp] just in case more details are needed on the old intent, this will stay in for a few more commits
# waves = WT.RandomWaves( Tp = period, # Peak period
#                         Hs = waveheight, # Height
#                         d = depth, # Depth
#                         fp = 1./period, #peak Frequency
#                         bandFactor = 2.0, #fmin=fp/Bandfactor, fmax = Bandfactor * fp
#                         N = 101, #No of frequencies for signal reconstruction
#                         mwl = inflowHeightMean, # Sea water level
#                         waveDir = waveDir, # waveDirection
#                         g = g, # Gravity vector, defines the vertical
#                         gamma=1.0, #Pierson Moskowitz spectum for gamma=1.0
#                         spec_fun = WT.JONSWAP)


# ----- TIME STEPPING & VELOCITY----- #

T = opts.T * period
dt_fixed = period / opts.dt_fixed_fraction
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
                 [tank_dim[0], slope_y1]],]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=sloped_shore)

# # Test Tank
# tank = st.Tank2D(domain=domain,
#                  dim=tank_dim)

# ----- GAUGES ----- #

point_gauge_locations = []
column_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:
    number_of_gauges = tank_dim[0] / opts.gauge_dx + 1
    for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):
        if gauge_x < slope_x0:
            gauge_y = 0.0
        elif gauge_x > slope_x1:
            gauge_y = slope_y1
        else:
            gauge_y = slope_y1 * (gauge_x - slope_x0) \
                      / (slope_x1 - slope_x0) + 0.001

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

he = tank_dim[0] / refinement_level
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

def theta(x,t):
    return 2. * pi * ( wavelength *x[0] - (1./period) * t)
#
# def z(x):
#     return x[1] - inflowHeightMean
#
# sigma = omega - k*inflowVelocityMean[0]
# h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0
#
def waveHeightValue(x,t):
    return waveheight + amplitude*cos(theta(x,t))
#
# def waveVelocity_u(x,t):
#     return sigma*amplitude*cosh(k*(z(x)+h))*cos(theta(x,t))/sinh(k*h)
#
# def waveVelocity_v(x,t):
#     return sigma*amplitude*sinh(k*(z(x)+h))*sin(theta(x,t))/sinh(k*h)
#
# #solution variables
#
def wavePhi(x,t):
    return x[1] - waveHeightValue(x,t)
#
# def waveVF(x,t):
#     return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))
#
# def twpflowVelocity_u(x,t):
#     waterspeed = waveVelocity_u(x,t)
#     H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
#     u = H*windVelocity[0] + (1.0-H)*waterspeed
#     return u
#
# def twpflowVelocity_v(x,t):
#     waterspeed = 0.0
#     H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
#     return H*windVelocity[1]+(1.0-H)*waterspeed
#
# def twpflowFlux(x,t):
#     return -twpflowVelocity_u(x,t)
#
# def outflowVF(x,t):
#     return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)
