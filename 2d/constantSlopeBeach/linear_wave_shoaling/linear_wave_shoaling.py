"""
Linear Wave Shoaling
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
    ("water_level", 1.0, "Height of (mean) free surface above bottom"),
    ("free_slip", True, "Free slip BC's enforced (otherwise, no slip)"),
    # tank
    ("tank_dim", (35.0, 1.50), "Dimensions (x,y) of the tank"),
    ("slope_start", 15.0, "x coordinate of start of slope"),
    ("slope_length", 6.0, "Horizontal length of the slope"),
    ("cot_slope", 10.0, "Cotangent of the slope of the bottom."),
    # absorption / generation zones
    ("generation", True, "Generate waves on the left boundary with a sponge layer"),
    ("absorption", True, "Absorb waves on the right boundary"),
    ("sponge_dim_wavelength", (1.,2.), "Dimensions of the sponge layers (if "
                                       "any) by number of wavelengths."),
    # waves
    ("wave_period", 1.94, "Period of the waves"),
    ("wave_height", 0.025, "Height of the waves"),
    ("wave_depth", 1.0, "Wave depth"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 5., "Wavelength"),
    # gauges
    ("point_gauge_output", False, "Produce point gauge data"),
    ("column_gauge_output", False, "Produce column gauge data"),
    ("gauge_dx", 0.25, "x spacing between gauges/columns of gauges"),
    # refinement
    ("refLevel", 100, "Refinement level (w/respect to wavelength)"),
    ("cfl", 0.9, "Target cfl"),
    # run time
    ("T", 40.0, "Simulation time (in numbers of wave_period's)"),
    ("dt_init", 0.1, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")])

# ----- CONTEXT ------ #

# tank
tank_dim = opts.tank_dim
slope_start_x = opts.slope_start
slope_end_x = slope_start_x + opts.slope_length
slope_height = (1./opts.cot_slope) * opts.slope_length

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
timeDiscretization='be'  #'vbdf'/ 'be'/ 'flcbdf'
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
                              wavelength=wavelength)

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

sloped_shore = [[[slope_start_x, 0],
                 [slope_end_x, slope_height],
                 [tank_dim[0], slope_height]],]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=sloped_shore)

# # Test Tank
# tank = st.Tank2D(domain=domain,
#                  dim=tank_dim)

# ----- SPONGE LAYERS ----- #
tank.setSponge(opts.sponge_dim_wavelength[0] * wavelength,
               opts.sponge_dim_wavelength[1] * wavelength)

if opts.generation:
    tank.setGenerationZones(waves=waves,
                            x_n=opts.generation)
if opts.absorption:
    tank.setAbsorptionZones(x_p=opts.absorption)

# ----- GAUGES ----- #

point_gauge_locations = []
column_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:
    number_of_gauges = tank_dim[0] / opts.gauge_dx + 1
    for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):
        if gauge_x < slope_start_x:
            gauge_y = 0.0
        elif gauge_x > slope_end_x:
            gauge_y = slope_height
        else:
            gauge_y = (1./opts.cot_slope) * (gauge_x - slope_start_x)

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
if not opts.generation:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waves,
                                                   wind_speed=windVelocity)

tank.BC['sponge'].setNonMaterial()

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


# Wave parameters
omega = 2.0 * pi / period
k = 2.0 * pi / wavelength

def theta(x, t):
    return k * x[0] - omega * t + pi / 2.0


# def z(x):
#     return x[1] - inflowHeightMean
#
#
# sigma = omega - k * inflowVelocityMean[0]
# h = inflowHeightMean  # - transect[0][1] if lower left hand corner is not at z=0


def waveHeight(x, t):
    return waterLine_z + amplitude * cos(theta(x, t))

#
# def waveVelocity_u(x, t):
#     return sigma * amplitude * cosh(k * (z(x) + h)) * cos(theta(x, t)) / sinh(
#         k * h)
#
#
# def waveVelocity_v(x, t):
#     return sigma * amplitude * sinh(k * (z(x) + h)) * sin(theta(x, t)) / sinh(
#         k * h)


#solution variables

def wavePhi(x, t):
    return x[1] - waveHeight(x, t)
#
#
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
#     waterspeed = waveVelocity_v(x, t)
#     H = smoothedHeaviside(epsFact_consrv_heaviside * he,
#                           wavePhi(x, t) - epsFact_consrv_heaviside * he)
#     return H * windVelocity[1] + (1.0 - H) * waterspeed
#
#
# def twpflowFlux(x, t):
#     return -twpflowVelocity_u(x, t)
#
#
# outflowHeight = inflowHeightMean
#
#
# def outflowVF(x, t):
#     return smoothedHeaviside(epsFact_consrv_heaviside * he,
#                              x[1] - outflowHeight)
#
#
# def outflowPhi(x, t):
#     return x[1] - outflowHeight
#
#
# def outflowPressure(x, t):
#     if x[1] > inflowHeightMean:
#         return (L[1] - x[1]) * rho_1 * abs(g[1])
#     else:
#         return (L[1] - inflowHeightMean) * rho_1 * abs(g[1]) + (
#                                                                inflowHeightMean -
#                                                                x[
#                                                                    1]) * rho_0 * abs(
#             g[1])
#
#
#         #p_L = L[1]*rho_1*g[1]
#         #phi_L = L[1] - outflowHeight
#         #phi = x[1] - outflowHeight
#         #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
#         #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
#
#
# def twpflowVelocity_w(x, t):
#     return 0.0
#
#
# def zeroVel(x, t):
#     return 0.0
