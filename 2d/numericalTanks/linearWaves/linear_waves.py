"""
Linear Wave Theory
"""
from __future__ import division
from builtins import str
from past.utils import old_div
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
    ("water_level", 1., "Water level from y=0"),
    # tank
    ("tank_dim", (15., 1.5,), "Dimensions of the operational domain of the tank in m (l x h)"),
    ("generation", True, "Generate waves at the left boundary (True/False)"),
    ("absorption", True, "Absorb waves at the right boundary (True/False)"),
    ("tank_sponge", (5., 10.), "Length of generation/absorption zone in m (left, right)"),
    ("free_slip", True, "Should tank walls have free slip conditions "
                        "(otherwise, no slip conditions will be applied)."),
    #gravity 
    ("g", [0, -9.81, 0], "Gravity vector in m/s^2"),
    # waves
    ("wave_period", 1.94, "Period of the waves in s"),
    ("wave_height", 0.025, "Height of the waves in m"),
    ("depth", 1., "Wave depth in m"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 5., "Wavelength in m"),
    # probe dx
    ("point_gauge_output", True, "Generate point gauge output"),
    ("column_gauge_output", True, "Generate column gauge output"),
    ("gauge_dx", 0.25, "Horizontal spacing of point gauges/column gauges in m"),
    # refinement
    ("refLevel", 100, "Refinement level (w/respect to wavelength)"),
    ("cfl", 0.33, "Target cfl"),
    # run time
    ("T", 0.1, "Simulation time in s"),
    ("dt_init", 0.001, "Initial time step in s"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("useHex", False, "Use (hexahedral) structured mesh"),
    ("structured", False, "Use (triangular/tetrahedral) structured mesh"),
    ("nperiod", 10., "Number of time steps to save per period"),
    ])

# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

# waves
period = opts.wave_period
height = opts.wave_height
mwl = opts.water_level
depth = opts.depth
direction = opts.wave_dir
wave = wt.MonochromaticWaves(period, height, mwl, depth, np.array(opts.g), direction)

# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

##########################################
#     Discretization Input Options       #
##########################################

# ----- From Context.Options ----- #
refinement_level = opts.refLevel
genMesh = opts.gen_mesh
useHex = opts.useHex
structured = opts.structured

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
dt_out = old_div(opts.wave_period, opts.nperiod)
nDTout = int(round(old_div(T, dt_out)))

# ----- MISC ----- #

weak_bc_penalty_constant = old_div(10, nu_0)
nLevels = 1
backgroundDiffusionFactor = 0.01

##########################################
#              Mesh & Domain             #
##########################################

# ----- DOMAIN ----- #

#[temp] an attempt to match the intentions of refLevel instead of refinement
#[temp] (wavelength based instead of dimension based)
refinement_x = int(math.ceil(refinement_level * (tank_dim[0] + sum(tank_sponge))
                        / opts.wavelength))
refinement_y = int(math.ceil(refinement_level * (tank_dim[1])
                        / opts.wavelength))
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
he = old_div(opts.wavelength, refinement_level)
smoothing = he*3.

# ----- TANK ------ #

tank = st.Tank2D(domain, tank_dim)
omega = 2.*math.pi/period
dragAlpha = 5.*omega/1e-6

# ----- GENERATION / ABSORPTION LAYERS ----- #

tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])

if opts.generation:
    tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing=smoothing)
if opts.absorption:
    tank.setAbsorptionZones(x_p=True, dragAlpha=dragAlpha)

# ----- BOUNDARY CONDITIONS ----- #

# waves
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

# open top
tank.BC['y+'].setAtmosphere()

if opts.free_slip:
    tank.BC['y-'].setFreeSlip()
    tank.BC['x+'].setFreeSlip()
    if not opts.generation:
        tank.BC['x-'].setFreeSlip()
else:  # no slip
    tank.BC['y-'].setNoSlip()
    tank.BC['x+'].setNoSlip()
    if not opts.generation:
        tank.BC['x-'].setNoSlip()

# sponge
tank.BC['sponge'].setNonMaterial()

# ----- GAUGES ----- #

column_gauge_locations = []
point_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:
    gauge_y = waterLevel - 0.5 * depth
    number_of_gauges = old_div(tank_dim[0], opts.gauge_dx) + 1
    for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):
        point_gauge_locations.append((gauge_x, gauge_y, 0), )
        column_gauge_locations.append(((gauge_x, 0., 0.),
                                       (gauge_x, tank_dim[1], 0.)))

if opts.point_gauge_output:
    tank.attachPointGauges('twp',
                           gauges=((('p',), point_gauge_locations),),
                           fileName='pressure_gaugeArray.csv')

if opts.column_gauge_output:
    tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),
                                  fileName='column_gauges.csv')

# ----- MESH CONSTRUCTION ----- #


domain.MeshOptions.he = he
st.assembleDomain(domain)

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
    smoothedHeaviside_integral(ecH * domain.MeshOptions.he, phi_L)
    - smoothedHeaviside_integral(ecH * domain.MeshOptions.he, phi)))

tank.BC['y+'].p_dirichlet.uOfXT = twpflowPressure_init
