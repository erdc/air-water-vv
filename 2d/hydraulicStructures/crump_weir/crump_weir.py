"""
Crump Weir
"""
from __future__ import division
from builtins import str
from past.utils import old_div
import numpy as np
from math import sqrt
from proteus import (Domain, Context,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent

opts = Context.Options([
    # water
    ("inflow_level", 1.5, "Upstream water level in m"),
    ("outflow_level", 0.5, "Downstream water level in m"),
    ("inflow_velocity", 1.345, "Inflow velocity in m/s"),
    # tank
    ("tank_dim", (8.0, 3.75), "Dimensions (x,y) of the domain (excluding relaxation zones)"),
    ("sponge_layers", False, "Use sponge layers"),
    ("tank_sponge", (2.,2.), "Length of (generation, absorption) zones, if any"),
    # waves
    ("waves", False, "Use waveTools"),
    # weir
    ("obstacle_dim", (3.0, 0.5), "Width and height of the triangular obstacle in m"),
    ("obstacle_x_start", 3.0, "x coordinate of the start of the obstacle in m"),
    ("cot_upstream_slope", 2.0, "Cotangent of upstream slope of obstacle"),
    # gauges
    ("point_gauge_output", True, "Produce point gauge data"),
    ("column_gauge_output", True, "Produce column gauge data"),
    ("gauge_dx", 0.25, "Horizontal spacing of gauges/gauge columns in m"),
    # refinement
    ("refinement", 40, "Refinement level he = tank_dim[0] / float(4 * refinement - 1)"),
    ("cfl", 0.75, "Target cfl"),
    ("variable_refine_borders", None, "List of vertical borders between "
                                    "refinement regions (include 0 and "
                                    "tank_dim[0] if you add sponge layers "
                                    "and want to differentiate them)"),
    ("variable_refine_levels", None, "List of refinement levels in each region"
                                   " (should have 1 more value than "
                                   "variable_refine_borders as a result)."),
    # run time
    ("T", 30.0, "Simulation time in s"),
    ("dt_fixed", 0.25, "Fixed time step in s"),
    ("dt_init", 0.1, "Minimum initial time step (otherwise dt_fixed/10) in s"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ])

# ----- CONTEXT ------ #

# water
inflow_level = opts.inflow_level
outflow_level = opts.outflow_level

# flow
inflow_velocity = opts.inflow_velocity

# tank
tank_dim = opts.tank_dim
obstacle_dim = opts.obstacle_dim
obstacle_height = obstacle_dim[1]
obstacle_x_start = opts.obstacle_x_start
obstacle_x_end = obstacle_x_start + obstacle_dim[0]
obstacle_x_highest = (obstacle_x_start
                      + obstacle_height * opts.cot_upstream_slope)

##########################################
#     Discretization Input Options       #
##########################################

# ----- From Context.Options ----- #
refinement = opts.refinement
genMesh = opts.gen_mesh

# ----- Structured Meshes ----- #
useHex = False
structured = False

# ----- Parallel Options ----- #
parallelPartitioningType = mt.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

# ---- SpaceOrder & Tool Usage ----- #
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
applyRedistancing = True
timeDiscretization='be'  #'vbdf'#'be','flcbdf'
applyCorrection = True

# ----- INPUT CHECKS ----- #
if spaceOrder not in [1,2]:
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
backgroundDiffusionFactor = 0.01

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
g = [0., -9.81, 0.]

# wind
windVelocity = (0.0, 0.0)

# ----- TIME STEPPING & VELOCITY----- #

T = opts.T
dt_fixed = opts.dt_fixed
dt_init = min(0.1 * dt_fixed, opts.dt_init)
runCFL = opts.cfl
nDTout = int(round(old_div(T, dt_fixed)))

##########################################
#              Mesh & Domain             #
##########################################

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
he = old_div(tank_dim[0], float(4 * refinement - 1))

# ----- TANK ----- #

weir = [[[obstacle_x_start, 0],
         [obstacle_x_highest, obstacle_height],
         [obstacle_x_end, 0]]]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=weir)

# ----- SPONGE & WAVES ----- #
if opts.sponge_layers:
    omega = 1.
    dragAlpha = 5.*omega/nu_0
    tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])
    tank.setAbsorptionZones(x_n=True, dragAlpha=dragAlpha)
    tank.setAbsorptionZones(x_p=True, dragAlpha=dragAlpha)

if opts.waves:
    # TODO: for now this is an actual wave, which we don't want.  We want a placid
    # wave such that our generation and absorption zones enforce a steady flow.
    omega = 2*np.pi/2.
    dragAlpha = 5.*omega/nu_0
    wave = wt.MonochromaticWaves(
        period = 2,
        waveHeight = 0.0,
        mwl = inflow_level,
        depth = inflow_level,
        g = np.array(g),
        waveDir = (1.,0.,0.),
        wavelength = 0.5,
        meanVelocity = np.array([inflow_velocity, 0., 0.])
    )
    tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])
    tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing=3.*he)
    tank.setAbsorptionZones(x_p=True, dragAlpha=dragAlpha)
    

# ----- VARIABLE REFINEMENT ----- #

if opts.variable_refine_borders or opts.variable_refine_levels:
    refinement_borders = opts.variable_refine_borders
    refinement_levels  = opts.variable_refine_levels

    if refinement_borders == None or refinement_levels == None:
        raise ValueError("For variable refinement, variable_refine_borders "
                         "and variable_refine_levels must both be defined.")

    if len(refinement_borders) + 1 != len(refinement_levels):
        raise ValueError("The tank is split into {0} regions, but {1} levels "
                         "of refinement have been "
                         "specified.".format(len(refinement_borders) + 1,
                                             len(refinement_levels)))

    refinement_borders = ([tank.x0 - opts.tank_sponge[0]]
                          + refinement_borders
                          + [tank.x1 + opts.tank_sponge[1]])
    #TODO: Horizontal Variable Refinement
    # Refinement borders should now contain one more element than refinement_levels
    # The borders can be zipped together with the levels:
    #           refinement_level[0] is bordered on the left by refinement_borders[0]
    #                               and on the right by refinement_borders[1]
    #                               and so on (each is bordered by i and i+1)
    # The y borders are just the tank dimensions.
    # This should hold all data necessary in an easy package for the final
    # GMSH box refinement interface.
    raise NotImplementedError("So you can find this unfinished point easier.")

# ----- GAUGES ----- #

column_gauge_locations = []
point_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:

    number_of_gauges = old_div(tank_dim[0], opts.gauge_dx) + 1

    for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):

        if obstacle_x_start <= gauge_x < obstacle_x_highest:
            gauge_y = (obstacle_height
                       / (obstacle_x_highest - obstacle_x_start)
                       * (gauge_x - obstacle_x_start))
        elif obstacle_x_highest <= gauge_x < obstacle_x_end:
            gauge_y = (obstacle_height
                       + obstacle_height
                       / (obstacle_x_end - obstacle_x_highest)
                       * (gauge_x - obstacle_x_highest))
        else:
            gauge_y = 0.
        point_gauge_locations.append((gauge_x, obstacle_height, 0),)
        column_gauge_locations.append(((gauge_x, gauge_y, 0.),
                                       (gauge_x, tank_dim[1], 0.)))

if opts.point_gauge_output:
    tank.attachPointGauges('twp',
                           gauges=((('u','v'), point_gauge_locations),
                                   (('p',), point_gauge_locations)),
                           fileName='combined_gauge_0_0.5_sample_all.txt')

if opts.column_gauge_output:
    tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),
                                  fileName='column_gauge.csv')
    column_over_crest = []
    column_over_crest.append(((opts.obstacle_x_start+old_div(opts.obstacle_dim[0],3.0),opts.obstacle_dim[1],0.0),
                              (opts.obstacle_x_start+old_div(opts.obstacle_dim[0],3.0),opts.tank_dim[1],0.0)))
    tank.attachLineGauges('twp',
                          gauges=((('u'), column_over_crest),),
                          fileName='u_over_crest.csv')                          

# ----- EXTRA BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=outflow_level,
                                                    rhoUp=rho_1,
                                                    rhoDown=rho_0,
                                                    g=g,
                                                    refLevel=tank_dim[1],
                                                    smoothing=3.0*he,
                                                    U=[0.,0.,0.],
                                                    )

tank.BC['x-'].setTwoPhaseVelocityInlet(U=[inflow_velocity,0.,0.],
                                       waterLevel=inflow_level,
                                       smoothing=3.0*he,
                                       )
pAdv = - inflow_velocity
tank.BC['x-'].p_advective.uOfXT = lambda x, t: pAdv

# ----- MESH CONSTRUCTION ----- #

he = he
domain.MeshOptions.he = he
st.assembleDomain(domain)

##########################################
# Numerical Options and Other Parameters #
##########################################

# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# ----- NUMERICAL PARAMETERS ----- #

if useMetrics:
    ns_shockCapturingFactor = 0.75
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.75
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.5
    vof_shockCapturingFactor = 0.75
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor = 0.75
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
#            Signed Distance             #
##########################################

def signedDistance(x):
    if x[0] < obstacle_x_start:
        phi_z = x[1] - inflow_level
    elif x[0] < obstacle_x_end:
        phi_z = x[1] - inflow_level \
        + ((inflow_level - outflow_level) / (obstacle_x_end - obstacle_x_start)
           * (x[0] - obstacle_x_start)
           )
    else:
        phi_z = x[1] - outflow_level

    return phi_z
