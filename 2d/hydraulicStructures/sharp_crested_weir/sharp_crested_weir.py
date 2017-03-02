"""
A Sharp Crested Weir
"""
import numpy as np
from math import sqrt
from proteus import (Domain, Context,
                     FemTools as ft,
                     #SpatialTools as st,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent

opts = Context.Options([
    # test options
    ("waves", False, "Generate waves - uses sponge layers."),
    ("air_vent", False, "Include an air vent in the obstacle."),
    # water
    ("water_level", 1.4, "Height of (mean) free surface above bottom"),
    ("water_width_over_obst", 0.02, "Initial width of free surface relative to"
                                   " the obstacle location"),
    ("outflow_level", -1., "Height of (mean) free surface of water outflow "
                          "give a negative number if no initial outflow."),
    ("inflow_velocity", 0.34, "Wave or steady water inflow velocity"),
    ("outflow_velocity", 0.0, "Initial wave or steady water outflow velocity"),
    # tank
    ("tank_dim", (7.5, 1.8), "Dimensions (x,y) of the tank"),
    ("tank_sponge", (0., 0.), "Length of (generation, absorption) zones, if any"),
    ("obstacle_dim", (0.01, 1.00), "Dimensions (x,y) of the obstacle."),
    ("obstacle_x_start", 5.0, "x coordinate of the start of the obstacle"),
    # gauges
    ("point_gauge_output", False, "Produce gauge data"),
    ("column_gauge_output", False, "Produce column gauge data"),
    ("gauge_dx", 0.5, "Horizontal spacing of gauges/gauge columns."),
    ("point_gauge_y", 0.5, "Height of point gauge placement"),
    # refinement
    ("refinement", 40, "Refinement level"),
    ("cfl", 0.9, "Target cfl"),
    ("variable_refine_borders", None, "List of vertical borders between "
                                    "refinement regions (include 0 and "
                                    "tank_dim[0] if you add sponge layers "
                                    "and want to differentiate them)"),
    ("variable_refine_levels", None, "List of refinement levels in each region"
                                   " (should have 1 more value than "
                                   "variable_refine_borders as a result)."),
    # run time
    ("T", 10.0, "Simulation time"),
    ("dt_fixed", 0.02, "Fixed time step"),
    ("dt_init", 0.001, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")])

# ----- CONTEXT ------ #

# water
waterLine_z = opts.water_level
waterLine_x = opts.water_width_over_obst

if opts.outflow_level < 0.0:
    outflow_level = -(opts.tank_dim[0] ** 2) - (opts.tank_dim[1] ** 2)
else:
    outflow_level = opts.outflow_level

# flow
inflow_velocity = opts.inflow_velocity
outflow_velocity = opts.outflow_velocity

# tank
tank_dim = opts.tank_dim
obstacle_dim = opts.obstacle_dim
obstacle_x_start = opts.obstacle_x_start
obstacle_x_end = obstacle_x_start + obstacle_dim[0]
obstacle_height = obstacle_dim[1]

# air vent
if opts.air_vent:
    air_vent = True
    airvent_y1 = 2.0 * obstacle_height / 4.0
    airvent_y2 = 3.0 * obstacle_height / 4.0
else:
    air_vent = False

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
nDTout = int(round(T / dt_fixed))

##########################################
#              Mesh & Domain             #
##########################################

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #

if air_vent:
    weir = [[[obstacle_x_start, 0], [obstacle_x_start, obstacle_height],
             [obstacle_x_end, obstacle_height], [obstacle_x_end, airvent_y2],
             [obstacle_x_end, airvent_y1], [obstacle_x_end, 0]]]
    vent = {'airvent': [[obstacle_x_end, airvent_y2]]}
else:
    weir = [[[obstacle_x_start, 0], [obstacle_x_start, obstacle_height],
             [obstacle_x_end, obstacle_height], [obstacle_x_end, 0]]]
    vent = None

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=weir,
                              special_boundaries=vent)

# ----- WAVES ----- #
if opts.waves:

    wave = wt.MonochromaticWaves(
        period = 2,
        waveHeight =0.018,
        mwl = waterLine_z,
        depth = waterLine_z,
        g = np.array(g),
        waveDir = (1.,0.,0.),
        wavelength = 0.5,
        meanVelocity = np.array([inflow_velocity, 0., 0.])
    )
    tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])

    tank.setGenerationZones(x_n=True, waves=wave)
    tank.setAbsorptionZones(x_p=True)

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

point_gauge_locations = []
column_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:

    number_of_gauges = tank_dim[0] / opts.gauge_dx + 1

    for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):
        if gauge_x < obstacle_x_start or gauge_x > obstacle_x_end:
            point_gauge_locations.append((gauge_x, opts.point_gauge_y, 0.),)
            column_gauge_locations.append(((gauge_x, 0., 0.),
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

# ----- EXTRA BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()

tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=outflow_level,
                                                    rhoUp=rho_1,
                                                    rhoDown=rho_0,
                                                    g=g,
                                                    refLevel=tank_dim[1])

if not opts.waves:
    tank.BC['x-'].setTwoPhaseVelocityInlet(U=[inflow_velocity,0.],
                                           waterLevel=waterLine_z)

# import pdb
# pdb.set_trace()
if air_vent:
    tank.BC['airvent'].p_dirichlet.uOfXT = lambda x, t: (tank_dim[1] - x[1]) \
                                                        * rho_1 * abs(g[1])
    tank.BC['airvent'].u_diffusive.uOfXT = lambda x, t: 0
    tank.BC['airvent'].v_dirichlet.uOfXT = lambda x, t: 0
    tank.BC['airvent'].v_diffusive.uOfXT = lambda x, t: 0
    tank.BC['airvent'].vof_dirichlet.uOfXT = lambda x, t: 1
    #tank.BC['airvent'].setTank()  #  unique boundary conditions in obstacle tanks don't have _b_or or _b_i setting yet (not sure how to pass that through intuitively) and thus cannot have setTank.  It could be set manually... but it's not really a big issue for the normal hydraulicStructures
    #[temp] check against report - different set of conditions than in the code, which might solve issues if issues need solving

# ----- MESH CONSTRUCTION ----- #

he = tank_dim[0] / float(4 * refinement - 1)
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
#            Signed Distance             #
##########################################


def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    phi_z_outflow = x[1] - outflow_level
    if phi_x <= 0.0:
        if phi_z < 0.0:
            return max(phi_x, phi_z)
        else:
            return phi_z
    else:
        if phi_z_outflow < 0.0:
            return phi_z_outflow
        else:
            if phi_z < 0.0:
                return min(phi_x, phi_z_outflow)
            else:
                return min(sqrt(phi_x ** 2 + phi_z ** 2), phi_z_outflow)