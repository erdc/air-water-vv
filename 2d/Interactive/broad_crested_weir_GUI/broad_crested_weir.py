"""
A Broad Crested Weir
"""
import numpy as np
from math import sqrt
import proteus.MeshTools
from proteus import (Context,
                     Domain,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside


opts = Context.Options([
    # test options
    ("waves", False, "Generate waves - uses sponge layers."),
    ("air_vent", False, "Include an air vent in the obstacle."),
    # air vent position
    ("airvent_y1",0.25,"Vertical distance from bottom to the air ventilation boundary patch"),
    ("airvent_dim",0.05,"Dimension of the air boundary patch"),
    # water
    ("water_level", 0.54, "Height of (mean) free surface above bottom"),
    ("water_width_over_obst",1.02, "Initial width of free surface relative to"
                                   " the obstacle location"),
    ("outflow_level", 0.04, "Height of (mean) free surface of water outflow "
                          "give a negative number if no initial outflow."),
    ("inflow_velocity", 0.139, "Wave or steady water inflow velocity"),
    ("outflow_velocity", 0.0, "Initial wave or steady water outflow velocity"),
    # tank
    ("tank_dim", (2.5, 1.0), "Dimensions (x,y) of the tank"),
    ("absorption", True, "Use generation/absorption zones"),
    ("tank_sponge", (0.5,0.5), "Length of (generation, absorption) zones, if any"),
    ("obstacle_dim", (0.5, 0.401), "Dimensions (x,y) of the obstacle."),
    ("obstacle_x_start", 1.0, "x coordinate of the start of the obstacle"),
    # gauges
    ("gauge_output", True, "Produce gauge data"),
    # refinement
    ("refinement", 10, "Refinement level"),
    ("he", 0.03, "Maximum size"),
    ("cfl", 0.75, "Target cfl"),
    ("variable_refine_borders", None, "List of vertical borders between "
                                    "refinement regions (include 0 and "
                                    "tank_dim[0] if you add sponge layers "
                                    "and want to differentiate them)"),
    ("variable_refine_levels", None, "List of refinement levels in each region"
                                   " (should have 1 more value than "
                                   "variable_refine_borders as a result)."),
    # run time
    ("T", 5, "Simulation time"),
    ("dt_fixed", 0.025, "Fixed time step"),
    ("dt_init", 0.0001, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")],
    mutable=True)


def Update_Model():
    global he
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
    tank_sponge = opts.tank_sponge
    obstacle_dim = opts.obstacle_dim
    obstacle_x_start = opts.obstacle_x_start
    obstacle_x_end = obstacle_x_start + obstacle_dim[0]
    obstacle_height = obstacle_dim[1]

    # air vent
    if opts.air_vent:
        air_vent = True
        airvent_y1 = opts.airvent_y1
        airvent_y2 = airvent_y1 + opts.airvent_dim
    else:
        air_vent = False

    # sanity checks
    if waterLine_z > tank_dim[1]:
        raise ValueError("ERROR: Water (level: %s) overflows height of tank (%s)"
                         % (waterLine_z, tank_dim[1]))
    if outflow_level > tank_dim[1]:
        raise ValueError("ERROR: Water (outflow level: %s) overflows height of tank (%s)"
                         % (outflow_level, tank_dim[1]))
    if obstacle_x_end > tank_dim[0] or obstacle_height > tank_dim[1]:
        raise ValueError("ERROR: Obstacle (height: %s, width: %s, start: %s) lies "
                         " outside of tank (height: %s, width: %s)"
                         % (obstacle_dim[1], obstacle_dim[0], obstacle_x_start,
                            tank_dim[1], tank_dim[0]))
    if waterLine_x + obstacle_dim[0] > tank_dim[0]:
        raise ValueError("ERROR: Water starts outside of tank at x = %s (tank: %s)"
                         % (waterLine_x+obstacle_dim[0], tank_dim[0]))
    if opts.air_vent:
        if airvent_y2 > obstacle_height:
            raise ValueError("ERROR: Air ventilation (%s) exceeds the obstacle (%s)"
                             % (airvent_y2, obstacle_height))

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
    useSuperlu = not opts.parallel
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
    omega = 1.
    if opts.waves:
        omega=2*np.pi/2.
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

    if opts.absorption:
        dragAlpha = 5.*omega/nu_0
        #tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])
        #tank.setAbsorptionZones(x_n=True, dragAlpha = dragAlpha)
        tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])
        tank.setAbsorptionZones(x_p=True, dragAlpha = dragAlpha)
    
    # ------ he -------------#
    he=opts.he


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

    if opts.gauge_output:

        tank.attachLineGauges(
            'twp',
            gauges=((('p','u','v'), (((2.0, 0.0, 0.0),
                                      (2.0, 0.5, 0.0)),
                                     )),),
            activeTime = None,
            sampleRate = 0,
            fileName = 'p_u_gauges.csv'
        )



    # ----- EXTRA BOUNDARY CONDITIONS ----- #

    # Open Top
    tank.BC['y+'].setAtmosphere()

    # Free Slip Tank
    tank.BC['y-'].setFreeSlip()

    # Outflow
    tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel = 0.0,
                                                        rhoUp = rho_1,
                                                        rhoDown = rho_0,
                                                        g = g,
                                                        refLevel = tank_dim[1],
                                                        smoothing = 1.5 * he,
                                                        U = [0.0, 0., 0.],
                                                        Uwind = [0.0, 0.0] )



    if opts.absorption:
        tank.BC['sponge'].setNonMaterial()

    # Inflow / Sponge
    if not opts.waves:
        tank.BC['x-'].setTwoPhaseVelocityInlet(U=[inflow_velocity,0.,0.],
                                               waterLevel=waterLine_z,
                                               smoothing=3.0*he)
    if air_vent:
        tank.BC['airvent'].setHydrostaticPressureOutletWithDepth(seaLevel=0.0,
                                                                 rhoUp=rho_1,
                                                                 rhoDown=rho_0,
                                                                 g=g,
                                                                 refLevel=tank_dim[1],
                                                                 smoothing=1.5*he)

    # ----- MESH CONSTRUCTION ----- #
    triangleOptions="VApq30Dena%8.8f" % ((opts.he**2)/2.0,)
    parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
    nLayersOfOverlapForParallel = 0
    structured=False
    #   he = tank_dim[0] / float(4 * refinement - 1)
    domain.MeshOptions.he = he
#    #domain.MeshOptions.triangleOptions="VApq30Dena%8.8f" % ((opts.he**2)/2.0,)
#    domain.MeshOptions.triangleOptions="VApq30Dena%8.8f" % (opts.he)
#    domain.MeshOptions.parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
#    domain.MeshOptions.nLayersOfOverlapForParallel = 0
    st.assembleDomain(domain)

    ##########################################
    # Numerical Options and Other Parameters #
    ##########################################

    # ----- STRONG DIRICHLET ----- #

    ns_forceStrongDirichlet = False
    weak_bc_penalty_constant = 10.0/nu_0#Re

    # ----- NUMERICAL PARAMETERS ----- #

    if useMetrics:
        ns_shockCapturingFactor = 0.75
        ns_lag_shockCapturing = True
        ns_lag_subgridError = True
        ls_shockCapturingFactor = 0.75
        ls_lag_shockCapturing = True
        ls_sc_uref = 1.0
        ls_sc_beta = 1.50
        vof_shockCapturingFactor = 0.75
        vof_lag_shockCapturing = True
        vof_sc_uref = 1.0
        vof_sc_beta = 1.50
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
        ls_sc_beta = 1.5
        vof_shockCapturingFactor = 0.9
        vof_lag_shockCapturing = True
        vof_sc_uref = 1.0
        vof_sc_beta = 1.5
        rd_shockCapturingFactor = 0.9
        rd_lag_shockCapturing = False
        epsFact_density = epsFact_viscosity = epsFact_curvature \
            = epsFact_vof = ecH = epsFact_consrv_dirac \
            = 3.0
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

    ns_nl_atol_res = max(1.0e-10,0.001*he**2)
    vof_nl_atol_res = max(1.0e-10,0.001*he**2)
    ls_nl_atol_res = max(1.0e-10,0.001*he**2)
    rd_nl_atol_res = max(1.0e-10,0.005*he)
    mcorr_nl_atol_res = max(1.0e-10,0.001*he**2)
    kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
    dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

    # ----- TURBULENCE MODELS ----- #
    #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

    if useRANS == 1:
        ns_closure = 3
    elif useRANS == 2:
        ns_closure = 4
    else:
        ns_closure = 0
    
    
    domain.MeshOptions.he = he
    st.assembleDomain(domain)
    globals().update(locals())

##########################################
#            Signed Distance             #
##########################################

def wavePhi(x,t):
    return x[1] - waterLine_z

def outflowPhi(x,t):
    return x[1] - outflow_level

def twpflowVelocity_u(x,t):
    waterspeed = inflow_velocity
    H = smoothedHeaviside(ecH*he,wavePhi(x,t)-ecH*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_u_D(x, t):
    waterspeed = outflow_velocity
    H = smoothedHeaviside(ecH * he, outflowPhi(x, t) - ecH * he)
    u = H * windVelocity[0] + (1.0 - H) * waterspeed
    return u


Update_Model()
    
def signedDistance(x,
                   waterLine_x=waterLine_x,
                   waterLine_z=waterLine_z,
                   outflow_level=outflow_level):
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
