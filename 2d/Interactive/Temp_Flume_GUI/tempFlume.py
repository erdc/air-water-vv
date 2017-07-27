"""
A Broad Crested Weir
"""
import numpy as np

from proteus import (Domain, Context,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans import BodyDynamics as bd
from math import *
import pandas as pd


opts = Context.Options([
    # test options
    ("waves", False, "Generate waves - uses sponge layers."),
    ("waveHeight", 0.018 , "Height of waves over water level"),
    # water
    ("water_level", 0.36576, "Height of (mean) free surface above bottom"),
     # caisson
    ("caisson", True, "caisson"),
    ("caisson_scale", 21, "Caisson scale ratio 1:scale"),
    ("caisson_mass", 1000, "Caisson original mass"),
    ("caisson_coords",None, "Coord of the caisson"),
    ("caisson_Yoffset",-0.005524, "Vertical offset from the original possition of the caisson"),
    ("caisson_width", 1., "Width of the caisson"),
    ("caisson_BC", 'noslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (0.0, 0.0, 0.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 0.0), "Rotational DOFs"),
    ("VCG", None, "vertical position of the barycenter of the caisson"),
    ('scheme', 'Forward_Euler', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Forward_Euler)'),
    ("springs", not True, "Switch on/off soil module"),
    ("Kx", 0.0, "Horizontal stiffness in N/m"),
    ("Ky", 0.0, "Vertical stiffness in N/m"),
    ("Krot", 0.0, "Rotational stiffness in N*m"),
    ("Cx", 0.0, "Damping factor in N/m*s "),
    ("Cy", 0.0, "Damping factor in N/m*s "),
    ("Crot", 0.0, "Rotational damping factor in N*m*s "),
#    ("caisson_mass", 125., "Mass of the caisson"),
#    ("caisson_inertia", 4.05, "Inertia of the caisson"),
    ("rotation_angle", 0 ,"Initial rotation angle (in degrees)"),
    #Initial Velocities
    ("inflow_velocity", 0.27432, "Wave or steady water inflow velocity"),
    ("outflow_velocity", 0.27432, "Initial wave or steady water outflow velocity"),
    ("wind_velocity", (0.,0.), "Initial wind velocity"),
    # Turbulence
    ("useRANS", 0, "Switch ON turbulence models: 0-None, 1-K-Epsilon, 2-K-Omega1998, 3-K-Omega1988"),
    ("ns_closure", 4, "1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega"),
    ("c_mu", 0.09, "mu coefficient for the turbulence model"),
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("sigma_e", 1.0, "sigma_e coefficient for the turbulence model"),
    ("Y", 0.01, "Y used for y+ calculation"),
    # tank
    ("tank_dim", (1.5,0.5), "Dimensions (x,y) of the tank"),
    ("absorption", True, "Use generation/absorption zones"),
    ("tank_sponge", (0.25,0.25), "Length of (generation, absorption) zones, if any"),
    # gauges
    ("gauge_output", True, "Produce gauge data"),
    # refinement
    ("refinement", 1, "Refinement level"),
    ("he", 0.05, "Set characteristic element size"),
    ("he_caisson", 0.05, "Set maximum characteristic element size on caisson boundary"),
    ("cfl", 0.85, "Target cfl"),
    ("variable_refine_borders", None, "List of vertical borders between "
                                    "refinement regions (include 0 and "
                                    "tank_dim[0] if you add sponge layers "
                                    "and want to differentiate them)"),
    ("variable_refine_levels", None, "List of refinement levels in each region"
                                   " (should have 1 more value than "
                                   "variable_refine_borders as a result)."),
    # run time
    ("T", 1, "Simulation time"),
    ("dt_fixed", 0.001, "Fixed time step"),
    ("dt_init", 0.0001, "Minimum initial time step (otherwise dt_fixed/10)"),
    # run details
    ('movingDomain', False, "Moving domain and mesh option"),
    ("gen_mesh", True, "Generate new mesh"),
    ("parallel", True, "Run in parallel")],
    mutable=True)


def Update_Model():

    # flow
    inflow_velocity = opts.inflow_velocity
    outflow_velocity = opts.outflow_velocity

    # tank
    tank_dim = opts.tank_dim
    tank_sponge = (opts.tank_dim[0]/10,opts.tank_dim[1]/10)

    # Initial condition
    waterLine_x = 2*tank_dim[0]
    waterLine_z = opts.water_level


    # sanity checks
    if waterLine_z > tank_dim[1]:
        raise ValueError("ERROR: Water (level: %s) overflows height of tank (%s)"
                         % (waterLine_z, tank_dim[1]))

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
    movingDomain = opts.movingDomain
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
    nu_1 = 2*1.500e-5

    rho_b=(rho_0+rho_1)/2

    # Surface Tension
    sigma_01 = 0.0

    # Gravity
    g = [0., -9.81, 0.]

    # wind
    windVelocity = (opts.wind_velocity[0],opts.wind_velocity[1], 0.)

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


    tank = st.Tank2D(domain=domain,dim=tank_dim)

        #------ CAISSON DEFINITION -----#
    rotation_angle = np.radians(opts.rotation_angle)

    if opts.caisson is True:
        # Read caisson options
        scale = 1/float(opts.caisson_scale)
        mass=scale*opts.caisson_mass
        VCG = opts.VCG
        rotation = np.radians(opts.rotation_angle)
        caisson_width  = opts.caisson_width
        he_min=opts.he_caisson
        if opts.caisson_coords is None:
            caisson_coords = [tank_dim[0]/2., waterLine_z+opts.caisson_Yoffset]
        else:
            caisson_coords = opts.caisson_coords

        if opts.he_caisson:
            he_caisson = opts.he_caisson
        else:
            he_caisson = opts.he

        # Read Caisson dimensions

        xl = pd.ExcelFile("Ribbon Bridge Section Line Segment.xlsx")
        df = xl.parse(0)
        x=np.asarray(df['x'].tolist())
        y=np.asarray(df['y'].tolist())
        # Offset to (0,0)
        x-=0.5*(max(x)+min(x))
        y-=min(y)

        #Real dimensions compensation (Need to correct real values on excel)
        prescale=0.6459851273
        x=x*prescale
        y=y*prescale
        # Scale to experimental model
        xp=x*scale
        yp=y*scale




        dim=[max(xp)-min(xp),max(yp)-min(yp)]

        # Interpolate points to he_caisson
        xd =np.diff(xp)
        yd = np.diff(yp)
        dist = np.sqrt(xd**2+yd**2)
        u = np.cumsum(dist)
        u = np.hstack([[0],u])

        t = np.linspace(0,u.max(),int(u.max()/he_min))
        xBase = np.interp(t, u, xp)
        yBase = np.interp(t, u, yp)

        xTop = np.linspace(xp.min() + he_min, xp.max(), int((xp.max() - xp.min()) / he_min), endpoint=False)
        yTop=np.full(len(xTop),yp[-1])
        x=np.hstack([xBase,xTop])
        y=np.hstack([yBase,yTop])


        # Calculate Parameters

        if VCG is None:
            VCG = dim[1]/2.
        I = ((dim[0]**2)+(dim[1]**2))*mass/12.
        barycenter = (0, 0, 0.)

        vertices = []
        vertexFlags = []
        segments = []
        segmentFlags = []

        v_start = 0
        for i in range(len(x)):
            v_start = len(vertices)
            v = [[x[i],y[i]]]
            if v_start >= 1:
                s = [[v_start-1, v_start]]
            else:
                s = []
            vertices += v
            vertexFlags += [1]*len(v)
            segments += s+[[len(vertices)-1, len(vertices)]]
            segmentFlags += [1]*len(s)+[1]
        segments[-1][1] = 0  # last segment links to vertex 0
        boundaryTags = {'caisson': 1}
        caisson = st.CustomShape(domain, barycenter=barycenter,
                                vertices=vertices, vertexFlags=vertexFlags,
                                segments=segments, segmentFlags=segmentFlags,
                                boundaryTags=boundaryTags)
        facet = []
        for i, vert in enumerate(caisson.vertices):
            facet += [i]
        caisson.facets = np.array([[facet]])
        caisson.facetFlags = np.array([1])
        caisson.regionFlags = np.array([1])

        ang = rotation_angle
        caisson.setHoles([[0., dim[1]/2]])
        caisson.holes_ind = np.array([0])
        caisson.translate([caisson_coords[0], caisson_coords[1]])
        caisson.rotate(ang, pivot=caisson.barycenter)
        for bc in caisson.BC_list:
            if opts.caisson_BC == 'noslip':
                bc.setNoSlip()
            if opts.caisson_BC == 'freeslip':
                bc.setFreeSlip()

    ## --- Body properties setup
    if opts.movingDomain:
        bridge2D = bd.RigidBody(shape=caisson)
        free_x = opts.free_x # Translational DOFs
        free_r = opts.free_r  # Rotational DOFs
        bridge2D.setMass(mass)
        bridge2D.It= I/bridge2D.mass/caisson_width
        bridge2D.setConstraints(free_x=free_x, free_r=free_r)
        bridge2D.setNumericalScheme(scheme=opts.scheme)
        bridge2D.setRecordValues(filename='bridge2D', all_values=True)

    #----- Spring setup
        if opts.springs:
            Kx = opts.Kx
            Ky = opts.Ky
            Krot = opts.Krot
            Cx = opts.Cx
            Cy = opts.Cy
            Crot = opts.Crot
            bridge2D.setSprings(opts.springs, Kx, Ky, Krot, Cx, Cy, Crot)
###########################################################################################################################################################################
# ----- Turbulence ----- #
###########################################################################################################################################################################

    d = dim[0]*scale #maxLength of Caisson
    L = 0. #tank_dim[0]
    U0 = inflow_velocity # freestream velocity
    c_mu = opts.c_mu
    Re0 = U0*d/nu_0
    ReL = U0*L/nu_0

    Uwind = windVelocity # opts.meanVelocity
    dwind = tank_dim[1] - opts.water_level #d
    Re0wind = Uwind[0]*dwind/nu_1

    # Schlichting
    cf = 0.045*(Re0**(-1./4.))
    cfwind = 0.
    if Re0wind >0.: cfwind = 0.045*(Re0wind**(-1./4.))

    Y_ = opts.Y

    Ut = U0*sqrt(cf/2.)
    kappaP = (Ut**2)/sqrt(c_mu)
    dissipationP = (Ut**3)/(0.41*Y_)

    # air
    Utwind = Uwind[0]*sqrt(cfwind/2.)
    kappaPwind = (Utwind**2)/sqrt(c_mu)
    if kappaPwind == 0.: kappaPwind = 1e-10
    dissipationPwind = (Utwind**3)/(0.41*Y_)
    if dissipationPwind == 0.: dissipationPwind = 1e-10

    # ke or kw
    useRANS = opts.useRANS  # 0 -- None
                            # 1 -- K-Epsilon
                            # 2 -- K-Omega, 1998
                            # 3 -- K-Omega, 1988
    if useRANS >= 2:
        # in kw model w = e/k
        dissipationP = dissipationP/kappaP
        dissipationP = dissipationPwind/kappaPwind

    # inlet values
    kInflow = kappaP #* 0.0001 # None
    dissipationInflow = dissipationP #* 0.0001 # None

    kInflowAir = kappaPwind #* 0.0001 # None
    dissipationInflowAir = dissipationPwind #* 0.0001 # None

#############################################################################################################################################################################################################################################################################################################################################################################################            ##############################################################################################################################################################################################################
    # ----- WAVES ----- #
    omega = 1.
    if opts.waves:
        omega=2*np.pi/2.
        wave = wt.MonochromaticWaves(
            period = 2,
            waveHeight =opts.waveHeight,
            mwl = waterLine_z,
            depth = waterLine_z,
            g = np.array(g),
            waveDir = (1.,0.,0.),
            wavelength = 0.5,
            meanVelocity = np.array([inflow_velocity, 0., 0.])
        )

    if opts.absorption:
        dragAlpha = 5.*omega/nu_0
        tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])
        #tank.setSponge(x_n = opts.tank_sponge[0])
        tank.setAbsorptionZones(x_n=True, dragAlpha = dragAlpha)
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

#    if opts.gauge_output:
#
#        tank.attachLineGauges(
#            'twp',
#            gauges=((('p','u','v'), (((2.0, 0.0, 0.0),
#                                      (2.0, 0.5, 0.0)),
#                                     )),),
#            activeTime = None,
#            sampleRate = 0,
#            fileName = 'p_u_gauges.csv'
#        )



    # ----- EXTRA BOUNDARY CONDITIONS ----- #

    # Open Top
    #tank.BC['y+'].setAtmosphere()
    tank.BC['y+'].setFreeSlip()

    # Free Slip Tank
    tank.BC['y-'].setFreeSlip()

    # Outflow
    tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=waterLine_z,
                                                        rhoUp=rho_1,
                                                        rhoDown=rho_0,
                                                        g=g,
                                                        refLevel=tank_dim[1],
                                                        smoothing=1.0*he,
                                                        U=[inflow_velocity,0.,0.],
                                                        Uwind=Uwind )

    if opts.absorption:
        tank.BC['sponge'].setNonMaterial()

    # Inflow / Sponge
    if not opts.waves:
    #    tank.BC['x-'].setTwoPhaseVelocityInlet(U=[inflow_velocity,0.,0.],
    #                                           waterLevel=waterLine_z,
    #                                           smoothing=3.0*he)
        tank.BC['x-'].setTwoPhaseVelocityInlet(U=[inflow_velocity,0.,0.],
                                               waterLevel=waterLine_z,
                                               smoothing=he*3.0,
                                               Uwind=Uwind,
                                               kInflow=kInflow,
                                               dissipationInflow=dissipationInflow,
                                               kInflowAir=kInflowAir,
                                               dissipationInflowAir=dissipationInflowAir)

    if opts.caisson:
        # let gmsh know that the caisson is IN the tank
        tank.setChildShape(caisson, 0)

    if opts.movingDomain==True:
        for tb in tank.BC_list:
            tb.setFixedNodes()
        for bc in caisson.BC_list:
            bc.hx_dirichlet.uOfXT, bc.hy_dirichlet.uOfXT, bc.hz_dirichlet.uOfXT = None, None, None
            bc.u_stress.uOfXT, bc.v_stress.uOfXT, bc.w_stress.uOfXT = None, None, None



    domain.MeshOptions.he = he
    st.assembleDomain(domain)


    triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
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
    mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
    # ----- TURBULENCE MODELS ----- #
    #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
    

    if useRANS == 1:
        ns_closure = 3
    elif useRANS == 2:
        ns_closure = 4
    else:
        ns_closure = 0
    
   
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
    #H = smoothedHeaviside(ecH * he, outflowPhi(x, t) - ecH * he)
    H = smoothedHeaviside(ecH * he, wavePhi(x, t) - ecH * he)
    u = H * windVelocity[0] + (1.0 - H) * waterspeed
    return u


Update_Model()
    
def signedDistance(x,waterLine_z=waterLine_z):
    phi_x = x[0]-waterLine_x
    phi_z = x[nd-1]-waterLine_z

    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)
