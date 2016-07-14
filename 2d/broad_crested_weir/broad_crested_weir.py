"""
A Broad Crested Weir
"""
import numpy as np
from math import *
from collections import namedtuple
import proteus.MeshTools
from proteus import Domain, WaveTools
from proteus.mprans import SpatialTools
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Context
from weir_tank import ShapeWeirTank2D

opts = Context.Options([
    # options
    ("waves", False, "Generate waves"),
    ("open_top", True, ""),
    # [temp] I think this is straightforward, but I'm not sure if it does apply that way
    ("air_vent", False, ""),
    # [temp] sends to AV conditions? But also the whole airVent bit?
    ("sponge_layers", False, ""),  # [temp] sets up the two sponge layers
    ("variable_mesh", False, ""),  # [temp] sends to VM conditions?  Probably could be reduced to refinement_x_borders/refinement_levels (as long as you have a good check for similar purpose, so users can't just put in one)
    ("refinement_level", 40, ""),
    ("generate_mesh", True, "Generate new mesh"),
    ("cfl", 0.9, "Target CFL"),
    # [temp] runCFL. I don't know what this means by target, but other context options label it as such (and take it to runCFL)
    # refinement
    ("refinement_x_borders", None, "list of cuts for variable meshes"),
    # [temp] I don't like this description. Get a better one.
    ("refinement_levels", None,
     "list of refinement levels for each variable mesh zone"),
    # dimensions & time stepping
    ("tank_dim", (3.5, 0.7), "(width,height) of the tank"),
    ("obstacle_dim", (0.5, 0.401), "(width,height) of the weir/ obstacle"),
    ("obstacle_start", 2.0, "x coordinate of start of obstacle"),
    ("sim_time", 10.0, "Simulation time"),
    ("dt_fixed", 0.02, ""),
    ("min_init_dt", 0.001, "Minimum initial time step (otherwise dt_fixed/10)"),
    # water
    ("init_water_height", 0.462,
     "initial inflow height of water (or mean height, for wave problems)"),
    ("init_water_over_obst", 0.1,
     "initial x overflow of water past obstacle_start"),
    ("init_water_outflow_height", -1,
     "initial outflow water height (or mean height, for wave problems)"),
    ("inflow_velocity", 0.047, ""),
    # [temp] I have a basic idea, but not sure I can fully write up the explanation
    ("outflow_velocity", 0, "")
    # [temp] unsure about this one's use. It's only an AV thing.
])

# ----------------------------------------------------
# discretization
# ----------------------------------------------------

# discretization input options
refinement = opts.refinement_level
variable_mesh = opts.variable_mesh
genMesh = opts.generate_mesh
useOldPETSc = False
useSuperlu = False
spaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = 0  # 0 -- None, 1 -- K-Epsilon, 2 -- K-Omega
timeDiscretization = 'be'  # 'vbdf','be','flcbdf'
applyCorrection = True
applyRedistancing = True
movingDomain = False

# discretization input checks
if spaceOrder not in [1, 2]:
    print "INVALID: spaceOrder" + str(spaceOrder)
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + str(useRBLES)
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

# discretization
nd = 2
if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 4)

# turbulence
# (1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega)
# [temp] Low priority, but it would be nice if RANS/turbulence was more readable
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4
else:
    ns_closure = 2

# parallel
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

# sponge layers
sponge_layers = opts.sponge_layers
#[temp] I currently believe the only references to this are sponge_layers itself, the two zone lengths, and the porosity/dragAlpha/dragBeta - the centers, epsFact, and so on are no longer part of this code.
if sponge_layers:
    GenerationZoneLength = 0.5  # [temp] candidates for being moved to Context or somewhere else more visible, if important and changing
    AbsorptionZoneLength = 0.5
    # generation zone
    left_sponge = GenerationZoneLength
    xRelaxCenter = left_sponge / 2.0
    epsFact_solid = GenerationZoneLength / 2.0
    # absorption zone
    right_sponge = opts.tank_dim[0] - AbsorptionZoneLength
    xRelaxCenter_2 = 0.5 * (right_sponge + opts.tank_dim[0])
    epsFact_solid_2 = AbsorptionZoneLength / 2.0
    #[temp] these should be equivalent to the arrays below, but let's track it carefully
    porosity = 1.0
    dragAlpha = 0.5 / 1.004e-6
    dragBeta = 0.0
    #[temp] end of new values
    porosityTypes = np.array([1.0,
                              1.0,
                              1.0,
                              1.0,
                              1.0,
                              1.0])
    dragAlphaTypes = np.array([0.0,
                               0.5 / 1.004e-6,
                               0.5 / 1.004e-6,
                               0.0,
                               0.0,
                               0.0])
    dragBetaTypes = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    epsFact_solidTypes = np.array(
        [0.0, epsFact_solid, epsFact_solid_2, 0.0, 0.0, 0.0])

# ----------------------------------------------------
# physical properties
# ----------------------------------------------------
# wind
if opts.waves:
    windVelocity = (0.0, 0.0)

# water
rho_0 = 998.2
nu_0 = 1.004e-6

# air
rho_1 = 1.205
nu_1 = 1.500e-5

# surface tension
sigma_01 = 0.0

# gravity
g = [0.0, -9.8]

# inflow/outflow
inflow_velocity = opts.inflow_velocity
outflow_velocity = opts.outflow_velocity

air_vent = opts.air_vent

# ----------------------------------------------------
# domain & mesh details
# ----------------------------------------------------

# tank
L = opts.tank_dim

# weir / obstacle
(obst_width, obst_height) = opts.obstacle_dim
obst_x_start = opts.obstacle_start
obst_x_end = obst_x_start + obst_width
obstacle = (obst_x_start, obst_height, obst_x_end)
obstacle_intersects = (obst_x_start, obst_x_end)
obstacle_points = ((obst_x_start, obst_height), (obst_x_end, obst_height))

# waterline
waterLine_x = obst_x_start + opts.init_water_over_obst
waterLine_z = opts.init_water_height
waterLine_z_outflow = opts.init_water_outflow_height
if waterLine_z_outflow < 0.0:
    waterLine_z_outflow = -(L[0] ** 2) - (L[1] ** 2)  # this guarantees it will always be bigger than any distance inside the tank
inflowHeightMean = waterLine_z
outflowHeightMean = waterLine_z_outflow

# weir/water size sanity checks
if not 0 <= obst_x_start <= L[0]:
    raise ValueError("OUT OF BOUNDS: Obstacle starts outside of tank.")
elif not 0 <= obst_x_end <= L[0]:
    raise ValueError("OUT OF BOUNDS: Obstacle ends outside of tank.")
if obst_height > L[1]:
    raise ValueError("OUT OF BOUNDS: Obstacle is taller than tank.")
if not 0 <= waterLine_x <= L[0]:
    raise ValueError("OUT OF BOUNDS: Water ends outside of tank.")
if waterLine_z > L[1]:
    raise ValueError("OUT OF BOUNDS: Inflow water is taller than tank.")
if waterLine_z_outflow > L[1]:
    raise ValueError("OUT OF BOUNDS: Outflow water is taller than tank.")

# air vent #[temp] low priority, but maybe could be neater/more explanatory
airvent_y1 = 2.5 * obst_height / 4.0
airvent_y2 = 3.5 * obst_height / 4.0
if air_vent:
    yp_airvent_point = (obst_x_end, airvent_y2)
    yn_airvent_point = (obst_x_end, airvent_y1)
    obstacle_points += (yp_airvent_point, yn_airvent_point)

# refinement
he = L[0] / float(4 * refinement - 1)
if variable_mesh:
    x_refine = opts.refinement_x_borders
    x_refine_level = opts.refinement_levels
    weak_bc_penalty_constant = 100.0
else:
    he *= 0.5

# mesh inputs
structured = False

# ----------------------------------------------------
# gauges
# ----------------------------------------------------
#
# # ##### Still a bit of a hack job, but based on the Proteus Gauges.
# #       The goal will be to pull these options to the context and
# #       let the different types of gauges be generated based on
# #       individual problem setup.
# point_gauges = PointGauges(gauges=((('p', 'u', 'v'), ((0.05, 0.65, 0.0),)),),
#                            activeTime=None,
#                            sampleRate=0,
#                            fileName='point_gauge_1.csv')
# # LineGauges(gauges=((('u0',),(((0,0,0),(1,1,1)),)),),
# line_gauges = LineGauges(
#     gauges=((('p', 'u', 'v'), (((3.4, 0.0, 0.0), (3.4, 0.1, 0.0)),)),),
#     activeTime=None,
#     sampleRate=0,
#     fileName='line_gauge_1.csv')
# line_gauges_phi = LineGauges(
#     gauges=((('phid',), (((3.4, 0.0, 0.0), (3.4, 0.1, 0.0)),)),),
#     activeTime=None,
#     sampleRate=0,
#     fileName='line_gauge_1_phi.csv')
# # #### What is/where is phi? Gauges can't find it.  For now we're switching to 'phid' as an attribute actually considered by redist and all.  However, this can easily be changed later on as the system is refit - this is just testing that the basic syntax can be fit into the code like this

# ----------------------------------------------------
# mesh
# ----------------------------------------------------

#[temp]
generating_waves = #needs to be set [the waves for AV]
trivial_waves = #also needs to be set [constant flow, if possible]

if useHex:
    nny = ceil(L[1]/he) + 1
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    if structured:
        nnx = ceil(L[0]/he) + 1
        nny = ceil(L[1]/he) + 1
        domain = Domain.RectangularDomain(L)
    else:
        domain = Domain.PlanarStraightLineGraphDomain()

        tank = ShapeWeirTank2D(domain = domain,
                               dim = L,
                               obstacle_intersects = obstacle_intersects,
                               obstacle_points = obstacle_points,
                               floating_obstacles = None,
                               floating_centers = None,
                               special_boundaries = {[yp_airvent_point,
                                                      yn_airvent_point]: 'airvent'},
                               coords = None,
                               from_0 = True)
        tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave = trivial_waves,
                                                       wind_speed = windVelocity)
        tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel = waterLine_z,
                                                            rhoUp = rho_1,
                                                            rhoDown = rho_0,
                                                            g = g,
                                                            refLevel= L[1])
        tank.setSponge(left=GenerationZoneLength,
                       right=AbsorptionZoneLength)
        tank.setGenerationZones(waves=generating_waves,
                                windSpeed=windVelocity,
                                left=True,
                                dragAlpha=dragAlpha,
                                dragBeta=dragBeta,
                                porosity=porosity)
        tank.setAbsorptionZones(right=True,
                                dragAlpha=dragAlpha,
                                dragBeta=dragBeta,
                                porosity=porosity) #[temp] should porosity and drags be zone specific (different between gen and abs)?
        tank.attachPointGauges('twp',
                               gauges=((('p', 'u', 'v'), ((0.05, 0.65, 0.0),)),),
                               activeTime=None,
                               sampleRate=0,
                               fileName='point_gauge_1.csv')
        tank.attachLineGauges('twp',
                              gauges=((('p', 'u', 'v'), (((3.4, 0.0, 0.0), (3.4, 0.1, 0.0)),)),),
                              activeTime=None,
                              sampleRate=0,
                              fileName='line_gauge_1.csv')
        tank.attachLineGauges('redist',
                              gauges=((('phid',), (((3.4, 0.0, 0.0), (3.4, 0.1, 0.0)),)),),
                              activeTime=None,
                              sampleRate=0,
                              fileName='line_gauge_1_phi.csv')
        tank.setHorizontalVariableMesh(region_boundaries = x_refine)
        #[temp] here we will set the airvent condition.  Each variable change that we need can be made like so:
        #[temp] tank.BC['airvent'].u_diffusive.uOfXT = lambda x, t: 0.
        #[temp] (with different names and different lambda function, of course)
        tank.BC['airvent'].p_dirichlet.uOfXT = lambda x, t: (L[1]-x[1])*rho_1*abs(g[1])
        ## actually, outflow pressure, which is  (L[1]-x[1])*rho_1*abs(g[1])
        tank.BC['airvent'].u_diffusive.uOfXT = lambda x, t: 0
        tank.BC['airvent'].v_dirichlet.uOfXT = lambda x, t: 0
        tank.BC['airvent'].v_diffusive.uOfXT = lambda x, t: 0
        tank.BC['airvent'].vof_dirichlet.uOfXT = lambda x, t: 1
        tank.BC['airvent'].setTank()
        #[temp] end of airvent conditions
        tank.constructShape()
        SpatialTools.assembleDomain(domain)
        trigArea = 0.5 * he ^ 2
        domain.regionConstraints = map(lambda level: trigArea * level ** 2,
                                       x_refine_level)

# ----------------------------------------------------
# time stepping
# ----------------------------------------------------
T = opts.sim_time
dt_fixed = opts.dt_fixed
dt_init = min(0.1 * dt_fixed, opts.min_init_dt)
runCFL = opts.cfl
nDTout = int(round(T / dt_fixed))

# ----------------------------------------------------
# numerical parameters
# ----------------------------------------------------

backgroundDiffusionFactor = 0.01  # [temp] move this somewhere (it's currently the default value - but may be useful to reset for people)
ns_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor = 0.5  # [temp] AV = 0.75 (BASE and VM match, but AV differs on a few, noted here)
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.25  # [temp] AV = 0.75
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0  # [temp] AV = 1.5
    vof_shockCapturingFactor = 0.25  # [temp] AV = 0.75
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0  # [temp] AV = 1.5
    rd_shockCapturingFactor = 0.25  # [temp] AV = 0.75
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_consrv_dirac \
        = epsFact_viscosity \
        = epsFact_curvature \
        = epsFact_density \
        = epsFact_consrv_heaviside \
        = epsFact_vof \
        = epsFact_density \
        = 3.0
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True
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
    epsFact_density = epsFact_consrv_dirac = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_density
    epsFact_consrv_heaviside = epsFact_vof = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

# tolerances
# [temp] I'd like these to be more readable.  Clearer names, if possible,
# [temp] -as well as nicer looking max statements that help us figure it out more
# [temp] -but this isn't high priority.
ns_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.005 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)


# ----------------------------------------------------
# functions & misc.
# ----------------------------------------------------


def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    phi_z_outflow = x[1] - waterLine_z_outflow
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


def outflowPressure(x, t):
    return (L[1] - x[1]) * rho_1 * abs(g[1]) #[temp] I might have taken the only use of it above, and thus made this declaration irrelevant


# solution variables
# [temp] AV only functions.  That was _AV's name for the set.

def wavePhi(x, t):
    return x[1] - inflowHeightMean


def waveVF(x, t):
    return smoothedHeaviside(epsFact_consrv_heaviside * he, wavePhi(x, t))


def twpflowVelocity_u(x, t):
    waterspeed = inflow_velocity
    H = smoothedHeaviside(epsFact_consrv_heaviside * he, wavePhi(x, t)
                          - epsFact_consrv_heaviside * he)
    u = H * windVelocity[0] + (1.0 - H) * waterspeed
    return u


def twpflowVelocity_v(x, t):
    waterspeed = 0.0
    H = smoothedHeaviside(epsFact_consrv_heaviside * he, wavePhi(x, t)
                          - epsFact_consrv_heaviside * he)
    return H * windVelocity[1] + (1.0 - H) * waterspeed


def twpflowVelocity_w(x, t):
    return 0.0


def twpflowFlux(x, t):
    return -twpflowVelocity_u(x, t)


def outflowPhi(x, t):
    return x[1] - outflowHeightMean


def outflowVF(x, t):
    return smoothedHeaviside(epsFact_consrv_heaviside * he, outflowPhi(x, t))


def outflowVel(x, t):
    waterspeed = outflow_velocity
    H = smoothedHeaviside(epsFact_consrv_heaviside * he, outflowPhi(x, t)
                          - epsFact_consrv_heaviside * he)
    u = (1.0 - H) * waterspeed
    return u


def zeroVel(x, t):
    return 0.0

#[temp] deprecated wave stuff in comment
# relaxation zone

# RelaxationZone = namedtuple("RelaxationZone", "center_x sign u v w")
#
#
# class RelaxationZoneWaveGenerator(AV_base):
#     """ Prescribe a velocity penalty scaling in a material zone via a Darcy-Forchheimer penalty
#
#     :param zones: A dictionary mapping integer material types to Zones, where a Zone is a named tuple
#     specifying the x coordinate of the zone center and the velocity components
#     """
#
#     def __init__(self, zones):
#         assert isinstance(zones, dict)
#         self.zones = zones
#
#     def calculate(self):
#         for l, m in enumerate(self.model.levelModelList):
#             for eN in range(m.coefficients.q_phi.shape[0]):
#                 mType = m.mesh.elementMaterialTypes[eN]
#                 if self.zones.has_key(mType):
#                     for k in range(m.coefficients.q_phi.shape[1]):
#                         t = m.timeIntegration.t
#                         x = m.q['x'][eN, k]
#                         m.coefficients.q_phi_solid[eN, k] = self.zones[
#                                                                 mType].sign * (
#                                                             self.zones[
#                                                                 mType].center_x -
#                                                             x[0])
#                         m.coefficients.q_velocity_solid[eN, k, 0] = self.zones[
#                             mType].u(x, t)
#                         m.coefficients.q_velocity_solid[eN, k, 1] = self.zones[
#                             mType].v(x, t)
#                         # m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
#         m.q['phi_solid'] = m.coefficients.q_phi_solid
#         m.q['velocity_solid'] = m.coefficients.q_velocity_solid
#
#
# if sponge_layers:  # [----] should be set to first and last (or pref, sponge layer specific) region flags, rather than hardcoded 1 and 4
#     rzWaveGenerator = RelaxationZoneWaveGenerator(
#         zones={1: RelaxationZone(xRelaxCenter,
#                                  -1.0,
#                                  twpflowVelocity_u,
#                                  twpflowVelocity_v,
#                                  twpflowVelocity_w),
#                4: RelaxationZone(xRelaxCenter_2,
#                                  1.0,
#                                  outflowVel,
#                                  zeroVel,
#                                  zeroVel)})
