"""
A Broad Crested Weir
"""
import numpy as np
from math import sqrt
from proteus import (Domain, Context,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

opts = Context.Options([
    # physical parameters
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0, -9.805, 0]), "Gravity"),

    #waves and ventilation
    ("waves", False, "Generate waves - uses sponge layers."),
    ("air_vent", False, "Include an air vent in the obstacle."),

    # air vent position
    ("airvent_y1",0.25,"Vertical distance from bottom to the lower vertex of the air ventilation boundary in m"),
    ("airvent_dim",0.1,"Dimension of the air boundary patch in m"),

    # water
    ("mwl", 0.54, "Mean level at inflow in m"),
    ("waterLine_z", 0.54, "used to calculate signed distance)"),
    ("water_level", 0.54, "used to define pressure"),

    ("waterLine_x",1.02, "Domain length upstream of the obstacle in m"),
    ("outflow_level", 0.2, "Estimated mean water level at the outlet in m "),
    ("U", np.array([0.139,0.,0.]), "Water inflow velocity in m/s"),
    ("outflow_velocity", 3.0, "Estimated water outflow velocity in m/s"),

    # tank
    ("tank_dim", (6., 1.0), "Dimensions (x,y) of the tank in m"),
    ("tank_sponge", (2.,2.), "Length of (generation, absorption) zones in m, if any"),
    ("obstacle_dim", (0.5, 0.401), "Dimensions (x,y) of the obstacle. in m"),
    ("obstacle_x_start", 1.0, "x coordinate of the start of the obstacle in m"),

    # gauges
    ("gauge_output", True, "Produce gauge data"),

    # refinement
    ("refinement", 125, " for manually defining element size he: he=horizontal tank dimension/refinement, he=0.02 for ref.125"),
    ("ecH", 3., "smoothing=ech*he"),
    ("cfl", 0.75, "Target cfl"),
    # run time
    ("dt_init", 0.001, "Minimum initial time step (otherwise dt_fixed/10) in s"),
    ("rampTime",3.,"duration in which the velocity magnitude is applied"),
    ("dt_output", 0.1, "output stepping"),
    ("Tend", 50., "Simulation time in s "),
    ("fract", 1.,"fraction of duration, Duration= opts.Tend/opts.fract"),
    ("dt_fixed", 0.025, "Fixed time step in s"),
    ("omega",1. ,"used to define dragAlpha")

    ])


# --- Domain
domain = Domain.PlanarStraightLineGraphDomain()
he = opts.tank_dim[0]/opts.refinement

# --- Current steady velocity                              
current=wt.SteadyCurrent(U=opts.U,
                        mwl=opts.mwl,
                        rampTime=opts.rampTime)

# ---  Water Level DownStream

if opts.outflow_level < 0.0:
    outflow_level = -(opts.tank_dim[0] ** 2) - (opts.tank_dim[1] ** 2)
else:
    outflow_level = opts.outflow_level

# --- Geometry of the obstacle
obstacle_x_end = opts.obstacle_x_start + opts.obstacle_dim[0]
obstacle_height = opts.obstacle_dim[1]

# --- Geometry of ventilated tailwater
if opts.air_vent:
    air_vent = True
    airvent_y1 = opts.airvent_y1
    airvent_y2 = airvent_y1 + opts.airvent_dim
else:
    air_vent = False

# ---  Sanity checks
if opts.waterLine_z > opts.tank_dim[1]:
    raise ValueError("ERROR: Water (level: %s) overflows height of tank (%s)"
                     % (opts.waterLine_z, opts.tank_dim[1]))
if outflow_level > opts.tank_dim[1]:
    raise ValueError("ERROR: Water (outflow level: %s) overflows height of tank (%s)"
                     % (outflow_level, opts.tank_dim[1]))
if obstacle_x_end > opts.tank_dim[0] or obstacle_height > opts.tank_dim[1]:
    raise ValueError("ERROR: Obstacle (height: %s, width: %s, start: %s) lies "
                     " outside of tank (height: %s, width: %s)"
                     % (opts.obstacle_dim[1], opts.obstacle_dim[0], opts.obstacle_x_start,
                        opts.tank_dim[1], opts.tank_dim[0]))
if opts.waterLine_x + opts.obstacle_dim[0] > opts.tank_dim[0]:
    raise ValueError("ERROR: Water starts outside of tank at x = %s (tank: %s)"
                     % (opts.waterLine_x+opts.obstacle_dim[0], opts.tank_dim[0]))
if opts.air_vent:
    if airvent_y2 > obstacle_height:
        raise ValueError("ERROR: Air ventilation (%s) exceeds the obstacle (%s)"
                         % (airvent_y2, obstacle_height))
           

# --- Tank Set up

if air_vent:
    weir = [[[opts.obstacle_x_start, 0], [opts.obstacle_x_start, obstacle_height],
             [obstacle_x_end, obstacle_height], [obstacle_x_end, airvent_y2],
             [obstacle_x_end, airvent_y1], [obstacle_x_end, 0]]]
    vent = {'airvent': [[obstacle_x_end, airvent_y2]]}
else:
    weir = [[[opts.obstacle_x_start, 0], [opts.obstacle_x_start, obstacle_height],
             [obstacle_x_end, obstacle_height], [obstacle_x_end, 0]]]
    vent = None

tank = st.TankWithObstacles2D(domain=domain,
                              dim=opts.tank_dim,
                              obstacles=weir,
                              special_boundaries=vent)
 
# --- Genetation & Absorption Zones
he=opts.tank_dim[0]/opts.refinement
smoothing=opts.ecH*he
dragAlpha = 5.*opts.omega/opts.nu_0

tank.setSponge(x_n = opts.tank_sponge[0], x_p = opts.tank_sponge[1])


tank.setGenerationZones(x_n=True, 
 			waves=current,
 			dragAlpha=dragAlpha,
			smoothing = smoothing)

#tank.setAbsorptionZones(x_n=True, dragAlpha = dragAlpha)
tank.setAbsorptionZones(x_p=True, dragAlpha = dragAlpha)

# --- Boundary Conditions

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=opts.outflow_level,
                                                    rhoUp=opts.rho_1,
                                                    rhoDown=opts.rho_0,
                                                    g=np.array([0.,-9.805,0]),
                                                    refLevel=1.,
                                                    smoothing=smoothing
 
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=current, 
                                               vert_axis=None,
					       smoothing=3.0*he
                                               )

tank.BC['x-'].p_advective.uOfXT = lambda x, t: - opts.U[0]
tank.BC['sponge'].setNonMaterial()
    
if air_vent:
    tank.BC['airvent'].reset()
    tank.BC['airvent'].p_dirichlet.uOfXT = lambda x, t: (opts.tank_dim[1] - x[1])*opts.rho_1*abs(opts.g[1])
    tank.BC['airvent'].v_dirichlet.uOfXT = lambda x, t: 0.0
    tank.BC['airvent'].vof_dirichlet.uOfXT = lambda x, t: 1.0
    tank.BC['airvent'].u_diffusive.uOfXT = lambda x, t: 0.0
    tank.BC['airvent'].v_diffusive.uOfXT = lambda x, t: 0.0
    

#--- Initial Conditions

def wavePhi(x,t):
    return x[1] - opts.waterLine_z

def outflowPhi(x,t):
    return x[1] - outflow_level

def twpflowVelocity_u(x,t):
    waterspeed = opts.u[0]
    H = smoothedHeaviside(ecH*he,wavePhi(x,t)-ecH*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_u_D(x, t):
    waterspeed = opts.outflow_velocity
    H = smoothedHeaviside(ecH * he, outflowPhi(x, t) - ecH * he)
    u = H * windVelocity[0] + (1.0 - H) * waterspeed
    return u

def signedDistance(x):
    phi_x = x[0] - opts.waterLine_x
    phi_z = x[1] - opts.waterLine_z
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
class P_IC:
    def __init__(self):
        self.mwl=opts.water_level
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_dim[1] - opts.water_level)*opts.rho_1*opts.g[1] - (opts.water_level - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_dim[1] - opts.water_level)*opts.rho_1*opts.g[1]
class AtRest:
    def uOfXT(self, x, t):
        return 0.0
class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH*he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

                
"""
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

"""   

# --- Two Phase Flow

dt_init = opts.dt_init
nDTout = int(round(opts.Tend / opts.dt_fixed))
Duration= opts.Tend/opts.fract

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest()}
initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()
            
outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                                   dt_init=opts.dt_init,
                                                   # cfl=cfl,
                                                   dt_output=opts.dt_output,
                                                   nDTout=None,
                                                   dt_fixed=None)
            
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                            ls_model=None,
                                            nd=domain.nd,
                                            cfl=opts.cfl,
                                            outputStepping=outputStepping,
                                            structured=False,
                                            he=he,
                                            nnx=None,
                                            nny=None,
                                            nnz=None,
                                            domain=domain,
                                            initialConditions=initialConditions,
                                            boundaryConditions=None, 
                                             )
params = myTpFlowProblem.Parameters
            
myTpFlowProblem.useSuperLu=False#True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma_01
            
# index in order of
m = params.Models
m.moveMeshElastic.index=0
m.rans2p.index = 1
m.vof.index = 2
m.ncls.index = 3
m.rdls.index = 4
m.mcorr.index = 5
            
# Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
                
# assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)
