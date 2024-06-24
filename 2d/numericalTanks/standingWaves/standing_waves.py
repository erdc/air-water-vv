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
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

opts = Context.Options([
    # Geometry
    ("mwl", 1., "Water level "),
    ("tank_dim", (5., 1.5,), "Dimensions of the operational domain of the tank in m (l x h)"),
    ("tank_sponge", (5., 0.), "Length of generation/absorption zone in m (left, right)"),
    ("free_slip", True, "True/False slip condition in walls"),

    #Physical Properties 
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0, -9.805, 0]), "Gravity"),

    # Waves
    #("waveType","Fenton" ,"Linear/Fenton"),
    ("T", 1.94, "Period of the waves in s"),
    ("wave_height", 0.2, "Height of the waves in m"),
    ("wave_length",5.,"wave length"),
    ("depth", 1., "Wave depth in m"),
    ("waveDir", np.array([1., 0., 0.]), "Direction of the waves (from left boundary)"),
    ("fract", 1, "total simulation time/ chosen simulation time"),

    
    # gauges
    ("point_gauge_output", True, "Generate point gauge output"),
    ("column_gauge_output", True, "Generate column gauge output"),
    ("gauge_dx", 0.25, "Horizontal spacing of point gauges/column gauges in m"),

    # Numerical Options
    ("refinement_level", 200., "he=wavelength/refinement_level"),
    ("cfl", 0.33, "Target cfl"),
    ("Tend", 100, "Simulation time in s"),
    ("dt_init", 0.001, "Initial time step in s"),
    ("gen_mesh", True, "Generate new mesh"),
    ("useHex", False, "Use (hexahedral) structured mesh"),
    ("structured", False, "Use (triangular/tetrahedral) structured mesh"),
    ("nperiod", 10., "Number of output points per period"),
    ("ecH",3,"Smoothing Coefficient"),
    #("Nf",8,"Fenton Fourier COmponents"), 
    ("Np", 15 ," Output points per period Tp/Np" )])

# Domain
tank_dim=opts.tank_dim
domain = Domain.PlanarStraightLineGraphDomain()

#Generation/ Absorption
Lgen = np.array([opts.tank_sponge[0], 0., 0.])

# Wave Input
wave = wt.MonochromaticWaves(period=opts.T,
                            waveHeight=opts.wave_height,
                            mwl=opts.mwl,
                            depth=opts.depth,
                            g=opts.g,
                            waveDir=opts.waveDir 
			    )


tank = st.Tank2D(domain, tank_dim)

#  Sponge

tank_sponge = opts.tank_sponge
tank.setSponge(x_n=opts.tank_sponge[0], x_p=None)

# Mesh Refinement
he=opts.wave_length/opts.refinement_level
ecH=opts.ecH
smoothing=ecH*he


# Boundary Conditions
# Tank
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
tank.BC['sponge'].setNonMaterial()

# ABSORPTION ZONE BEHIND PADDLE  
dragAlpha = 5*(2*np.pi/opts.T)/1e-6

tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing = smoothing)

# Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)

waterLine_x=10000
waterLine_z = opts.mwl
def signedDistance(x):
    phi_x = x[0]- waterLine_x
    phi_y = x[1] - opts.mwl
    if phi_x < 0.0:
        if phi_y < 0.0:
            return max(phi_x, phi_y)
        else:
            return phi_y
    else:
        if phi_y < 0.0:
            return phi_x
        else:
            return np.sqrt(phi_x ** 2 + phi_y ** 2)


class P_IC:
    def __init__(self):
        self.mwl=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1] - (opts.mwl - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1]

class AtRest:
    def uOfXT(self, x, t):
        return 0.0

class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

# Numerics

Duration= opts.Tend/opts.fract
dt_output = opts.T/opts.Np

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()
myTpFlowProblem.domain=domain

# --- Timestepping for output
myTpFlowProblem.outputStepping.final_time = Duration
myTpFlowProblem.outputStepping.dt_output = dt_output
myTpFlowProblem.outputStepping.dt_init = opts.dt_init

myTpFlowProblem.SystemNumerics.cfl=opts.cfl
myTpFlowProblem.SystemNumerics.useSuperlu=False

myTpFlowProblem.SystemPhysics.setDefaults()

myTpFlowProblem.SystemPhysics.useDefaultModels(flowModel=0,interfaceModel=0)

myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['p']=P_IC()
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['u']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['v']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['vof'].p.initialConditions['vof'] = VOF_IC()
myTpFlowProblem.SystemPhysics.modelDict['ncls'].p.initialConditions['phi'] = LS_IC()
myTpFlowProblem.SystemPhysics.modelDict['rdls'].p.initialConditions['phid'] = LS_IC()
myTpFlowProblem.SystemPhysics.modelDict['mcorr'].p.initialConditions['phiCorr'] = AtRest()

params = myTpFlowProblem.SystemPhysics

params['rho_0'] = opts.rho_0  # water
params['rho_1'] = opts.rho_1  # air
params['nu_0'] = opts.nu_0  # water
params['nu_1'] = opts.nu_1  # air
params['surf_tension_coeff'] = opts.sigma_01

m = params.modelDict

m['flow'].auxiliaryVariables += domain.auxiliaryVariables['twp']

"""

 GAUGES 

column_gauge_locations = []
point_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:
    gauge_y = opts.water_level - 0.5*opts.depth
    number_of_gauges = tank_dim[0] / opts.gauge_dx + 1
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
"""
def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd - 1] - waterLevel
    phi = x[nd - 1] - waterLevel
    return p_L - g[nd - 1] * (rho_0 * (phi_L - phi) + (rho_1 - rho_0) * (
    smoothedHeaviside_integral(ecH * domain.MeshOptions.he, phi_L)
    - smoothedHeaviside_integral(ecH * domain.MeshOptions.he, phi)))



