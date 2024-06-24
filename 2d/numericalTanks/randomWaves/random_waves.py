from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
import numpy as np
from proteus.mprans import BodyDynamics as bd
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import proteus.TwoPhaseFlow.utils.Parameters as Parameters
from proteus import Gauges as ga



# Options (Geometry, Physical Properties, Waves, Numerical Options, Obstacle Dimentions)

opts=Context.Options([
    
    
    # Geometry
    ("tank_dim", (15.0, 1.5), "(x,y) dimensions of the tank domain"),
    
    # Physical Properties
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0., -9.805, 0.]), "gravity"),

    # Waves
    ("Tstart", 0, "Start time"),
    ("fract", 1, "total simulation time/ chosen simulation time"),
    ("Ntotalwaves",200,"totalnumber of waves"),
    ("x0", np.array([0.,0.,0.]), "Position vector for the tinme series"),
    ("Tp", 1.94, "Peak wave period"),
    ("Hs", 0.15, "Significant wave height"),
    ("mwl", 0.8, "Mean Water level"),
    ("depth", 0.8 , "Water depth"),
    ("waveDir", np.array([1.,0.,0.]),"Direction of the waves"),
    ("N", 2000, "Number of frequency components"),
    ("bandFactor", 2.0 ,"Spectal Band Factor"),
    ("spectName", "JONSWAP","Name of Spectral Distribution"),
    ("spectral_params",{"gamma": 3.3, "TMA":False,"depth": 0.4} ,"Spectral Distribution Characteristics"),
    ("seed", 420,"Seed for random phases"),
    ("Nwaves", 15, "Number of waves per window"),
    ("Nfreq",32 , "Number of fourier components per window"),

    # Genetation/Absorption zone
    ("tank_sponge", (5., 5.), "Length of generation/absorption zone in m (left, right)"),

    # gauges
    #("gauge_output", True, "Places Gauges in tank (10 per wavelength)"),
    ("point_gauge_output", False, "Produce point gauge output"),
    ("column_gauge_output", True, "Produce column gauge output"),
    ("gauge_dx", 0.25, "Horizontal spacing of point gauges/column gauges before structure [m]"),
   

   # Numerical Options
    ("refinement_level", 50.,"he=wavelength/refinement_level"),
    ("cfl", 0.5,"Target cfl"),
    ("ecH", 3,"Smoothing Coefficient"),
    ("Np", 15 ," Output points per period Tp/Np" ),
    ("dt_init", 0.001 , "initial time step" )
	])
    

# Domain
tank_dim=opts.tank_dim
domain = Domain.PlanarStraightLineGraphDomain()

#Generation/ Absorption
Lgen = np.array([opts.tank_sponge[0], 0., 0.])

# Wave Input

np.random.seed(opts.seed)
phi = 2*np.pi*np.random.rand(opts.N)
Tend=opts.Ntotalwaves*opts.Tp/1.1
wave = wt.RandomWavesFast(Tstart=opts.Tstart,
                         Tend=Tend,
                         x0=opts.x0,
                         Tp=opts.Tp,
                         Hs=opts.Hs,
                         mwl=opts.mwl,
                         depth=opts.depth,
                         waveDir=opts.waveDir,
                         g=opts.g,
                         N=opts.N,
                         bandFactor=opts.bandFactor,
                         spectName=opts.spectName,
                         spectral_params=opts.spectral_params,
                         phi=phi,
                         Lgen=Lgen,
                         Nwaves=opts.Nwaves,
                         Nfreq=opts.Nfreq,
                         checkAcc=True,
                         fast=True)


# Script on wave length
wave_length=wave.wavelength           

# Domain
tank = st.Tank2D(domain, tank_dim)

#shapes
tank.setSponge(x_n=opts.tank_sponge[0], x_p=opts.tank_sponge[1])

# Mesh Refinement
he=wave_length/opts.refinement_level
ecH=opts.ecH
smoothing=ecH*he
                  
# Boundary Conditions
# Tank
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
tank.BC['sponge'].setNonMaterial()

#Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)


# ABSORPTION ZONE BEHIND PADDLE  
dragAlpha = 5*(2*np.pi/opts.Tp)/1e-6

tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing = smoothing)
tank.setAbsorptionZones(x_p=True, dragAlpha = dragAlpha)

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

Duration= Tend/opts.fract
dt_output = opts.Tp/opts.Np

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

#Gauges 
column_gauge_locations = []

#if opts.point_gauge_output or opts.column_gauge_output:
gauge_y = opts.mwl - 0.5 * opts.depth
import math
number_of_gauges = math.ceil(tank_dim[0] / opts.gauge_dx + 1)
for gauge_x in np.linspace(0, tank_dim[0], number_of_gauges):column_gauge_locations.append(((gauge_x, 0., 0.),
                                       								(gauge_x, tank_dim[1], 0.)))
tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),

                                  fileName='column_gauges.csv')

"""
if opts.column_gauge_output:
tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),
                                  fileName='column_gauges.csv')

"""









