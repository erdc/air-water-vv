from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
import numpy as np
from proteus.mprans import BodyDynamics as bd
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import Gauges as ga



# Options (Geometry, Physical Properties, Waves, Numerical Options, Obstacle Dimentions)

opts=Context.Options([
    
    
    # Geometry
    ("Ly",0.7,"Vertical Dimention of the tank"),
    ("Lback",15.0,"Horizontal Dimention of overtopping collection tank"),

    # Physical Properties
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0., -9.805, 0.]), "gravity"),

    

    # Waves
    ("Tstart", 0, "Start time"),
    ("fract", 1, "fraction of duration"),
    ("Ntotalwaves",1000,"totalnumber of waves"),
    ("x0", np.array([0.,0.,0.]), "Position vector for the tinme series"),
    ("Tp", 3.5, "Peak wave period"),
    ("Hs", 0.096, "Significant wave height"),
    ("mwl", 0.4, "Mean Water level"),
    ("depth", 0.4 , "Water depth"),
    ("waveDir", np.array([1.,0.,0.]),"Direction of the waves"),
    ("N", 2000, "Number of frequency components"),
    ("bandFactor", 2.0 ,"Spectal Band Factor"),
    ("spectName", "JONSWAP","Name of Spectral Distribution"),
    ("spectral_params",{"gamma": 3.3, "TMA":False,"depth": 0.4} ,"Spectral Distribution Characteristics"),
    ("seed", 420,"Seed for random phases"),
    ("Lgen",np.array([0.,0.,0.,]) , "Length of the generation zone"),
    ("Nwaves", 15, "Number of waves per window"),
    ("Nfreq",32 , "Number of fourier components per window"),

   # Numerical Options
    ("refinement_level", 200.,"he=walength/refinement_level"),
    ("cfl", 0.5,"Target cfl"),
    ("ecH", 1.5,"Smoothing Coefficient"),
    ("Np", 10 ," Output points per period Tp/Np" ),
    ("dt_init", 0.001 , "initial time step" ),
    
    
    # Obstacle Dimensions 
    ("structure_slope", 4, "1/slope"),
    ("structureCrestLevel", 0.5, "elevation of structure crest. Equal to Water depth + Rc (crest freeboard)")
    ])



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
                         Lgen=opts.Lgen,
                         Nwaves=opts.Nwaves,
                         Nfreq=opts.Nfreq,
                          checkAcc=True,
                          fast=True)


# Script on wave length
wave_length=wave.wavelength     
           

# Domain
domain = Domain.PlanarStraightLineGraphDomain()
tank_dim=[2*wave_length+opts.structureCrestLevel*opts.structure_slope+opts.Lback,opts.Ly]


obstacle = [
           [ [2*wave_length,0.],
            [2*wave_length+opts.structureCrestLevel*opts.structure_slope,opts.structureCrestLevel],
            [2*wave_length+opts.structureCrestLevel*opts.structure_slope,0.]
            ]
            ]
                     

tank = st.TankWithObstacles2D(domain = domain, dim = tank_dim, obstacles = obstacle, hole = True)

# Mesh Refinement

he=wave_length/opts.refinement_level
ecH=opts.ecH
smoothing=ecH*he

# Boundary Conditions
boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'sponge' : 5,
                'porousLayer' : 6,
                'moving_porousLayer' : 7}

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'porousLayer': None,
                        'moving_porousLayer': None}


# Tank
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.BC['sponge'].setNonMaterial()


########################################################################################################################################################################################################################################################################################################################################################
# -----  ABSORPTION ZONE BEHIND PADDLE  ----- #
########################################################################################################################################################################################################################################################################################################################################################

tank.setSponge(x_n=wave_length)
omega = 2.*np.pi/opts.Tp

dragAlpha = 5.*omega/opts.nu_0

tank.setGenerationZones(x_n=True,
                        waves=wave,
                        dragAlpha=dragAlpha,
                        smoothing=smoothing)

# Initial Conditions

def signedDistance(x):
    waterLine_x = 2*wave_length+opts.structure_slope*opts.mwl
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
        self.waterLevel=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.Ly - self.waterLevel)*opts.rho_1*opts.g[1] - (self.waterLevel - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.Ly - self.waterLevel)*opts.rho_1*opts.g[1]
class AtRest:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest()}
class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()



# Numerics

Duration= Tend/opts.fract
dt_output = opts.Tp/opts.Np

outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                       dt_init=opts.dt_init,
                                       # cfl=cfl,
                                       dt_output=dt_output,
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
                                             boundaryConditions=None, # set with SpatialTools,
)

params = myTpFlowProblem.Parameters
myTpFlowProblem.useSuperlu=True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma_01
# index in order of
m = params.Models
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4
m.rdls.n.maxLineSearches = 0

############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
"""
gauge_x = (2*wave_length+opts.structureCrestLevel*opts.structure_slope+opts.Lback/2)



column_gauge_location = []

#column_gauge_location.append(((gauge_x, 0., 0.),
#                             (gauge_x, tank_dim[1], 0.)))


column_gauge_location.append(range([0],2*wave_length,[wave_length/10]),
                            range([2*wave_length+opts.structure_slope*opts.structureCrestLevel],[2*wave_length+opts.structure_slope*opts.structureCrestLevel+15],0.1))


tank.attachLineIntegralGauges('vof',
                              gauges=([['vof'],column_gauge_location]),
                              fileName='column_gauges.csv')    

"""
#assembling domain
domain.MeshOptions.he = he
st.assembleDomain(domain)




