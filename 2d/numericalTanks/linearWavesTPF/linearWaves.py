import numpy as np
from proteus import (Domain, Context,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow


# general parameters

rho_0 = 998.2
nu_0 = 1.004e-6
rho_1 = 1.205
nu_1 = 1.5e-5
sigma_01 = 0.
g = np.array([0., -9.81, 0.])
he = 0.02
cfl = 0.5


# wave

water_level = mwl = depth = 0.5
direction = np.array([1., 0., 0.])
height = 0.1
period = 1.2
wave = wt.MonochromaticWaves(period=period,
                             waveHeight=height,
                             mwl=mwl,
                             depth=depth,
                             g=g,
                             waveDir=direction,
                             waveType='Linear')
wavelength = wave.wavelength  # retrieve actual wavelength


# domain

domain = Domain.PlanarStraightLineGraphDomain()

tank_dim = [2.*wavelength, 2.*water_level]
tank = st.Tank2D(domain, dim=tank_dim)
tank.setSponge(x_n=wavelength, x_p=2*wavelength)

domain.MeshOptions.he = he
st.assembleDomain(domain)


# Boundary Conditions

smoothing = he*1.5
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.BC['x+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['y+'].setAtmosphere()
tank.BC['sponge'].setNonMaterial()

omega = 2.*np.pi/period
dragAlpha = 5.*omega/1e-6

tank.setGenerationZones(x_n=True,
                        waves=wave,
                        dragAlpha=dragAlpha,
                        smoothing=smoothing)
tank.setAbsorptionZones(x_p=True,
                        dragAlpha=dragAlpha)

# Initial Conditions

nd = domain.nd
g_ind = nd-1

ecH = 1.5

def signedDistance(x):
    phi_z = x[g_ind]-water_level
    return phi_z

from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[g_ind] - water_level
    phi = x[g_ind] - water_level
    return p_L -g[g_ind]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(ecH*he,phi_L)
                                                        -smoothedHeaviside_integral(ecH*he,phi)))
class P_IC:
    def uOfXT(self, x, t):
        return twpflowPressure_init(x, t)

class AtRest:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest()}
class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()



# numerics

T = 10.
dt_init = 0.001
dt_output = 0.1
dt_fixed = None
cfl = 0.5
outputStepping = TpFlow.OutputStepping(final_time=T,
                                       dt_init=dt_init,
                                       # cfl=cfl,
                                       dt_output=dt_output,
                                       nDTout=None,
                                       dt_fixed=dt_fixed)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=cfl,
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

params.physical.densityA = rho_0  # water
params.physical.densityB = rho_1  # air
params.physical.kinematicViscosityA = nu_0  # water
params.physical.kinematicViscosityB = nu_1  # air
params.physical.surf_tension_coeff = sigma_01

# index in order of
m = params.Models
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4

