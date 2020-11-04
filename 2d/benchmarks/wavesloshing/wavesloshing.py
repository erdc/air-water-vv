"""
Wavesloshing 
"""
import numpy as np
from math import cos
from proteus import Domain
from proteus import Context 
from proteus import FemTools as ft
from proteus import MeshTools as mt
from proteus import WaveTools as wt
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import proteus.TwoPhaseFlow.utils.Parameters as Parameters



opts=Context.Options([
    # water
    ("water_depth_fraction", 0.5, "Nondimensional water level normalised to the tank height"),
    ("amplitude", 0.005, "Nondimensional amplitude normalised to tank height"),
    # tank
    ("tank_dim", [0.1 , 0.1], "Dimensions of the tank (length, height) in m"),
    #physical parameters
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 0., "Water kinematic viscosity m/sec^2"),
    ("rho_1", 1.205, "Air Densiy"),
    ("nu_1", 0., "Air kinematic viscosity m/sec^2"),
    ("sigma_01", 0., "Surface tension"),
    ("g",np.array([0,-9.81,0]), "Gravity vector"),
    # gauges
    ("gauge_output", True, "Produce gauge data. A free-surface gauge is located at the far left boundary"),
    # mesh refinement and boundary condition
    ("he", 0.001,"Characteristic element size in m"),
    ("cfl", 0.5 ,"Target cfl"),
    ("ecH",1.5 , "/ "),
    ("opentop", True,"Enabling atmosphere boundary at the top"),
    ("movingDomain", False,"True: domain (mesh) is moving"),
    ("genMesh", True ,"Generate mesh (True/False)"),
    ("constantRefinement", False, "constant refinement"),


    # run time options
    ("T", 5. ,"Simulation time in s"),
    ("nsave", 0., "number of saves per second"),
    ("dt_init", 0.001 ,"Maximum initial time step in s"),
    ("dt_fixed", None, "fixed time step for proteus (scale with period) in s"),
    ("dt_output", 0.05, "number of saves per second"),
    ("runCFL", 0.5 ,"Target CFL value"),
    ("cfl", 0.5 ,"Target CFL value"),
    ("structured", False, "Use a (triangular) structured mesh"),
    ("useHex", False, "Use a hexahedral structured mesh"),
    ("gen_mesh", True ,"Generate new mesh"),
    ("useSuperlu", False, "useSuperlu")
    ])

# Tank
water_depth = opts.tank_dim[1]*opts.water_depth_fraction
wavelength = opts.tank_dim[0]*2.
k = 2*np.pi/wavelength
h = opts.water_depth_fraction*opts.tank_dim[1]
eps = k*opts.amplitude
water_depth=h
g=np.array([0,-9.81,0])


# Domain

if opts.useHex or opts.structured:
    domain = Domain.RectangularDomain(opts.tank_dim)
    boundaryTags = domain.boundaryTags
else:
    domain = Domain.PlanarStraightLineGraphDomain()

# Tank

tank = Tank2D(domain,opts.tank_dim)
tank.facets = np.array([[[0, 1, 2, 3]]])
tank.facetFlags = np.array([1])

# Mesh Construction

he = opts.he
domain.MeshOptions.he = he
domain.MeshOptions.genMesh = opts.genMesh
st.assembleDomain(domain)

if opts.useHex:
    domain.MeshOptions.nnx = 4 * refinement + 1
    domain.MeshOptions.nny = 2*refinement+1
elif opts.structured:
    domain.MeshOptions.nnx = 4 * refinement
    domain.MeshOptions.nny = 2 * refinement


# Boundary Conditions

if opts.opentop is True:
    tank.BC['y+'].setAtmosphere()
else:
    tank.BC['y+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

# Gauges
import math
gauge_dx = opts.tank_dim[0]/100.
probes=np.linspace(0., opts.tank_dim[0], math.ceil(opts.tank_dim[0]/gauge_dx+1))
PG=[]
PG2=[]
zProbes=water_depth*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)
    PG2.append((i, water_depth, 0.),)

if opts.gauge_output:

    tank.attachPointGauges(
        'twp',
        gauges = ((('p',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_pressure.csv'
    )
    tank.attachPointGauges(
        'ls',
        gauges = ((('phi',), PG2),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_levelset.csv'
    )


class PerturbedSurface_p:
    def __init__(self,waterdepth,amplitude):
        self.waterdepth=waterdepth
        self.amplitude=amplitude
    def uOfXT(self,x,t):
        d = signedDistance(x, 0.)
        if d <= 0:
            return pressure(x[0], x[1]-self.waterdepth, t, h, eps, opts.rho_0, opts.g, k, 0.)+(opts.tank_dim[1]-(self.waterdepth+eta(x[2], 0.)))*opts.rho_1*(-opts.g[1])
            # return (ct.tank_dim[1]-(self.waterdepth+ct.eta(x)))*ct.rho_1*(-ct.opts.g[1])+((self.waterdepth+ct.eta(x))-x[1])*ct.rho_0*(-ct.g[1])
        else:
            return (opts.tank_dim[1] - x[1])*opts.rho_1*(-opts.g[1])

from proteus.ctransportCoefficients import smoothedHeaviside
class PerturbedSurface_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * opts.he, signedDistance(x, 0.))

class PerturbedSurface_phi:
    def uOfXT(self,x,t):
        return signedDistance(x, t)

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0


# Solution

def eta0(x, t):
    eta_0 = np.sin(t)*np.cos(x)
    return eta_0

def phi0(x, y, t, h, w0):
    phi_0 = w0/(np.sinh(h))*np.cos(t)*np.cos(x)*np.cosh(y+h)
    return phi_0
    
def w0(h):
    w_0 = np.sqrt(np.tanh(h))
    return w_0
    
    
def eta1(x, t, w0):
    eta_1 = 1./8.*((w0**2-w0**(-2))+(w0**(-2)-3*w0**(-6))*np.cos(2.*t))*np.cos(2.*x)
    return eta_1

def phi1(x, y, t, h, w0):
    beta0 = 0  # arbitrary
    phi_1 = beta0+1./8.(w0-w0**(-3))*t-1./16.*(3*w0+w0**(-3))*np.sin(2*t)-3./(16.*np.cosh(2.*h))*(w0-w0**(-7))*np.sin(2*t)*np.cos(2*x)*np.cosh(2*(y+h))
    return phi_1
    
w1 = 0.

def eta2(x, t, w0):
    b11 = 1./32.*(3.*w0**(-8)+6.*w0**(-4)-5.+2.*w0**4)
    b13 = 3./128.*(9.*w0**(-8)+27.*w0**(-4)-15.+w0**4+2*w0**8)
    b31 = 1./128.*(3.*w0**(-8)+18.*w0**(-4)-5.)
    b33 = 3./128.*(-9.*w0**(-12)+3.*w0**(-8)-3.*w0**(-4)+1)
    eta_2 = b11*np.sin(t)*np.cos(x)+b13*np.sin(t)*np.cos(3*x)+b31*np.sin(3*t)*np.cos(x)+b33*np.sin(3*t)*np.cos(3*x)
    return eta_2

def phi2(x, y, t, h, w0):
    beta13 = 1./(128.*np.cosh(3.*h))*(1+3*w0**4)*(3*w0**(-9)-5.*w0**(-1)+2*w0**3)
    beta31 = 1./(128.*np.cosh(h))*(9.*w0**(-9)+62.*w0**(-5)-31.*w0**(-1))
    beta33 = 1./(128.*np.cosh(3.*h))*(1+3.*w0**4)*(-9.*w0**(-13)+22.*w0**(-9)-13.*w0**(-5))
    beta2 = 0.  # arbitrary
    phi_2 = beta2+beta13*np.cos(t)*np.cos(3*x)*np.cosh(3*(y+h))+beta31*np.cos(3*t)*np.cos(x)*np.cosh(y+h)+beta33*np.cos(3*t)*np.cos(3*x)*np.cosh(3*(y+h))
    return phi_2
        
def eps_eta(x, t, h, eps):
    w_0 = w0(h)
    eta_0 = eta0(x, t)
    eta_1 = eta1(x, t, w_0)
    eta_2 = eta2(x, t, w_0)
    epseta = eps*eta_0+eps**2*eta_1+0.5*eps**3*eta_2
    return epseta

def w2(w0):
    w_2 = 1./32.*(9.*w0**(-7)-12.*w0**(-3)-3*w0-2*w0**5)
    return w_2
    
def omega(h, eps):
    w_0 = w0(h)
    w_2 = w2(w_0)
    w = w_0+0.5*eps**2*w_2
    return w
    

def d_phi0_d_t(x, y, t, h, w0):
    return -w0/np.sinh(h)*np.sin(t)*np.cos(x)*np.cosh(y+h)

def d_phi0_d_x(x, y, t, h, w0):
    return -w0/np.sinh(h)*np.cos(t)*np.sin(x)*np.cosh(y+h)

def d_phi0_d_y(x, y, t, h, w0):
    return w0/np.sinh(h)*np.cos(t)*np.cos(x)*np.sinh(y+h)

def d_phi1_d_t(x, y, t, h, w0):
    return 1./8.*(w0-w0**(-3))-2./16.*(3.*w0+w0**(-3))*np.cos(2*t)-2.*3./(16.*np.cosh(2.*h))*(w0-w0**(-7))*np.cos(2.*t)*np.cos(2.*x)*np.cosh(2.*(y+h))

def d_phi1_d_x(x, y, t, h, w0):
    return 2.*3./(16.*np.cosh(2.*h))*(w0-w0**(-7))*np.sin(2.*t)*np.sin(2.*x)*np.cosh(2.*(y+h))

def d_phi1_d_y(x, y, t, h, w0):
    return 2.*3./(16.*np.cosh(2.*h))*(w0-w0**(-7))*np.sin(2.*t)*np.cos(2.*x)*np.sinh(2.*(y+h))
    
def d_phi2_d_t(x, y, t, h, w0):
    beta13 = 1./(128.*np.cosh(3.*h))*(1+3*w0**4)*(3*w0**(-9)-5.*w0**(-1)+2*w0**3)
    beta31 = 1./(128.*np.cosh(h))*(9.*w0**(-9)+62.*w0**(-5)-31.*w0**(-1))
    beta33 = 1./(128.*np.cosh(3.*h))*(1+3.*w0**4)*(-9.*w0**(-13)+22.*w0**(-9)-13.*w0**(-5))
    return -beta13*np.sin(t)*np.cos(3.*x)*np.cosh(3.*(y+h))-3.*beta31*np.sin(3.*t)*np.cos(x)*np.cosh(y+h)-3.*beta33*np.sin(3.*t)*np.cos(3.*x)*np.cosh(3.*(y+h))

def d_phi2_d_x(x, y, t, h, w0):
    beta13 = 1./(128.*np.cosh(3.*h))*(1+3*w0**4)*(3*w0**(-9)-5.*w0**(-1)+2*w0**3)
    beta31 = 1./(128.*np.cosh(h))*(9.*w0**(-9)+62.*w0**(-5)-31.*w0**(-1))
    beta33 = 1./(128.*np.cosh(3.*h))*(1+3.*w0**4)*(-9.*w0**(-13)+22.*w0**(-9)-13.*w0**(-5))
    return -3.*beta13*np.cos(t)*np.sin(3.*x)*np.cosh(3.*(y+h))-beta31*np.cos(3.*t)*np.sin(x)*np.cosh(y+h)-3.*beta33*np.cos(3.*t)*np.sin(3.*x)*np.cosh(3.*(y+h))

def d_phi2_d_y(x, y, t, h, w0):
    beta13 = 1./(128.*np.cosh(3.*h))*(1+3*w0**4)*(3*w0**(-9)-5.*w0**(-1)+2*w0**3)
    beta31 = 1./(128.*np.cosh(h))*(9.*w0**(-9)+62.*w0**(-5)-31.*w0**(-1))
    beta33 = 1./(128.*np.cosh(3.*h))*(1+3.*w0**4)*(-9.*w0**(-13)+22.*w0**(-9)-13.*w0**(-5))
    return 3.*beta13*np.cos(t)*np.cos(3.*x)*np.sinh(3.*(y+h))+beta31*np.cos(3.*t)*np.cos(x)*np.sinh(y+h)+3.*beta33*np.cos(3.*t)*np.cos(3.*x)*np.sinh(3.*(y+h))   

def pressure(x, y, t, h, eps, rho, g, k, p0):
    x_ = x*k
    y_ = y*k
    h_ = h*k
    w_ = omega(h_, eps)
    t_ = (t+2*np.pi/(w_*np.sqrt(k*(-opts.g[1])))*0.25)*(w_*np.sqrt(k*(-opts.g[1])))
    w_0 = w0(h_)
    dt = eps*d_phi0_d_t(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_t(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_t(x_, y_, t_, h_, w_0)
    dx = eps*d_phi0_d_x(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_x(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_x(x_, y_, t_, h_, w_0)
    dy = eps*d_phi0_d_y(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_y(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_y(x_, y_, t_, h_, w_0)
    p = (-y_-dt*w_-0.5*dx**2-0.5*dy**2)*rho*abs(opts.g[1])/k+p0
    return p

def eta(x, t):
    x_ = x*k
    h_ = h*k
    w_ = omega(h_, eps)
    t_ = (t+2*np.pi/(w_*np.sqrt(k*(-opts.g[1])))*0.25)*(w_*np.sqrt(k*(-opts.g[1])))
    eta = eps_eta(x_, t_, h_, eps)/k
    return eta

def signedDistance(x, t):
    #d=abs(x[1] - water_depth  - opts.amplitude * cos(x[0]))
    d = x[1]-(water_depth+eta(x[0], t))
    return d
    
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()
myTpFlowProblem.domain=domain

# --- Timestepping for output
myTpFlowProblem.outputStepping.final_time = opts.T
myTpFlowProblem.outputStepping.dt_output = opts.dt_output
myTpFlowProblem.outputStepping.dt_init = opts.dt_init
myTpFlowProblem.outputStepping.dt_fixed = opts.dt_fixed

myTpFlowProblem.SystemNumerics.cfl=opts.cfl
myTpFlowProblem.SystemNumerics.useSuperlu=opts.useSuperlu

myTpFlowProblem.SystemPhysics.setDefaults()
myTpFlowProblem.SystemPhysics.movingDomain = opts.movingDomain

myTpFlowProblem.SystemPhysics.addModel(Parameters.ParametersModelMoveMeshElastic,'move')
myTpFlowProblem.SystemPhysics.useDefaultModels(flowModel=0,interfaceModel=0)

myTpFlowProblem.SystemPhysics.modelDict['move'].p.initialConditions['hx']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['move'].p.initialConditions['hy']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['p']=PerturbedSurface_p(water_depth, opts.amplitude)
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['u']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['v']=AtRest()
myTpFlowProblem.SystemPhysics.modelDict['vof'].p.initialConditions['vof'] = PerturbedSurface_H()
myTpFlowProblem.SystemPhysics.modelDict['ncls'].p.initialConditions['phi'] =  PerturbedSurface_phi()
myTpFlowProblem.SystemPhysics.modelDict['rdls'].p.initialConditions['phid'] =  PerturbedSurface_phi()
myTpFlowProblem.SystemPhysics.modelDict['mcorr'].p.initialConditions['phiCorr'] = AtRest()

params = myTpFlowProblem.SystemPhysics

params['rho_0'] = opts.rho_0  # water
params['rho_1'] = opts.rho_1  # air
params['nu_0'] = opts.nu_0  # water
params['nu_1'] = opts.nu_1  # air
params['surf_tension_coeff'] = opts.sigma_01

m = params.modelDict

m['flow'].auxiliaryVariables += domain.auxiliaryVariables['twp']
