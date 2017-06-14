"""
Wavesloshing Problem
"""
import numpy as np
from math import cos
from proteus import (Domain, Context,
                     FemTools as ft,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D

# predefined options
opts=Context.Options([
    # water
    ("water_depth_fraction", 0.5, "Depth of the water relative to the tank height"),
    ("amplitude", 0.005, "Amplitude of the sloshing relative to tank height"),
    # tank
    ("tank_dim", (0.1 , 0.1), "Dimensions of the tank"),
    #gravity
    ("g",(0,-9.81,0), "Gravity vector"),
    # gauges
    ("gauge_output", True, "Produce gauge data."),
    # mesh refinement and boundary condition
    ("he", 0.001,"Characteristic element size"),
    ("cfl", 0.5 ,"Target cfl"),
    ("opentop", True,"BC"),
    # run time options
    ("T", 0.1 ,"Simulation time"),
    ("nsave", 0., "number of saves per second"),
    ("dt_init", 0.001 ,"Maximum initial time step"),
    ("dt_fixed", None, "fixed time step for proteus (scale with period)"),
    ("structured", False, "Use a (triangular) structured mesh"),
    ("useHex", False, "Use a hexahedral structured mesh"),
    ("gen_mesh", True ,"Generate new mesh"),
    ])

# ----- CONTEXT ------ #

# tank
tank_dim = opts.tank_dim

m = 1
km = m*np.pi/tank_dim[0]
water_amplitude = opts.amplitude
water_depth = h = tank_dim[1]*opts.water_depth_fraction

wavelength = tank_dim[0]*2.
k = 2*np.pi/wavelength
h = opts.water_depth_fraction*opts.tank_dim[1]
eps = k*water_amplitude

##########################################
#     Discretization Input Options       #
##########################################

#[temp] temporary location
backgroundDiffusionFactor = 0.01 #[temp]

genMesh = opts.gen_mesh
movingDomain = False
checkMass = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = 1
useHex = opts.useHex
structured = opts.structured
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega

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
    hFactor=1.0
    if useHex:
        basis=ft.C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd,2)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd-1,2)
    else:
    	 basis=ft.C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = ft.SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd-1,3)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
        basis=ft.C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd-1,4)
    else:
        basis=ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd-1,4)

##########################################
# Numerical Options and Other Parameters #
##########################################

weak_bc_penalty_constant = 100.0
nLevels = 1

# ----- PHYSICAL PROPERTIES ----- #

# Water
rho_0 = 998.2
nu_0 = 0.#1.004e-6

# Air
rho_1 = 1.205
nu_1 = 0.#1.500e-5

# Surface Tension
sigma_01 = 0.0

# Gravity
g = opts.g

# ----- TIME STEPPING & VELOCITY----- #

T = opts.T
dt_init = opts.dt_init
runCFL = opts.cfl
nDTout = int(round(T*opts.nsave))
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
dt_fixed = opts.dt_fixed

# ----- DOMAIN ----- #

if useHex:
    nnx = 4 * refinement + 1
    nny=2*refinement+1
    hex=True
    domain = Domain.RectangularDomain(tank_dim)
elif structured:
    nnx = 4 * refinement
    nny = 2 * refinement
    domain = Domain.RectangularDomain(tank_dim)
    boundaryTags = domain.boundaryTags
else:
    domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #

tank = Tank2D(domain, tank_dim)
tank.facets = np.array([[[0, 1, 2, 3]]])
tank.facetFlags = np.array([1])

# ----- GAUGES ----- #
gauge_dx = tank_dim[0]/100.
probes=np.linspace(0., tank_dim[0], tank_dim[0]/gauge_dx+1)
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

# ----- EXTRA BOUNDARY CONDITIONS ----- #

if opts.opentop is True:
    tank.BC['y+'].setAtmosphere()
else:
    tank.BC['y+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

# ----- MESH CONSTRUCTION ----- #

he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)


# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# ----- NUMERICAL PARAMETERS ----- #
if useMetrics:
    ns_shockCapturingFactor  = 0.25
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.25
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.25
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.25
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_viscosity = epsFact_curvature \
        = epsFact_vof = ecH \
        = epsFact_consrv_dirac = epsFact_density \
        = 3.0
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = False
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_viscosity = epsFact_curvature \
        = epsFact_vof = ecH \
        = epsFact_consrv_dirac = epsFact_density \
        = 1.5
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

# ----- NUMERICS: TOLERANCES ----- #
ns_nl_atol_res = max(1.0e-8,0.01*he**2)
vof_nl_atol_res = max(1.0e-8,0.01*he**2)
ls_nl_atol_res = max(1.0e-8,0.01*he**2)
rd_nl_atol_res = max(1.0e-8,0.01*he)
mcorr_nl_atol_res = max(1.0e-8,0.01*he**2)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2)

# ----- TURBULENCE MODELS ----- #
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4

##########################################
#                Solution                #
##########################################
# 3rd order solution

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
    t_ = (t+2*np.pi/(w_*np.sqrt(k*(-g[1])))*0.25)*(w_*np.sqrt(k*(-g[1])))
    w_0 = w0(h_)
    dt = eps*d_phi0_d_t(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_t(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_t(x_, y_, t_, h_, w_0)
    dx = eps*d_phi0_d_x(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_x(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_x(x_, y_, t_, h_, w_0)
    dy = eps*d_phi0_d_y(x_, y_, t_, h_, w_0)+eps**2*d_phi1_d_y(x_, y_, t_, h_, w_0)+0.5*eps**3*d_phi2_d_y(x_, y_, t_, h_, w_0)
    p = (-y_-dt*w_-0.5*dx**2-0.5*dy**2)*rho*abs(g[1])/k+p0
    return p

def eta(x, t):
    x_ = x*k
    h_ = h*k
    w_ = omega(h_, eps)
    t_ = (t+2*np.pi/(w_*np.sqrt(k*(-g[1])))*0.25)*(w_*np.sqrt(k*(-g[1])))
    eta = eps_eta(x_, t_, h_, eps)/k
    return eta

def signedDistance(x, t):
    #d=abs(x[1] - water_depth  - water_amplitude * cos(x[0]))
    d = x[1]-(water_depth+eta(x[0], t))
    return d
