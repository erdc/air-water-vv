from math import *
from proteus import *
from proteus.default_p import *
from proteus import Context
from tank import *

ct = Context.get()
tank_dim = ct.L

eps=1.0e-8
def onLeft(x):
    return x[0] < eps
def onRight(x):
    return x[0] > tank_dim[0] - eps
def onBottom(x):
    return x[1] < eps
def onTop(x):
    return x[1] > tank_dim[1] - eps

def getPDBC(x,flag):
    if onLeft(x) or onRight(x):
        return np.array([0.0,round(x[1],5),0.0])

##### Set zero on top or bottom
def setZero(x,flag):
    if onTop(x) or onBottom(x):
        return lambda x,t: 0.0
def setZeroFluxOnPeriodic(x,flag):
    if not (onTop(x) or onBottom(x)):
        return lambda x,t: 0.0
def setAllZero(x,flag):
    return lambda x,t: 0.0

##### For no slip conditions, set the following zero
# Dirichlet: u, v, us, vs
# Advective Flux: p, pInit, pInc, vof, vos
# Diffusive Flux: pInc
    
def zero(x):
    return 0.0
def one(x, t):
    return 1.0
def small(x, t):
    return 1.0e-10
def smaller(x, t):
    return 1.0e-20


isPer = False

### Not using a turbulence model anyhow

# Dissipation
diss_parallelPeriodic = isPer
diss_periodic = None
diss_dirichlet = {0: lambda x, flag: small}
diss_advective = {0: lambda x, flag: None}
diss_diffusive = {0: {},
                  1: {1: lambda x, flag: zero}}

# Kappa
kapp_parallelPeriodic = isPer
kapp_periodic = None
kapp_dirichlet = {0: lambda x, flag: smaller}
kapp_advective = {0: lambda x, flag: None}
kapp_diffusive = {0: {},
                  1: {1: lambda x, flag: zero}}

##########################################################################
# Pressure
##########################################################################

pres_parallelPeriodic = False
pres_periodic = None
pres_dirichlet = {0: lambda x, flag: None}
pres_advective = {0: setZero}
pres_diffusive = {0:{0: lambda x, flag: None}}

##########################################################################
# Pressure Init
##########################################################################

pInt_parallelPeriodic = False
pInt_periodic = None
pInt_dirichlet = {0: lambda x, flag: None}
pInt_advective = {0: setZero}
pInt_diffusive = {0:{0: lambda x, flag: None}}

##########################################################################
# Pressure Increment
##########################################################################

pInc_parallelPeriodic = False
pInc_periodic = None#{0:getPDBC}
def pIncDirichlet(x,flag):
    if onRight(x):
        return lambda x,t: 0.0
pInc_dirichlet = {0: pIncDirichlet}

pInc_advective = {0: setAllZero}
def pIncFlux(x,flag):
    if onLeft(x):
        return lambda x,t: -1.0
    elif onRight(x):
        return None
    else:
        return lambda x,t: 0.0

pInc_diffusive = {0:{0:pIncFlux}}

##########################################################################
# 3P Navier Stokes Sed
# 0: u_s -- 1: v_s
##########################################################################
    
ns3P_parallelPeriodic = False

ns3P_periodic = None; dummy = {0:getPDBC,
                  1:getPDBC}
def u_dirichlet(x,flag):
    if onLeft(x) or onRight(x):
        return lambda x,t: 1.0
    else:
        return lambda x,t: 0.0

def v_dirichlet(x,flag):
    return lambda x,t: 0.0

ns3P_dirichlet = {0: u_dirichlet,
                  1: v_dirichlet}

def mom_flux(x,flag):
    return None

def mom_dflux(x,flag):
    if onRight(x):
        return lambda x,t: 0.0

ns3P_advective = {0: mom_flux,
                  1: mom_flux}
ns3P_diffusive = {0: {0: mom_dflux},
                  1: {1: mom_dflux}}

##########################################################################
# 2P Navier Stokes
# 0: u_f -- 1: v_f
##########################################################################

ns2P_parallelPeriodic = False

ns2P_periodic = None; dummy = {0:getPDBC,
                  1:getPDBC}
def us_dirichlet(x,flag):
    if onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0

def vs_dirichlet(x,flag):
    return lambda x,t: 0.0

ns2P_dirichlet = {0: us_dirichlet,
                  1: vs_dirichlet}

ns2P_advective = {0: mom_flux,
                  1: mom_flux}
ns2P_diffusive = {0: {0: mom_dflux},
                  1: {1: mom_dflux}}

##########################################################################
# Volume of Fluid
##########################################################################

vof_parallelPeriodic = False#True

vof_periodic = None#{0:getPDBC}
vof_dirichlet = {0: lambda x, flag: None}
vof_advective = {0:setZero}
vof_diffusive = {0: {}}

##########################################################################
# Volume of Solid
##########################################################################

vos_parallelPeriodic = False#

vos_periodic = None#{0:getPDBC}
vos_dirichlet = {0: lambda x, flag :None}
vos_advective = {0:setAllZero}
vos_diffusive = {0: {}}
