from math import *
from proteus import *
from proteus.default_p import *
from proteus import Context
from tank import *

ct = Context.get()
tank_dim = ct.L

eps=1.0e-8
def onLeft(x):
    return x[0] < eps and x[1] > eps and x[1] < tank_dim[1] - eps
def onRight(x):
    return x[0] > tank_dim[0] - eps and x[1] > eps and x[1] < tank_dim[1] - eps
def onBottom(x):
    return x[1] < eps
def onTop(x):
    return x[1] > tank_dim[1] - eps

def zero(x, t):
    return 0.0
def one(x, t):
    return 1.0
def smol(x, t):
    return 1.0e-10
def smoller(x, t):
    return 1.0e-20

def noneTop(x, t):
	if onLeft(x) or onRight(x) or onBottom(x):
		return zero(x, t)
	
def zeroTop(x, t):
	if onTop(x):
		return zero(x, t)
#	else:
#		return None
	
def oneTop(x, t):
	if onTop(x):
		return one(x, t)
#	else:
#		return None
	
isPer = False

# Dissipation
diss_parallelPeriodic = isPer
diss_periodic = None
diss_dirichlet = {0: lambda x, flag: smol}
diss_advective = {0: lambda x, flag: None}
diss_diffusive = {0: {},
                  1: {1: lambda x, flag: zero}}

# Kappa
kapp_parallelPeriodic = isPer
kapp_periodic = None
kapp_dirichlet = {0: lambda x, flag: smoller}
kapp_advective = {0: lambda x, flag: None}
kapp_diffusive = {0: {},
                  1: {1: lambda x, flag: zero}}

# Pressure
pres_parallelPeriodic = isPer
pres_periodic = None
pres_dirichlet = {0: lambda x, flag: zero} # pressure bc are explicitly set
pres_advective = {0: lambda x, flag: noneTop}

# Pressure Init
pInt_parallelPeriodic = isPer
pInt_periodic = None
pInt_dirichlet = {0: lambda x, flag: zeroTop}
pInt_advective = {0: lambda x, flag: noneTop}
pInt_diffusive = {0:{0: lambda x, flag: None}}

# Pressure Increment 
pInc_parallelPeriodic = isPer
pInc_periodic = None
pInc_dirichlet = {0: lambda x, flag: zeroTop}
pInc_advective = {0: lambda x, flag: noneTop}
pInc_diffusive = {0:{0: lambda x, flag: noneTop}}

# 3P Navier Stokes Sed
# 0: u_s -- 1: v_s
ns3P_parallelPeriodic = isPer
ns3P_periodic = None
ns3P_dirichlet = {0: lambda x, flag: zero,
                  1: lambda x, flag: zero}
ns3P_advective = {0: lambda x, flag: None,
                  1: lambda x, flag: None}
ns3P_diffusive = {0: {0: lambda x, flag: zeroTop},
                  1: {1: lambda x, flag: zeroTop}}
if nd == 3:
    ns3P_dirichlet[2] = lambda x, flag: zero
    ns3P_advective[2] = lambda x, flag: None
    ns3P_diffusive[2] = {2: lambda x, flag: zeroTop}

# 2P Navier Stokes
# 0: u_f -- 1: v_f
ns2P_parallelPeriodic = isPer

ns2P_periodic = None
ns2P_dirichlet = {0: lambda x, flag: zero,
                  1: lambda x, flag: zero}
ns2P_advective = {0: lambda x, flag: None,
                  1: lambda x, flag: None}
ns2P_diffusive = {0: {0: lambda x, flag: zeroTop},
                  1: {1: lambda x, flag: zeroTop}}
if nd == 3:
    ns2P_dirichlet[2] = lambda x, flag: zero
    ns2P_advective[2] = lambda x, flag: None
    ns2P_diffusive[2] = {2: lambda x, flag: zeroTop}

# Volume of Fluid
vof_parallelPeriodic = isPer
vof_periodic = None
vof_dirichlet = {0: lambda x, flag: oneTop}
vof_advective = {0: lambda x, flag: noneTop}
vof_diffusive = {0: {}}

# Volume of Solids
vos_parallelPeriodic = isPer
vos_periodic = None
vos_dirichlet = {0: lambda x, flag: zeroTop}
vos_advective = {0: lambda x, flag: noneTop}
vos_diffusive = {0: {}}