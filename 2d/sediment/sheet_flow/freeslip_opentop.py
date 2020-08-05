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
def small(x, t):
    return 1.0e-10
def smaller(x, t):
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

def getDBC_p(x,flag):
    if onTop(x):
        return lambda x,t: 0.0
    
def getAF_p(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
        
pres_parallelPeriodic = isPer
pres_periodic = None
pres_dirichlet = {0:getDBC_p} # pressure bc are explicitly set
pres_advective = {0:getAF_p}

##########################################################################
# Pressure Init
##########################################################################

def getDBC_pInit(x,flag):
    if onTop(x):
        return lambda x,t: 0.0
    
def getAF_pInit(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

def getDF_pInit(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
pInt_parallelPeriodic = isPer
pInt_periodic = None
pInt_dirichlet = {0:getDBC_pInit}
pInt_advective = {0:getAF_pInit}
pInt_diffusive = {0:{0:getDF_pInit}}

##########################################################################
# Pressure Increment
##########################################################################

def getDBC_pInc(x,flag):
    if onTop(x):
        return lambda x,t: 0.0
    
def getAF_pInc(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

def getDF_pInc(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

pInc_parallelPeriodic = isPer
pInc_periodic = None
pInc_dirichlet = {0:getDBC_pInc}
pInc_advective = {0:getAF_pInc}
pInc_diffusive = {0:{0:getDF_pInc}}

##########################################################################
# 3P Navier Stokes Sed
# 0: u_s -- 1: v_s
##########################################################################

def getDBC_us(x,flag):
     if onTop(x):
         return lambda x,t: 0.0

def getDBC_vs(x,flag):
     if onTop(x):
         return lambda x,t: 0.0

def getAF_us(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
def getAF_vs(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

def getDF_us(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
def getDF_vs(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
ns3P_parallelPeriodic = isPer
ns3P_periodic = None
ns3P_dirichlet = {0:getDBC_us,
                  1:getDBC_vs}
ns3P_advective = {0:getAF_us,
                  1:getAF_vs}
ns3P_diffusive = {0: {0:getDF_us},
                  1: {1:getAF_us}}
if nd == 3:
    ns3P_dirichlet[2] = lambda x, flag: zero
    ns3P_advective[2] = lambda x, flag: None
    ns3P_diffusive[2] = {2: lambda x, flag: zeroTop}

##########################################################################
# 2P Navier Stokes
# 0: u_f -- 1: v_f
##########################################################################

ns2P_parallelPeriodic = isPer

def getDBC_uf(x,flag):
     if onTop(x):
         return lambda x,t: 0.0

def getDBC_vf(x,flag):
     if onTop(x):
         return lambda x,t: 0.0

def getAF_uf(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
def getAF_vf(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

def getDF_uf(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0
    
def getDF_vf(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0


ns2P_periodic = None
ns2P_dirichlet = {0:getDBC_uf,
                  1:getDBC_vf}
ns2P_advective = {0:getAF_uf,
                  1:getAF_vf}
ns2P_diffusive = {0: {0:getDF_uf},
                  1: {1:getAF_uf}}
if nd == 3:
    ns2P_dirichlet[2] = lambda x, flag: zero
    ns2P_advective[2] = lambda x, flag: None
    ns2P_diffusive[2] = {2: lambda x, flag: zeroTop}

##########################################################################
# Volume of Fluid
##########################################################################

vof_parallelPeriodic = isPer

def getDBC_vof(x,flag):
    if onTop(x):
        return lambda x,t: 1.0

def getAF_vof(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

vof_periodic = None
vof_dirichlet = {0:getDBC_vof}
vof_advective = {0:getAF_vof}
vof_diffusive = {0: {}}

##########################################################################
# Volume of Solid
##########################################################################

vos_parallelPeriodic = isPer

def getDBC_vos(x,flag):
    if onTop(x):
        return lambda x,t: 0.0

def getAF_vos(x,flag):
    if onBottom(x) or onLeft(x) or onRight(x):
        return lambda x,t: 0.0

vos_periodic = None
vos_dirichlet = {0:getDBC_vos}
vos_advective = {0:getAF_vos}
vos_diffusive = {0: {}}
