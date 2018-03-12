from proteus import *
from proteus.default_p import *
from backstep import *
from proteus.mprans import RANS2P
from decimal import *
import math

LevelModelType = RANS2P.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2

dragAlphaTypes = numpy.array([0.0,
                              0.0,
                              0.0,
                              0.0])
dragBetaTypes = numpy.array([0.0,0.0,0.0,0.0])

coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_1,
                                   nu_1 = nu_1,
                                   g=g,
                                   nd=nd,
                                   VF_model=1,
                                   LS_model=LS_model,
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useVF=useVF,
				   useRBLES=useRBLES,
				   useMetrics=useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   dragAlphaTypes=dragAlphaTypes,
                                   dragBetaTypes=dragAlphaTypes,
                                   forceStrongDirichlet=ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ns_closure)

Umax = 1.5

def getDBC_p(x,flag):
    if flag==boundaryTags['back']:
        return lambda x,t: 0.0
    elif flag==boundaryTags['front1']:
        return None
    else:
        return None

def getDBC_u(x,flag):
    if flag == boundaryTags['front1']:
        return lambda x,t: 0.0
    elif flag in [boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top'], boundaryTags['front2']]:
        return lambda x,t: 0.0 
    elif flag ==boundaryTags['back']:
        return None
    else:
        return None

def getDBC_v(x,flag):
    if flag in [boundaryTags['front1']]:
        v = 8*Umax*(2*x[2]-1)*(1-x[2])
        return lambda x,t: v
    elif flag in [boundaryTags['bottom1'],boundaryTags['bottom2'], boundaryTags['top'],boundaryTags['front2']]:
        return lambda x,t: 0.0
    elif flag ==boundaryTags['back']:
        return None
    else:
        return None

def getDBC_w(x,flag):
    if flag == boundaryTags['front1']:
        return lambda x,t: 0.0
    elif flag in [ boundaryTags['bottom1'],boundaryTags['bottom2'], boundaryTags['top'],boundaryTags['front2']]:
        return lambda x,t: 0.0
    elif flag ==boundaryTags['back']:
        return None
    else:
        return None

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag==boundaryTags['front1']:
        v = 8*Umax*(2*x[2]-1)*(1-x[2])
        return lambda x,t: -v
    elif flag==boundaryTags['back']:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag==boundaryTags['back']:
        return None
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag==boundaryTags['back']:
        return None
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag == boundaryTags['back']:
        return None
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0
        
def getDFBC_w(x,flag):
    if flag == boundaryTags['back']:
        v = 2*nu_0*Umax*(1-2*x[2])
        return lambda x,t: -v
    elif flag in [boundaryTags['front1'],boundaryTags['front2'],boundaryTags['bottom1'],boundaryTags['bottom2'],boundaryTags['top']]:
        return None
    else:
        return lambda x,t: 0.0


advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Constant:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        v = 1.0
        return v

initialConditions = {0:AtRest(),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
