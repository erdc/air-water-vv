from proteus.default_p import *
from proteus.mprans import RANS3PF
from proteus import Context
name="rans3p"
ct = Context.get()
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh

LevelModelType = RANS3PF.LevelModel

coefficients = RANS3PF.Coefficients(epsFact=ct.epsFact_viscosity,
                                    sigma=0.0,
                                    rho_0 = ct.rho_0,
                                    nu_0 = ct.nu_0,
                                    rho_1 = ct.rho_1,
                                    nu_1 = ct.nu_1,
                                    g=ct.g,
                                    nd=ct.nd,
                                    ME_model=ct.V_model,
                                    PRESSURE_model=ct.PRESSURE_model,
                                    SED_model=None,
                                    VOF_model=ct.VOF_model,
                                    VOS_model=None,
                                    LS_model=ct.LS_model,
                                    Closure_0_model=None,
                                    Closure_1_model=None,
                                    epsFact_density=ct.epsFact_density,
                                    stokes=False,
                                    useVF=ct.useVF,
                                    useRBLES=ct.useRBLES,
                                    useMetrics=ct.useMetrics,
                                    eb_adjoint_sigma=1.0,
                                    eb_penalty_constant=10.0,
                                    forceStrongDirichlet=ct.ns_forceStrongDirichlet,
                                    turbulenceClosureModel=0,
                                    movingDomain=False,
                                    PSTAB=0.0)

def getDBC_u(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return lambda  x,t: 0.0

def getDBC_v(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return lambda  x,t: 0.0

def getDBC_w(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        if ct.ns_forceStrongDirichlet:
            return None
        else:
            return lambda  x,t: 0.0
    
dirichletConditions = {0:getDBC_u,
                       1:getDBC_v,
                       2:getDBC_w}

     
def getAFBC_u(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    return lambda x,t: 0.0
    
def getDFBC_v(x,flag):
    return lambda x,t: 0.0

def getDFBC_w(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v,
                                    2:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v},
                                   2:{2:getDFBC_w}}

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest(),
                     2:AtRest()}
