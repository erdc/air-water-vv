from proteus import *
from proteus.default_p import *
from tank3D import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

if spongeLayer or levee or slopingSpongeLayer:
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
                                        Closure_0_model=Closure_0_model,
                                        Closure_1_model=Closure_1_model,
                                        epsFact_density=epsFact_density,
                                        stokes=False,
                                        useVF=useVF,
                                        useRBLES=useRBLES,
                                        useMetrics=useMetrics,
                                        eb_adjoint_sigma=1.0,
                                        eb_penalty_constant=weak_bc_penalty_constant,
                                        forceStrongDirichlet=ns_forceStrongDirichlet,
                                        turbulenceClosureModel=ns_closure,
                                        movingDomain=movingDomain,
                                        porosityTypes=porosityTypes,
                                        dragAlphaTypes=dragAlphaTypes,
                                        dragBetaTypes=dragBetaTypes,
                                        epsFact_solid = epsFact_solidTypes)
else:
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
                                        Closure_0_model=Closure_0_model,
                                        Closure_1_model=Closure_1_model,
                                        epsFact_density=epsFact_density,
                                        stokes=False,
                                        useVF=useVF,
                                        useRBLES=useRBLES,
                                        useMetrics=useMetrics,
                                        eb_adjoint_sigma=1.0,
                                        eb_penalty_constant=weak_bc_penalty_constant,
                                        forceStrongDirichlet=ns_forceStrongDirichlet,
                                        turbulenceClosureModel=ns_closure,
                                        movingDomain=movingDomain)


def getDBC_p(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['right']:
        return outflowPressure
    elif flag == boundaryTags['back']:
        return outflowPressure
    
def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
          return twpflowVelocity_u
    elif flag == boundaryTags['right']:
        return u_current
    elif flag == boundaryTags['back']:
        return u_current
    elif flag == boundaryTags['front']:
        return u_current

def getDBC_v(x,flag):
    if flag == boundaryTags['left']: 
          return twpflowVelocity_v
    elif flag == boundaryTags['right']:
        return v_current
    elif flag == boundaryTags['back']:
        return v_current
    elif flag == boundaryTags['front']:
        return v_current


def getDBC_w(x,flag):
    if flag == boundaryTags['left']: 
          return twpflowVelocity_w
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['back']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['front']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: -twpflowVelocity_u(x,t)
    elif flag == boundaryTags['front']:
        return lambda x,t: -currentVelocity
    elif flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    
def getAFBC_u(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    
def getAFBC_v(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    
def getDFBC_u(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0  
    elif flag == boundaryTags['empty']:
        return lambda x,t: 0.0
    else:
        return None
        
def getDFBC_v(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0 
    elif flag == boundaryTags['empty']:
        return lambda x,t: 0.0
    else:
        return None


def getDFBC_w(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['empty']:
        return lambda x,t: 0.0
    else:
        return None
    

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[2] - self.waterLevel)*rho_1*g[2] - (self.waterLevel - x[2])*rho_0*g[2]
        else:
            return -(L[2] - self.waterLevel)*rho_1*g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class WaterVelocity_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):       
        return  u_current(x,0)

class WaterVelocity_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):        
        return  v_current(x,0)

initialConditions = {0:PerturbedSurface_p(waterLine_z),
                     1:WaterVelocity_u(),
                     2:WaterVelocity_v(),
                     3:AtRest()}
