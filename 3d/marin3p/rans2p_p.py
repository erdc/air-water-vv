from proteus.default_p import *
from proteus import Context
from proteus.mprans import RANS2P

name="rans2p"
ct = Context.get()
nd=ct.nd
domain=ct.domain
genMesh=ct.genMesh
LevelModelType = RANS2P.LevelModel
if ct.opts.useOnlyVF:
    LS_model = None
else:
    LS_model = 2

coefficients = RANS2P.Coefficients(epsFact=ct.epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = ct.rho_0,
                                   nu_0 = ct.nu_0,
                                   rho_1 = ct.rho_1,
                                   nu_1 = ct.nu_1,
                                   g=ct.g,
                                   nd=nd,
                                   ME_model=ct.V_model,
                                   VF_model=ct.VOF_model,
                                   LS_model=ct.LS_model,
                                   epsFact_density=ct.epsFact_density,
                                   stokes=False,
                                   useVF=ct.useVF,
				   useRBLES=ct.useRBLES,
				   useMetrics=1.0,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=ct.weak_bc_penalty_constant,
                                   forceStrongDirichlet=ct.ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ct.ns_closure,
                                   NONCONSERVATIVE_FORM=1.0)

def getDBC_p(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return lambda  x,t: 0.0

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
    
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

if ct.bcCoords:
    assert(False)
    eps=1.0e-4
    def getAFBC_u(x,flag):
        if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
            return lambda x,t: 0.0
        if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_v(x,flag):
        if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
            return lambda x,t: 0.0
        if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0


    def getAFBC_w(x,flag):
        if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
            return lambda x,t: 0.0
        if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0


    def getDFBC_u(x,flag):
        return lambda x,t: 0.0

    def getDFBC_v(x,flag):
        return lambda x,t: 0.0

    def getDFBC_w(x,flag):
        return lambda x,t: 0.0
else:
    def getAFBC_p(x,flag):
        if ct.openTop and flag == ct.boundaryTags['top']:
            return None
        else:
            return lambda x,t: 0.0

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
        if ct.openTop and flag == ct.boundaryTags['top']:
            return None
        else:
            return lambda x,t: 0.0

    def getDFBC_v(x,flag):
        if ct.openTop and flag == ct.boundaryTags['top']:
            return None
        else:
            return lambda x,t: 0.0

    def getDFBC_w(x,flag):
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if ct.signedDistance(x) < 0:
            return -(ct.L[2] - self.waterLevel)*ct.rho_1*ct.g[2] - (self.waterLevel - x[2])*ct.rho_0*ct.g[2]
        else:
            return -(ct.L[2] - self.waterLevel)*ct.rho_1*ct.g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(ct.waterLine_z),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
