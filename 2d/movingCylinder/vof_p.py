from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans import VOF
from proteus import Context
ct = Context.get()
genMesh = ct.genMesh
movingDomain = ct.movingDomain
L = ct.L
T = ct.T
nd = ct.nd
domain = ct.domain

LevelModelType = VOF.LevelModel

if ct.useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2

coefficients = VOF.Coefficients(LS_model=int(ct.movingDomain)+LS_model,
                                V_model=int(ct.movingDomain)+0,
                                RD_model=int(ct.movingDomain)+RD_model,
                                ME_model=int(ct.movingDomain)+1,
                                checkMass=True,
                                useMetrics=ct.useMetrics,
                                epsFact=ct.epsFact_vof,
                                sc_uref=ct.vof_sc_uref,
                                sc_beta=ct.vof_sc_beta,
                                movingDomain=ct.movingDomain)

def getDBC_vof(x,flag):
    if flag == ct.boundaryTags['top']:
        return lambda x,t: 1.0

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag != ct.boundaryTags['top']:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class VF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(ct.epsFact_consrv_heaviside*ct.he,x[1]-ct.waterLevel)

initialConditions  = {0:VF_IC()}
