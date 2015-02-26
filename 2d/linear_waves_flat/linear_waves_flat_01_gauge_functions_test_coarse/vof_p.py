from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
ct = Context.get()
from proteus.mprans import VOF
L = ct.L
T = ct.T
domain = ct.domain
nd = ct.nd

LevelModelType = VOF.LevelModel
if ct.useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2

coefficients = VOF.Coefficients(LS_model=LS_model,V_model=0,RD_model=RD_model,ME_model=1,
                                checkMass=False,useMetrics=ct.useMetrics,
                                epsFact=ct.epsFact_vof,sc_uref=ct.vof_sc_uref,sc_beta=ct.vof_sc_beta,movingDomain=movingDomain)

def getDBC_vof(x,flag):
   if flag == ct.boundaryTags['left']:
       return ct.waveVF
   elif flag == ct.boundaryTags['top']:# or x[1] >= L[1] - 1.0e-12:
       return lambda x,t: 1.0
   elif flag == ct.boundaryTags['right']:
       return  ct.outflowVF

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag == ct.boundaryTags['left']:
        return None
    elif flag == ct.boundaryTags['top']:# or x[1] >= L[1] - 1.0e-12:
        return None
    elif flag == ct.boundaryTags['right']:
        return None
    else:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(ct.epsFact_consrv_heaviside*ct.he,ct.signedDistance(x))
	    
initialConditions  = {0:PerturbedSurface_H()}
