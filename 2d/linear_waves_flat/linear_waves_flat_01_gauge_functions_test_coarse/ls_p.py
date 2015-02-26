from proteus import *
from proteus.default_p import *
ct = Context.get()
from proteus.mprans import NCLS
L = ct.L
T = ct.T
domain = ct.domain
nd = ct.nd

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=2,
                                 checkMass=False, useMetrics=ct.useMetrics,
                                 epsFact=ct.epsFact_consrv_heaviside,sc_uref=ct.ls_sc_uref,sc_beta=ct.ls_sc_beta,movingDomain=movingDomain)
 
def getDBC_ls(x,flag):
    if flag == ct.boundaryTags['left']:
        return ct.wavePhi
    else:
        return None

dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return ct.signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
