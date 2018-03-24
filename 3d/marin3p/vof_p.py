from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import Context
from proteus.mprans import VOF3P
ct = Context.get()
name="vof"
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh
LevelModelType = VOF3P.LevelModel
coefficients = VOF3P.Coefficients(LS_model=ct.LS_model,
                                  V_model=ct.V_model,
                                  RD_model=ct.RD_model,
                                  ME_model=ct.VOF_model,
                                  checkMass=True,
                                  useMetrics=1.0,
                                  epsFact=ct.epsFact_vof,
                                  sc_uref=ct.vof_sc_uref,
                                  sc_beta=ct.vof_sc_beta)

def getDBC_vof(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return lambda x,t: 1.0

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(ct.epsFact_consrv_heaviside*ct.he,
                                 ct.signedDistance(x))
	    
initialConditions  = {0:PerturbedSurface_H()}
