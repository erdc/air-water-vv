from proteus import *
from proteus.default_p import *
from proteus.mprans import NCLS3P
from proteus import Context
name="ls"
ct = Context.get()
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh
LevelModelType = NCLS3P.LevelModel

coefficients = NCLS3P.Coefficients(V_model=ct.V_model,
                                   RD_model=ct.RD_model,
                                   ME_model=ct.LS_model,
                                   checkMass=True,
                                   useMetrics=ct.useMetrics,
                                   epsFact=ct.epsFact_consrv_heaviside,
                                   sc_uref=ct.ls_sc_uref,
                                   sc_beta=ct.ls_sc_beta)

def getDBC_ls(x,flag):
    return None

dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return ct.signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
