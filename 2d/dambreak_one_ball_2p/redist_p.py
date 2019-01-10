from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import RDLS
from dambreak import *


LevelModelType = RDLS.LevelModel
coefficients = RDLS.Coefficients(epsFact=epsFact_redistance,
                                 nModelId=LS_model,
                                 rdModelId=RD_model,
                                 useMetrics=useMetrics)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_phi(object):       
    def uOfXT(self,x,t):
        return signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
