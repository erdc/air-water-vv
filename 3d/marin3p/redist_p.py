from proteus import *
from proteus.default_p import *
from math import *
from proteus.mprans import RDLS3P
from proteus import Context
name="redist"
ct = Context.get()
nd = ct.nd
genMesh = ct.genMesh
"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS3P.LevelModel

coefficients = RDLS3P.Coefficients(applyRedistancing=True,
                                   epsFact=ct.epsFact_redistance,
                                   nModelId=ct.LS_model,
                                   rdModelId=ct.RD_model,
		                   useMetrics=ct.useMetrics)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS3P.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return ct.signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
