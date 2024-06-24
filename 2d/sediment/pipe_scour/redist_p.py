from proteus import StepControl
from proteus import *
from proteus.default_p import *
from math import *
from proteus.mprans import RDLS
from proteus import Context

"""
The redistancing equation in the sloshbox test problem.
"""

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = RDLS.LevelModel

if ct.sedimentDynamics:
    LS_model=2
    RD_model=3
else:
    LS_model=1
    RD_model=2


coefficients = RDLS.Coefficients(applyRedistancing=ct.applyRedistancing,
                                   epsFact=ct.epsFact_redistance,
                                   nModelId=LS_model,
                                   rdModelId=RD_model,
                                   useMetrics=ct.useMetrics,
                                   backgroundDiffusionFactor=ct.backgroundDiffusionFactor)

def getDBC_rd(x,flag):
    pass

dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - ct.waterLevel

initialConditions  = {0: PHI_IC()}
