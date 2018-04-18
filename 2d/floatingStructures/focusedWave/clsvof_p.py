from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from proteus.mprans import CLSVOF
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(V_model=int(ct.movingDomain)+0,                                   
                                   ME_model=int(ct.movingDomain)+1,
                                   useMetrics=ct.useMetrics,
                                   timeOrder=2,
                                   epsFactHeaviside=ct.epsHeaviside_clsvof,
                                   epsFactDirac=ct.epsHeaviside_clsvof,
                                   lambdaFact=ct.lambdaFact_clsvof,
                                   outputQuantDOFs=True,
                                   computeMetrics=0)
coefficients.variableNames=['phi']
name="clsvof"

#####################
# INITIAL CONDITION #
#####################
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - ct.waterLevel
initialConditions = {0: PHI_IC()}

#######################
# BOUNDARY CONDITIONS #
#######################
dirichletConditions = {0: lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vof_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0: {}}



