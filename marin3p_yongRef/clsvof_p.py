from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from proteus.mprans import CLSVOF

ct = Context.get()
name="vof"
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh

LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(V_model=ct.V_model,
                                   ME_model=ct.CLSVOF_model,
                                   useMetrics=1.0,
                                   timeOrder=2,
                                   epsFactHeaviside=1.5,
                                   epsFactDirac=1.5,
                                   lambdaFact=10.0,
                                   outputQuantDOFs=True,
                                   computeMetrics=0)
coefficients.variableNames=['phi']
name="clsvof"
#
#####################
# INITIAL CONDITION #
#####################
class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return ct.signedDistance(x)    
initialConditions  = {0:PerturbedSurface_phi()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC_vof(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return lambda x,t: 1.0

def getAFBC_vof(x,flag):
    if ct.openTop and flag == ct.boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_vof}    
advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}
