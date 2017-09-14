from math import *
from proteus import *
from proteus.default_p import *
from tank import *
from proteus.mprans import PresInit

#domain = ctx.domain
#nd = ctx.nd
name = "pressureInitial"

coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=PINIT_model,
                                   fluidModelIndex=V_model,
                                   pressureModelIndex=PRESSURE_model)

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_pInit(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_pInit(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

def getDiffusiveFlux_pInit(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

class getIBC_pInit:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[1] - self.waterLevel)*rho_1*g[1] - (self.waterLevel - x[1])*rho_0*g[1]
        else:
            return -(L[1] - self.waterLevel)*rho_1*g[1]

initialConditions = {0:PerturbedSurface_p(waterLine_z)} #getIBC_pInit()}

dirichletConditions = {0:getDBC_pInit }
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_pInit}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_pInit}}
