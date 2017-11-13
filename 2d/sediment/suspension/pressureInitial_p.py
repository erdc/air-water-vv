from math import *
from proteus import *
from proteus.default_p import *
from proteus.mprans import PresInit
from tank import *
from proteus import Context
import tank_so

ct = Context.get()

#domain = ctx.domain
#nd = ctx.nd
name = "pressureInitial"

coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=tank_so.PINIT_model,
                                   fluidModelIndex=tank_so.FLOW_model,
                                   pressureModelIndex=tank_so.P_model)

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_pInit(x,flag):
    if flag == boundaryTags['y+']:
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_pInit(x,flag):
    if flag != boundaryTags['y+']:
        return lambda x,t: 0.0

def getDiffusiveFlux_pInit(x,flag):
    None
    #if flag != boundaryTags['top']:
    #    return lambda x,t: 0.0

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        self.vos = ct.vos_function(x)
        self.rho_a = rho_1#(1.-self.vos)*rho_1 + (self.vos)*ct.rho_s
        self.rho_w = rho_0#(1.-self.vos)*rho_0 + (self.vos)*ct.rho_s
        if signedDistance(x) < 0:
            return -(L[1] - self.waterLevel)*self.rho_a*g[1] - (self.waterLevel - x[1])*self.rho_w*g[1]
        else:
            return -(L[1] - x[1])*self.rho_a*g[1]

initialConditions = {0:PerturbedSurface_p(waterLine_z)} #getIBC_pInit()}
dirichletConditions = {0:getDBC_pInit}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_pInit}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_pInit}}
