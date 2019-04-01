from math import *
from proteus import *
from proteus.default_p import *
from circle_caisson2D_oscillationIB import *


#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

from proteus.mprans import PresInc
rho_s=rho_1
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_s,
                                  nd = nd,
                                  modelIndex=5,
                                  fluidModelIndex=4)

LevelModelType = PresInc.LevelModel

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_phi(x,flag):
    return None

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_qt(x,flag):
    if onBoundary(x):
        return lambda x,t: 0.0

def getDiffusiveFlux_phi(x,flag):
    if onBoundary(x):
        return lambda x,t: 0.0

class getIBC_phi:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}

dirichletConditions = {0:getDBC_phi}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
