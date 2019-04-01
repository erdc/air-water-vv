from math import *
from proteus import *
from proteus.default_p import *
from circle_caisson2D_oscillationIB import *
from proteus.mprans import Pres

name = "pressure"

coefficients=Pres.Coefficients(modelIndex=6,
                               fluidModelIndex=4,
                               pressureIncrementModelIndex=5,
                               useRotationalForm=False)

LevelModelType = Pres.LevelModel

def getDBC_p(x,flag):
    return None

def getFlux(x,flag):
    return lambda x,t: 0.0

class getIBC_p:
    def uOfXT(self,x,t):
        return -(L[1] - x[1])*rho_1*g[1]

initialConditions = {0:getIBC_p()}

dirichletConditions = {0:getDBC_p } # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0:getFlux}
