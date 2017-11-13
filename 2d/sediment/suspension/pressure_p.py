from math import *
from proteus import *
from proteus.default_p import *
from proteus.mprans import Pres
from tank import *
from proteus import Context
import tank_so

ct = Context.get()
name = "pressure"

coefficients=Pres.Coefficients(modelIndex=tank_so.P_model,
                               fluidModelIndex=tank_so.FLOW_model,
                               pressureIncrementModelIndex=tank_so.PINC_model,
                               useRotationalForm=False)
LevelModelType = Pres.LevelModel

def getDBC_p(x,flag):
    if flag == boundaryTags['y+'] and openTop:
        return lambda x,t: 0.0

def getFlux(x,flag):
    if not openTop or flag != boundaryTags['y+']:
        return lambda x,t: 0.0

class getIBC_p:
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

initialConditions = {0:getIBC_p(waterLine_z)}
dirichletConditions = {0:getDBC_p} # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0:getFlux}
