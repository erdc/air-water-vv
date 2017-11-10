from math import *
from proteus import *
from proteus.default_p import *
from proteus.mprans import Pres
from tank import *
from proteus import Context

ct = Context.get()
name = "pressure"

if ct.sedimentDynamics:
    V_model=6
    PINC_model=7
    PRESSURE_model=8

else:
    V_model=4
    PINC_model=5
    PRESSURE_model=6


coefficients=Pres.Coefficients(modelIndex=PRESSURE_model,
                               fluidModelIndex=V_model,
                               pressureIncrementModelIndex=PINC_model,
                               useRotationalForm=False)
LevelModelType = Pres.LevelModel

def getDBC_p(x,flag):
    if flag == boundaryTags['y+'] and openTop:
        return lambda x,t: 0.0

def getFlux(x,flag):
    if not(flag == boundaryTags['y+'] and openTop):
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
