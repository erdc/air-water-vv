from math import *
from proteus import *
from proteus.default_p import *
from proteus.mprans import Pres
from proteus import Context
name = "pressure"

ct = Context.get()
nd = ct.nd
domain = ct.domain

coefficients=Pres.Coefficients(modelIndex=ct.PRESSURE_model,
                               fluidModelIndex=ct.V_model,
                               pressureIncrementModelIndex=ct.PINC_model,
                               useRotationalForm=False)

def getDBC_p(x,flag):
    if flag == ct.boundaryTags['top'] and ct.openTop:
        return lambda x,t: 0.0
    
def getFlux(x,flag):
    if not(flag == ct.boundaryTags['top'] and ct.openTop):
        return lambda x,t: 0.0
    
class getIBC_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
       if ct.signedDistance(x) < 0:
           return -(ct.L[2] - self.waterLevel)*ct.rho_1*ct.g[2] - (self.waterLevel - x[2])*ct.rho_0*ct.g[2]
       else:
           return -(ct.L[2] - x[2])*ct.rho_1*ct.g[2]

initialConditions = {0:getIBC_p(ct.waterLine_z)}

dirichletConditions = {0:getDBC_p } # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0:getFlux}
