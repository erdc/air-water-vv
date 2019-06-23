from math import *
from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import Pres

name = "pressure"

coefficients=Pres.Coefficients(modelIndex=6,
                               fluidModelIndex=4,
                               pressureIncrementModelIndex=5,
                               useRotationalForm=False)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython()} # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0: lambda x,flag: domain.bc[flag].p_advective.init_cython()}


class getIBC_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[1] - self.waterLevel)*rho_1*g[1] - (self.waterLevel - x[1])*rho_0*g[1]
        else:
            return -(L[1] - x[1])*rho_1*g[1]


initialConditions = {0:getIBC_p(waterLine_z)}


