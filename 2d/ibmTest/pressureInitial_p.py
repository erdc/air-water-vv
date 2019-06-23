from math import *
from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import PresInit

#domain = ctx.domain
#nd = ctx.nd
name = "pressureInitial"

coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=7,
                                   fluidModelIndex=4,
                                   pressureModelIndex=6)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].pInit_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pInit_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0:{0: lambda x, flag: domain.bc[flag].pInit_diffusive.init_cython()}}


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
            return -(L[1] - x[1])*rho_1*g[1]


initialConditions = {0:getIBC_pInit()}
#initialConditions = {0:PerturbedSurface_p(waterLine_z)}


