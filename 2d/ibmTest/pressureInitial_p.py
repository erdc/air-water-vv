from math import *
from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import PresInit

#domain = ctx.domain
#nd = ctx.nd
name = "pressureInitial"

coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=PINIT_model,
                                   fluidModelIndex=V_model,
                                   pressureModelIndex=PRESSURE_model)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].pIni_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pIni_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0:{0: lambda x, flag: domain.bc[flag].pIni_diffusive.init_cython()}}


class getIBC_pInit:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_pInit()}


