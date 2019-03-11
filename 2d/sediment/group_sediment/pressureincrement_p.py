from math import *
from proteus import *
from proteus.default_p import *
from tank import *
import tank_so
from proteus.mprans import PresInc
from proteus import Context

ct = Context.get()

if ct.sedimentDynamics:
    VOS_model=0
    VOF_model=1
    V_model=6
    SED_model=5
    PINC_model=7
else:
    VOS_model=None
    VOF_model=0
    V_model=4
    PINC_model=5
    SED_model=None


#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

LevelModelType = PresInc.LevelModel
#from ProjectionScheme import PressureIncrement
#coefficients=PressureIncrement(rho_f_min = rho_1,
#                               rho_s_min = rho_s,
#                               nd = nd,
#                               modelIndex=PINC_model,
#                               fluidModelIndex=V_model)
from proteus.mprans import PresInc
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_s,
                                  nd = nd,
                                  modelIndex= PINC_model,
                                  fluidModelIndex= V_model,
                                  sedModelIndex= SED_model,
                                  VOF_model= VOF_model,
                                  VOS_model= VOS_model,
                                  fixNullSpace=fixNullSpace_PresInc, 
                                  INTEGRATE_BY_PARTS_DIV_U=ct.INTEGRATE_BY_PARTS_DIV_U_PresInc)


#pressure increment should be zero on any pressure dirichlet boundaries

#the advectiveFlux should be zero on any no-flow  boundaries

class getIBC_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}

dirichletConditions = {0: lambda x, flag: domain.bc[flag].pInc_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pInc_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0:{0: lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython()}}
