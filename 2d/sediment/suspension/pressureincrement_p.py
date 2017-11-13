from math import *
from proteus import *
from proteus.default_p import *
from tank import *
from proteus import Context
from proteus.mprans import PresInc
import tank_so

ct = Context.get()

#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

from proteus.mprans import PresInc
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_s,
                                  nd = nd,
                                  modelIndex=tank_so.PINC_model,
                                  fluidModelIndex=tank_so.FLOW_model,
                                  sedModelIndex=tank_so.SED_model,
                                  VOF_model=tank_so.VOF_model,
                                  VOS_model=tank_so.VOS_model,
                                  fixNullSpace=fixNullSpace_PresInc, 
                                  INTEGRATE_BY_PARTS_DIV_U=ct.INTEGRATE_BY_PARTS_DIV_U_PresInc)

LevelModelType = PresInc.LevelModel

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_phi(x,flag):
    if openTop and flag == boundaryTags['y+']:
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_qt(x,flag):
    return None

def getDiffusiveFlux_phi(x,flag):
    if not openTop or flag != boundaryTags['y+']:
        return lambda x,t: 0.0

class getIBC_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}
dirichletConditions = {0:getDBC_phi}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
