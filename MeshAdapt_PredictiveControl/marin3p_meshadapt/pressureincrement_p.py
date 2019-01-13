from proteus.default_p import *
from proteus.mprans import PresInc
from proteus import Context
name = "pressureincrement"
ct = Context.get()
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh

coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*ct.rho_1,
                                  rho_s_min = (1.0-1.0e-8)*ct.rho_1,
                                  nd = ct.nd,
                                  modelIndex=ct.PINC_model,
                                  fluidModelIndex=ct.V_model)

LevelModelType = PresInc.LevelModel

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_phi(x,flag):
    if flag == ct.boundaryTags['top'] and ct.openTop:
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_qt(x,flag):
    if not (flag == ct.boundaryTags['top'] and ct.openTop):
        return lambda x,t: 0.0

def getDiffusiveFlux_phi(x,flag):
    if not (flag == ct.boundaryTags['top'] and ct.openTop):
        return lambda x,t: 0.0

class getIBC_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}

dirichletConditions = {0:getDBC_phi }
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
