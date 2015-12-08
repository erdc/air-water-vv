from proteus.default_p import *
from proteus.mprans import NCLS
from proteus import Context
ct = Context.get()

domain = ct.domain
nd = domain.nd

genMesh = ct.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=int(ct.movingDomain)+0,
                                 RD_model=int(ct.movingDomain)+3,
                                 ME_model=int(ct.movingDomain)+2,
                                 checkMass=False,
                                 useMetrics=ct.useMetrics,
                                 epsFact=ct.epsFact_consrv_heaviside,
                                 sc_uref=ct.ls_sc_uref,
                                 sc_beta=ct.ls_sc_beta,
                                 movingDomain=ct.movingDomain)


def getDBC_ls(x,flag):
    if flag == ct.boundaryTags['left']:
        return wavePhi
    else:
        return None


dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {}}



class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - ct.waterLevel

initialConditions = {0: PHI_IC()}
