from proteus.default_p import *
from proteus.mprans import NCLS3P
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = NCLS3P.LevelModel

coefficients = NCLS3P.Coefficients(V_model=4,
                                 RD_model=2,
                                 ME_model=1,
                                 checkMass=True,
                                 useMetrics=ct.useMetrics,
                                 epsFact=3,
                                 sc_uref=ct.ls_sc_uref,
                                 sc_beta=ct.ls_sc_beta,
                                 movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: None}
advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0: {}}

class PerturbedSurface_phi:
    def uOfXT(self,x,t):
        return ct.signedDistance(x)

initialConditions  = {0:PerturbedSurface_phi()}
