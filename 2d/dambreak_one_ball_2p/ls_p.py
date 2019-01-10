from proteus.default_p import *
from proteus.mprans import NCLS
from dambreak import *

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=V_model,
                                 RD_model=RD_model,
                                 ME_model=LS_model,
                                 checkMass=False,
                                 useMetrics=useMetrics,
                                 epsFact=epsFact_consrv_heaviside,
                                 sc_uref=ls_sc_uref,
                                 sc_beta=ls_sc_beta,
                                 movingDomain=movingDomain)

dirichletConditions = {0: lambda x, flag: None}
advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0: {}}

class PerturbedSurface_phi:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions  = {0:PerturbedSurface_phi()}
