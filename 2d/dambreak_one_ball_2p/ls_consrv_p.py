from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MCorr
from dambreak import *

LevelModelType = MCorr.LevelModel
coefficients = MCorr.Coefficients(LSModel_index=LS_model,
                                  V_model=V_model,
                                  me_model=MC_model,
                                  VOFModel_index=VF_model,
                                  nd=nd,
                                  checkMass=False,
                                  useMetrics=useMetrics,
                                  epsFactHeaviside=epsFact_consrv_heaviside,
                                  epsFactDirac=epsFact_consrv_dirac,
                                  epsFactDiffusion=epsFact_consrv_diffusion)

class zero_phi(object):
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return 0.0

initialConditions  = {0:zero_phi()}
