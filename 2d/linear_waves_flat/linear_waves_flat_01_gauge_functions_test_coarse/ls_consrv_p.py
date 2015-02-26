from proteus import *
from proteus.default_p import *
ct = Context.get()
from proteus.mprans import MCorr
L = ct.L
T = ct.T
domain = ct.domain
nd = ct.nd

LevelModelType = MCorr.LevelModel

coefficients = MCorr.Coefficients(LSModel_index=2,V_model=0,me_model=4,VOFModel_index=1,
                                  applyCorrection=ct.applyCorrection,nd=nd,checkMass=False,useMetrics=ct.useMetrics,
                                  epsFactHeaviside=ct.epsFact_consrv_heaviside,
                                  epsFactDirac=ct.epsFact_consrv_dirac,
                                  epsFactDiffusion=ct.epsFact_consrv_diffusion)

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return 0.0

initialConditions  = {0:zero_phi()}




