from proteus import *
from proteus.default_p import *
from proteus.mprans import MCorr3P
from proteus import Context
name="ls_consrv"
ct = Context.get()
nd = ct.nd
domain = ct.domain
genMesh = ct.genMesh
LevelModelType = MCorr3P.LevelModel

coefficients = MCorr3P.Coefficients(LS_model=ct.LS_model,
                                    V_model=ct.V_model,
                                    ME_model=ct.MCORR_model,
                                    VOF_model=ct.VOF_model,
                                    applyCorrection=ct.applyCorrection,
                                    nd=nd,
                                    checkMass=True,
                                    useMetrics=1.0,
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
