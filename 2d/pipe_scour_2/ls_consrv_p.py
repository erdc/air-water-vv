from proteus import StepControl
from proteus import *
from proteus.default_p import *
from proteus.mprans import MCorr3P
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = MCorr3P.LevelModel

if ct.sedimentDynamics:
    VOS_model=0
    VOF_model=1
    LS_model=2
    MCORR_model=4
    V_model=6
else:
    VOS_model=None
    VOF_model=0
    LS_model=1
    MCORR_model=3
    V_model=4

coefficients = MCorr3P.Coefficients(LS_model=LS_model,
                                    V_model=V_model,
                                    ME_model=MCORR_model,
                                    VOF_model=VOF_model,
                                    VOS_model=VOS_model,
                                    applyCorrection=ct.applyCorrection,
                                    nd=nd,
                                    checkMass=False,
                                    useMetrics=ct.useMetrics,
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
