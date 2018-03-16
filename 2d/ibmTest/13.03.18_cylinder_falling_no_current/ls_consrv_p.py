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

coefficients = MCorr3P.Coefficients(LS_model= 1,
                                  V_model=4,
                                  ME_model= 3,
                                  VOF_model=0,
                                  applyCorrection=ct.applyCorrection,
                                  nd=nd,
                                  checkMass=True,
                                  useMetrics=ct.useMetrics,
                                  epsFactHeaviside=3,
                                  epsFactDirac=3,
                                  epsFactDiffusion=1)

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self, X):
        return 0.0
    def uOfXT(self, X, t):
        return 0.0

initialConditions = {0: zero_phi()}
