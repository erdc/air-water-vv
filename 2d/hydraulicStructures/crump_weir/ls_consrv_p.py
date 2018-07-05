from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import MCorr

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = ct.domain.MeshOptions

genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = MCorr.LevelModel

coefficients = MCorr.Coefficients(LSModel_index=2,
                                  V_model=0,
                                  me_model=4,
                                  VOFModel_index=1,
                                  applyCorrection=ct.applyCorrection,
                                  nd=nd,
                                  checkMass=False,
                                  useMetrics=ct.useMetrics,
                                  epsFactHeaviside=ct.ecH,
                                  epsFactDirac=ct.epsFact_consrv_dirac,
                                  epsFactDiffusion=ct.epsFact_consrv_diffusion)

class zero_phi(object):
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return 0.0

initialConditions  = {0:zero_phi()}




