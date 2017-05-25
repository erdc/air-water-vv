from proteus.default_p import *
from proteus import Context
from proteus.mprans import NCLS

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,
                                 RD_model=3,
                                 ME_model=2,
                                 checkMass=False, 
                                 useMetrics=ct.useMetrics,
                                 epsFact=ct.ecH,
                                 sc_uref=ct.ls_sc_uref,
                                 sc_beta=ct.ls_sc_beta,
                                 movingDomain=ct.movingDomain)

def getDBC_ls(x,flag):
    if flag == ct.tank.boundaryTags['x-']:
        return ct.wavePhi
#    elif flag == boundaryTags['right']:
#        return  outflowPhi
    else:
        return None

dirichletConditions = {0: getDBC_ls}
advectiveFluxBoundaryConditions =  {0: lambda x, flag: None}
diffusiveFluxBoundaryConditions = {0: {}}

class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return ct.signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
