from proteus.default_p import *
from proteus.mprans import NCLS
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
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
if ct.opts.lbc:
    bt = ct.domain.boundaryTags
    def phi_DBC(x,flag):
        if flag == bt['tank2D1_x-']:
            return ct.inlet_phi_dirichlet

    dirichletConditions = {0: phi_DBC}
else:
    def phi_DBC(x,flag):
        return None
    dirichletConditions = {0: phi_DBC}
            

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {}}

class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - ct.waterLevel

initialConditions = {0: PHI_IC()}
