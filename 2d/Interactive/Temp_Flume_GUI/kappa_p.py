from proteus.default_p import *
from proteus import Context
from proteus.mprans import Kappa
#import Context
ct = Context.get()
domain = ct.domain
nd = domain.nd

LevelModelType = Kappa.LevelModel
if ct.useOnlyVF:
    RD_model = None
    LS_model = None
    dissipation_model = 3
    ME_model = 2
else:
    RD_model = 3
    LS_model = 2
    ME_model = 5
    dissipation_model = 6
#
dissipation_model_flag = 1
if ct.useRANS == 2:
    dissipation_model_flag=2
coefficients = Kappa.Coefficients(V_model=0,
                                  ME_model=ME_model,
                                  LS_model=LS_model,
                                  RD_model=RD_model,
                                  dissipation_model=dissipation_model,
                                  dissipation_model_flag=dissipation_model_flag,#1 -- K-epsilon, 2 -- K-omega
                                  useMetrics=ct.useMetrics,
                                  rho_0=ct.rho_0,nu_0=ct.nu_0,
                                  rho_1=ct.rho_1,nu_1=ct.nu_1,
                                  g=ct.g,
                                  nd=ct.domain.nd,
                                  c_mu=0.09,
                                  sigma_k=1.0,
                                  sc_uref=ct.kappa_sc_uref,
                                  sc_beta=ct.kappa_sc_beta)

dirichletConditions = {0: lambda x, flag: domain.bc[flag].k_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].k_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].k_diffusive.init_cython()}}

class ConstantIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if ct.signedDistance(x,waterLine_z) < 0:
            return ct.kInflow
        else:
            return ct.kInflowAir
   
initialConditions  = {0:ConstantIC()}
