from proteus import *
from proteus.default_p import *
from tank import *
from proteus.mprans import Dissipation
from proteus import Context

ct = Context.get()

LevelModelType = Dissipation.LevelModel

dissipation_model_flag = 1
if ct.useRANS >= 2:
    dissipation_model_flag=2
coefficients = Dissipation.Coefficients(V_model=ct.V_model,
                                        ME_model=ct.EPS_model,
                                        LS_model=ct.LS_model,
                                        RD_model=ct.RD_model,
                                        kappa_model=ct.K_model,
                                        SED_model=ct.SED_model,
                                        dissipation_model_flag=dissipation_model_flag+int(ct.movingDomain),#1 -- K-epsilon, 2 -- K-omega
                                        useMetrics=useMetrics,
                                        rho_0=rho_0,
                                        nu_0=nu_0,
                                        rho_1=rho_1,
                                        nu_1=nu_1,
                                        g=g,
                                        c_mu=ct.opts.Cmu,sigma_e=ct.opts.sigma_e,
                                        nd = ct.nd,
                                        sc_uref=dissipation_sc_uref,
                                        sc_beta=dissipation_sc_beta,
                                        closure = ct.sedClosure )

kInflow=ct.dissipationInflow

dissipationInflow=ct.dissipationInflow


dirichletConditions = {0: lambda x, flag: domain.bc[flag].dissipation_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].dissipation_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: domain.bc[flag].dissipation_diffusive.init_cython()},
                                  }

 
class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval

initialConditions  = {0:ConstantIC(cval=dissipationInflow)}
