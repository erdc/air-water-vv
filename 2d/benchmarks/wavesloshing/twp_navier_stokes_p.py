from builtins import object
from proteus.default_p import *
from proteus.mprans import RANS2P
import numpy as np
from math import cos
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

LevelModelType = RANS2P.LevelModel
if ct.useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if ct.useRANS >= 1:
    Closure_0_model = 5
    Closure_1_model = 6
    if ct.useOnlyVF:
        Closure_0_model=2
        Closure_1_model = 3
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS2P.Coefficients(epsFact=ct.epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0=ct.rho_0,
                                   nu_0=ct.nu_0,
                                   rho_1=ct.rho_1,
                                   nu_1=ct.nu_1,
                                   g=ct.g,
                                   nd=nd,
                                   ME_model=int(ct.movingDomain)+0,
                                   VF_model=int(ct.movingDomain)+1,
                                   LS_model=int(ct.movingDomain)+LS_model,
                                   Closure_0_model=Closure_0_model,
                                   Closure_1_model=Closure_1_model,
                                   epsFact_density=ct.epsFact_density,
                                   stokes=False,
                                   useVF=ct.useVF,
                                   useRBLES=ct.useRBLES,
                                   useMetrics=ct.useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=ct.weak_bc_penalty_constant,
                                   forceStrongDirichlet=ct.ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ct.ns_closure,
                                   movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                       2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
                                   1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                   2: lambda x, flag: domain.bc[flag].v_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                   2: {2: lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}

class PerturbedSurface_p(object):
    def __init__(self,waterdepth,amplitude):
        self.waterdepth=waterdepth
        self.amplitude=amplitude
    def uOfXT(self,x,t):
        d = ct.signedDistance(x, 0.)
        if d <= 0:
            return ct.pressure(x[0], x[1]-self.waterdepth, t, ct.h, ct.eps, ct.rho_0, ct.g, ct.k, 0.)+(ct.tank_dim[1]-(self.waterdepth+ct.eta(x[2], 0.)))*ct.rho_1*(-ct.g[1])
            # return (ct.tank_dim[1]-(self.waterdepth+ct.eta(x)))*ct.rho_1*(-ct.g[1])+((self.waterdepth+ct.eta(x))-x[1])*ct.rho_0*(-ct.g[1])
        else:
            return (ct.tank_dim[1] - x[1])*ct.rho_1*(-ct.g[1])

class AtRest(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(ct.water_depth,
                                          ct.water_amplitude), #[temp] from main we need depth and amplitude
                     1:AtRest(),
                     2:AtRest()}
