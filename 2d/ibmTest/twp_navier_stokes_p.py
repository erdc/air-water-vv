from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import RANS3PF
name="rans3p"
LevelModelType = RANS3PF.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
                                    sigma=0.0,
                                    rho_0 = rho_0,
                                    nu_0 = nu_0,
                                    rho_1 = rho_1,
                                    nu_1 = nu_1,
                                    g=g,
                                    nd=nd,
                                    ME_model=4,
                                    PRESSURE_model=6,
                                    SED_model=None,
                                    VOF_model=0,
                                    VOS_model=None,
                                    LS_model=1,
                                    Closure_0_model=None,
                                    Closure_1_model=None,
                                    epsFact_density=epsFact_density,
                                    stokes=False,
                                    useVF=useVF,
                                    useRBLES=useRBLES,
                                    useMetrics=useMetrics,
                                    eb_adjoint_sigma=1.0,
                                    eb_penalty_constant=weak_bc_penalty_constant,
                                    forceStrongDirichlet=ns_forceStrongDirichlet,
                                    turbulenceClosureModel=ns_closure,
                                    movingDomain=movingDomain,
                                    dragAlpha=dragAlpha,
                                    PSTAB=1.0,
                                    nParticles=1,
                                    particle_epsFact=2.0,
                                    particle_alpha=1e6,
                                    particle_beta=1e6,
                                    particle_penalty_constant=1e16,
                                    particle_sdfList=[particle_sdf],
                                    particle_velocityList=[particle_vel])




dirichletConditions = {0: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}



advectiveFluxBoundaryConditions =  {0: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                    1: lambda x, flag: domain.bc[flag].v_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {0:lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                   1: {1:lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}

fluxBoundaryConditions = {0: 'mixedFlow',
                          1: 'mixedFlow'}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}

