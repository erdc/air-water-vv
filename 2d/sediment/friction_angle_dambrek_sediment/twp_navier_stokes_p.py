from proteus import *
from proteus.default_p import *
from tank import *
from proteus.mprans import RANS3PF
from proteus import Context

ct = Context.get()

LevelModelType = RANS3PF.LevelModel

if ct.sedimentDynamics:
    VOS_model=0
    VOF_model=1
    LS_model=2
    RD_model=3
    MCORR_model=4
    SED_model=5
    V_model=6
    PINC_model=7
    PRESSURE_model=8
    PINIT_model=9
else:
    VOS_model=None
    SED_model=None
    VOF_model=0
    LS_model=1
    RD_model=2
    MCORR_model=3
    V_model=4
    PINC_model=5
    PRESSURE_model=6
    PINIT_model=7
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
                                    ME_model=V_model,
                                    PRESSURE_model=PRESSURE_model,
                                    SED_model=SED_model,
                                    VOF_model=VOF_model,
                                    VOS_model=VOS_model,
                                    LS_model=LS_model,
                                    Closure_0_model=Closure_0_model,
                                    Closure_1_model=Closure_1_model,
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
                                    PSTAB=ct.opts.PSTAB,
                                    aDarcy = sedClosure.aDarcy,
                                    betaForch = sedClosure.betaForch,
                                    grain = sedClosure.grain,
                                    packFraction = sedClosure.packFraction,
                                    maxFraction = sedClosure.maxFraction,
                                    frFraction = sedClosure.frFraction,
                                    sigmaC = sedClosure.sigmaC,
                                    C3e = sedClosure.C3e,
                                    C4e = sedClosure.C4e,
                                    eR = sedClosure.eR,
                                    fContact = sedClosure.fContact,
                                    mContact = sedClosure.mContact,
                                    nContact = sedClosure.nContact,
                                    angFriction = sedClosure.angFriction,
                                    vos_function = ct.vos_function,
                                    vos_limiter = ct.sedClosure.vos_limiter,
                                    mu_fr_limiter = ct.sedClosure.mu_fr_limiter,
                                    )

dirichletConditions = {0: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                   1: lambda x, flag: domain.bc[flag].v_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                   1: {1: lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}

if nd == 3:
    dirichletConditions[2] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
    advectiveFluxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
    diffusiveFluxBoundaryConditions[2] = {2: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
