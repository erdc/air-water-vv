from proteus import *
from proteus.default_p import *
from circle_caisson2D_oscillationIB import *
from proteus.mprans import RANS3PF

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


mylist_of_sdfs=[];

coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
                                    sigma=0.0,
                                    USE_SUPG = True,
                                    ARTIFICIAL_VISCOSITY=2,
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
                                    dragAlpha=0.0,
                                    PSTAB=0.0,
                                    nParticles=1,
                                    particle_epsFact=0.33,
                                    particle_alpha=0,
                                    particle_beta=0,
                                    particle_penalty_constant=100.0,
                                    particle_sdfList=[],
                                    particle_velocityList=[],
                                    use_sbm=0,
                                    use_ball_as_particle = 1,
                                    ball_center=ball_center,
                                    ball_radius=ball_radius,
                                    ball_velocity=ball_velocity,
                                    ball_angular_velocity=ball_angular_velocity)
#===============================================================================
# BC
#===============================================================================

def getDBC_u(x,flag):
    if onBoundary(x):
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if onBoundary(x):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

def getAFBC_u(x,flag):
    return None

def getAFBC_v(x,flag):
    return None

def getDFBC_u(x,flag):
    return None

def getDFBC_v(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

class AtRest:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
