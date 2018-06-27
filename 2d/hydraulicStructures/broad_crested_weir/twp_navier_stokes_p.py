from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions
genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = RANS2P.LevelModel
if ct.useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if ct.useRANS >= 1:
    Closure_0_model = 5
    Closure_1_model = 6
    if ct.useOnlyVF:
        Closure_0_model = 2
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
                                   nd=domain.nd,
                                   ME_model=0,
                                   VF_model=1,
                                   LS_model=LS_model,
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
                                   NONCONSERVATIVE_FORM=1.0,
                                   movingDomain=ct.movingDomain)

if ct.ns_forceStrongDirichlet:
    walls = [ct.domain.boundaryTags[f] for f in ['tank2D1_y+','tank2D1_y-','tank2D1_obstacle1']]
    inflow = [ct.domain.boundaryTags['tank2D1_x-']]
    outflow = [ct.domain.boundaryTags['tank2D1_x+']]
    if ct.opts.air_vent:
        outflow.append(domain.boundaryTags['tank2D1_airvent'])
    def get_p_Dirichlet(x, flag):
        from math import fabs
        if flag in outflow:
            return lambda x,t: fabs(ct.g[1])*(
                ct.rho_1*(ct.tank_dim[1] - x[1])
                +
                (ct.rho_0 - ct.rho_1)*max(0,ct.opts.outflow_level - x[1]))

    def get_u_Dirichlet(x, flag):
        if flag in walls:
            return lambda x,t: 0.0
        elif flag in inflow:
            return ct.twpflowVelocity_u

    def get_v_Dirichlet(x, flag):
        if flag in walls:
            return lambda x,t: 0.0
        elif flag in inflow:
            return lambda x,t: 0.0
        elif flag in outflow:
            return lambda x,t: 0.0

    dirichletConditions = {
        0: get_p_Dirichlet,
        1: get_u_Dirichlet,
        2: get_v_Dirichlet
    }

    def get_mass_flux(x, flag):
        if flag in walls:
            return lambda x,t: 0.0
        elif flag in inflow:
            return lambda x,t: -ct.twpflowVelocity_u(x,t)
        elif flag == 0:
            return lambda x,t: 0.0

    def get_mom_u_adv_flux(x, flag):
        if flag == 0:
            return lambda x,t: 0.0
    
    def get_mom_v_adv_flux(x, flag):
        if flag == 0:
            return lambda x,t: 0.0

    advectiveFluxBoundaryConditions = {
        0: get_mass_flux,
        1: get_mom_u_adv_flux,
        2: get_mom_v_adv_flux
    }
    
    def get_mom_u_diff_flux(x, flag):
        if flag in outflow:
            return lambda x,t: 0.0
        elif flag == 0:
            return lambda x,t: 0.0

    def get_mom_v_diff_flux(x, flag):
        if flag == 0:
            return lambda x,t: 0.0

    diffusiveFluxBoundaryConditions = {
        0: {},
        1: {1: get_mom_u_diff_flux},
        2: {2: get_mom_v_diff_flux}
    }
else:
    dirichletConditions = {
        0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
        1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
        2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()
    }
    
    advectiveFluxBoundaryConditions = {
        0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
        1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
        2: lambda x, flag: domain.bc[flag].v_advective.init_cython()
    }
    
    diffusiveFluxBoundaryConditions = {
        0: {},
        1: {1: lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
        2: {2: lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}
    }
    
class PerturbedSurface_p:
    def __init__(self, waterLevel):
        self.waterLevel = waterLevel

    def uOfXT(self, x, t):
        if ct.signedDistance(x) < 0:
            return -(ct.tank_dim[1] - self.waterLevel) * ct.rho_1 * ct.g[1] \
                   - (self.waterLevel -x[1]) * ct.rho_0 * ct.g[1]
        else:
            return -(ct.tank_dim[1] - self.waterLevel) * ct.rho_1 * ct.g[1]
        

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
    

class initialVelocity_u:
    def __init__(self, waterLevel):
        self.waterLevel = waterLevel
    def uOfXT(self,x,t):
        if x[0]<ct.waterLine_x and x[nd-1]<self.waterLevel:
            return ct.opts.inflow_velocity
        elif x[0] < ct.tank_dim[0]+ct.tank_sponge[1]:
            return ct.twpflowVelocity_u_D(x, 0)
        else: 
            return 0.0


initialConditions = {0: PerturbedSurface_p(ct.waterLine_z),
                     1: AtRest(),
                     2: AtRest()}
