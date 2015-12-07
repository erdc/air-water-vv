from proteus.default_p import *
from proteus.mprans import RANS2P
import numpy as np
from proteus import Context
ct = Context.get()


domain = ct.domain
nd = domain.nd

genMesh = ct.genMesh
movingDomain = ct.movingDomain
T = ct.T  # might not be necessary


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
    if ct.movingDomain:
        Closure_0_model += 1
        Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

    
# for absorption zones (defined as regions)
# (!) should be done with regionFlags but all regions have different flags so far
porosityTypes = np.ones(len(ct.domain.regions)+1)
dragAlphaTypes = np.zeros(len(ct.domain.regions)+1)
dragBetaTypes = np.zeros(len(ct.domain.regions)+1)
epsFact_solid = np.zeros(len(ct.domain.regions)+1)
for auxvar in ct.domain.auxiliaryVariables:
    if hasattr(auxvar.shape, 'regions'):
        if auxvar.shape.regions is not None:
            auxvar_porosityTypes = np.ones(len(auxvar.shape.regions))
            auxvar_dragAlphaTypes = np.zeros(len(auxvar.shape.regions))
            auxvar_dragBetaTypes = np.zeros(len(auxvar.shape.regions))
            auxvar_epsFact_solid = np.zeros(len(auxvar.shape.regions))
            if hasattr(auxvar.shape, 'porosityTypes'):
                auxvar_porosityTypes[:] = auxvar.shape.porosityTypes
                auxvar_dragAlphaTypes[:] = auxvar.shape.dragAlphaTypes
                auxvar_dragBetaTypes[:] = auxvar.shape.dragBetaTypes
                auxvar_epsFact_solid[:] = auxvar.shape.epsFact_solid
            i0 = auxvar.shape._snr+1
            i1 = i0 + len(auxvar.shape.regions)
            porosityTypes[i0:i1] = auxvar_porosityTypes
            dragAlphaTypes[i0:i1] = auxvar_dragAlphaTypes
            dragBetaTypes[i0:i1] = auxvar_dragBetaTypes
            epsFact_solid[i0:i1] = auxvar_epsFact_solid
if not np.any(dragAlphaTypes):  # checking if all values are 0.
    porosityTypes = None
    dragAlphaTypes = None
    dragBetaTypes = None
    epsFact_solid = None

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
                                   movingDomain=ct.movingDomain,
                                   porosityTypes=porosityTypes,
                                   dragAlphaTypes=dragAlphaTypes,
                                   dragBetaTypes=dragBetaTypes,
                                   epsFact_solid=epsFact_solid,
                                   barycenters=ct.domain.barycenters)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet,
                       1: lambda x, flag: domain.bc[flag].u_dirichlet,
                       2: lambda x, flag: domain.bc[flag].v_dirichlet}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective,
                                   1: lambda x, flag: domain.bc[flag].u_advective,
                                   2: lambda x, flag: domain.bc[flag].v_advective}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: domain.bc[flag].u_diffusive},
                                   2: {2: lambda x, flag: domain.bc[flag].v_diffusive}}

if nd == 3:
    dirichletConditions[3] = lambda x, flag: domain.bc[flag].w_dirichlet
    advectiveFluxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_advective
    diffusiveFluxBoundaryConditions[3] = {3: lambda x, flag: domain.bc[flag].w_diffusive}

    
def signedDistance(x):
    phi_x = x[0]-ct.waterLine_x
    phi_z = x[1]-ct.waterLine_z 
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)
    
class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(ct.tank_dim[1] - self.waterLevel)*ct.rho_1*ct.g[1] - (self.waterLevel - x[1])*ct.rho_0*ct.g[1]
        else:
            return -(ct.tank_dim[1] - self.waterLevel)*ct.rho_1*ct.g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0


initialConditions = {0:PerturbedSurface_p(ct.waterLine_z),
                     1:AtRest(),
                     2:AtRest()}
