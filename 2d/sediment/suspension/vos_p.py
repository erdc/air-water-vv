from proteus import StepControl
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans import VOS3P
from proteus import Context
from proteus import *
import tank_so
ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

LevelModelType = VOS3P.LevelModel


coefficients = VOS3P.Coefficients(LS_model=None,
                                  V_model=tank_so.SED_model,
                                  SED_model=tank_so.SED_model,
                                  RD_model=None,
                                  ME_model=tank_so.VOS_model,
                                  VOF_model=tank_so.VOF_model,
                                  checkMass=False,
                                  useMetrics=ct.useMetrics,
                                  epsFact=ct.epsFact_vos,
                                  sc_uref=ct.vos_sc_uref,
                                  sc_beta=ct.vos_sc_beta,
                                  movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: domain.bc[flag].vos_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vos_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {}}


initialConditions  = {0:ct.Suspension}
