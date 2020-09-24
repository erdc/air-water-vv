from proteus import StepControl
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans import VOS3P
from proteus import Context
from proteus import *
import sheetflowBC as sfbc

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

LevelModelType = VOS3P.LevelModel

coefficients = VOS3P.Coefficients(V_model=ct.V_model,
                                  SED_model=ct.SED_model,
                                  ME_model=ct.VOS_model,
                                  checkMass=False,
                                  useMetrics=ct.useMetrics,
                                  epsFact=ct.epsFact_vos,
                                  sc_uref=ct.vos_sc_uref,
                                  sc_beta=ct.vos_sc_beta,
                                  movingDomain=ct.movingDomain,
                                  vos_function=ct.vos_function,
                                  global_max_u=ct.opts.vos_limiter,
                                  STABILIZATION_TYPE=4,
                                  LUMPED_MASS_MATRIX=True,
                                  cE = 1.,
                                  cK=0.
                                  
                                  )

manualbc = ct.manualbc

if manualbc == True:
	parallelPeriodic=sfbc.vos_parallelPeriodic
	periodicDirichletConditions 	= sfbc.vos_periodic
	dirichletConditions		= sfbc.vos_dirichlet
	advectiveFluxBoundaryConditions = sfbc.vos_advective
	diffusiveFluxBoundaryConditions = sfbc.vos_diffusive
else:
	dirichletConditions = {0: lambda x, flag: domain.bc[flag].vos_dirichlet.init_cython()}
	advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vos_advective.init_cython()}
	diffusiveFluxBoundaryConditions = {0: {}}

initialConditions  = {0:ct.Suspension}

def dbc(x,flag):
	if sfbc.onLeft(x) or sfbc.onRight(x):
		return ct.Suspension.uOfXT

dirichletConditions		= {0:dbc}

def abc(x,flag):
	if sfbc.onLeft(x) or sfbc.onRight(x):
		return None
	else:
		return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0: abc}
diffusiveFluxBoundaryConditions = {0: {}}
