from proteus.default_n import *
from rans2p_p import *
from marin import *

timeIntegration = BackwardEuler_cfl#VBDF
timeOrder=1
stepController  = Min_dt_cfl_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis,
	     3:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = SimpleNavierStokes3D

matrix = SparseMatrix

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*vof_nl_atol_res
nl_atol_res = ns_nl_atol_res
useEisenstatWalker = False#True
maxNonlinearIts = 100
maxLineSearches = 0
if not ns_forceStrongDirichlet:
    conservativeFlux = None#{0:'pwl-bdm-opt'}
auxiliaryVariables=[pressure_gauges]
