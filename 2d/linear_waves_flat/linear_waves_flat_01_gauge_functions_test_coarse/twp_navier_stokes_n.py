from proteus import *
from twp_navier_stokes_p import *
from proteus.default_n import *   
ct=Context.get()
triangleOptions = ct.triangleOptions
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
parallelPartitioningType=ct.parallelPartitioningType
nLevels = ct.nLevels
elementQuadrature=ct.elementQuadrature
elementBoundaryQuadrature=ct.elementBoundaryQuadrature
runCFL=ct.runCFL

if ct.timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
elif ct.timeDiscretization=='flcbdf':
    timeIntegration = FLCBDF
    #stepController = FLCBDF_controller_sys
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*ct.ns_nl_atol_res
    atol_u = {1:time_tol,2:time_tol}
    rtol_u = {1:time_tol,2:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:ct.basis,
	     1:ct.basis,
	     2:ct.basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ct.ns_lag_subgridError,hFactor=ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ct.ns_shockCapturingFactor,lag=ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None

linearSmoother    = SimpleNavierStokes2D

matrix = SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'rans2p_'
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}

auxiliaryVariables=[ct.pointGauges]
