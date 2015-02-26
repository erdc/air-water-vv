from proteus import *
from ls_p import *
from proteus.default_n import *   
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
    #stepController = FLCBDF_controller
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*ct.ls_nl_atol_res
    atol_u = {0:time_tol}
    rtol_u = {0:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:ct.basis}

massLumping       = False
conservativeFlux  = None
numericalFluxType = NCLS.NumericalFlux
subgridError      = NCLS.SubgridError(coefficients,nd)
shockCapturing    = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=ct.ls_shockCapturingFactor,lag=ct.ls_lag_shockCapturing)

fullNewtonFlag  = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

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

linear_solver_options_prefix = 'ncls_'
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac = 0.0
nl_atol_res = ct.ls_nl_atol_res

linTolFac = 0.0
l_atol_res = 0.1*ct.ls_nl_atol_res

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0

