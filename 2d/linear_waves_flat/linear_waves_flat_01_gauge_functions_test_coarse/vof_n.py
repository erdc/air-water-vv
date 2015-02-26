from proteus import *
ct = Context.get()
from vof_p import *
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
    time_tol = 10.0*ct.vof_nl_atol_res
    atol_u = {0:time_tol}
    rtol_u = {0:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:ct.basis}

massLumping       = False
numericalFluxType = VOF.NumericalFlux
conservativeFlux  = None
subgridError      = VOF.SubgridError(coefficients=coefficients,nd=nd)
shockCapturing    = VOF.ShockCapturing(coefficients,nd,shockCapturingFactor=ct.vof_shockCapturingFactor,lag=ct.vof_lag_shockCapturing)

fullNewtonFlag = True
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

linear_solver_options_prefix = 'vof_'
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac      = 0.0
nl_atol_res = ct.vof_nl_atol_res

linTolFac   = 0.0
l_atol_res = 0.1*ct.vof_nl_atol_res

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
