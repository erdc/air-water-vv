from proteus.default_n import *
from vof_p import *
triangleOptions = ct.triangleOptions
timeIntegration = VBDF
timeOrder = ct.timeOrder
stepController  = Min_dt_cfl_controller
runCFL = ct.runCFL
femSpaces = {0:ct.basis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

massLumping       = False
numericalFluxType = VOF3P.NumericalFlux
conservativeFlux  = None
subgridError      = VOF3P.SubgridError(coefficients=coefficients,
                                       nd=ct.nd)
shockCapturing    = VOF3P.ShockCapturing(coefficients,
                                         ct.nd,
                                         shockCapturingFactor=ct.vof_shockCapturingFactor,
                                         lag=ct.vof_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = LinearSolvers.SparseMatrix

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'vof_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac      = 0.0
linTolFac   = 0.0
l_atol_res = 0.1*ct.vof_nl_atol_res
nl_atol_res = ct.vof_nl_atol_res

maxNonlinearIts = 50
maxLineSearches = 0

nLevels=ct.nLevels
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
auxiliaryVariables=[ct.height_gauges]
