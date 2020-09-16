from proteus.default_n import *
from ls_p import *
triangleOptions = ct.triangleOptions
timeIntegration = VBDF
timeOrder = ct.timeOrder
stepController  = Min_dt_cfl_controller
femSpaces = {0:ct.basis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

massLumping       = False
conservativeFlux  = None
numericalFluxType = NCLS3P.NumericalFlux
subgridError      = NCLS3P.SubgridError(coefficients,
                                      ct.nd)
shockCapturing    = NCLS3P.ShockCapturing(coefficients,
                                        ct.nd,
                                        shockCapturingFactor=ct.ls_shockCapturingFactor,
                                        lag=ct.ls_lag_shockCapturing)

fullNewtonFlag  = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'ncls_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.01*ct.ls_nl_atol_res
nl_atol_res = ct.ls_nl_atol_res

maxNonlinearIts = 50
maxLineSearches = 0

parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
nLevels=ct.nLevels
