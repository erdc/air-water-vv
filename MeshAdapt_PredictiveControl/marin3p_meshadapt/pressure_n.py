from proteus.default_n import *
from pressure_p import *


triangleOptions = ct.triangleOptions

femSpaces = {0:ct.pbasis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

stepController=FixedStep

#matrix type
numericalFluxType = ConstantAdvection_exterior
matrix = SparseMatrix

linear_solver_options_prefix = 'pressure_'

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py
parallelPartitioningType = parallelPartitioningType
nLayersOfOverlapForParallel = nLayersOfOverlapForParallel
nonlinearSmoother = None
linearSmoother    = None

multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 0.1*ct.pressure_nl_atol_res
tolFac = 0.0
nl_atol_res = ct.pressure_nl_atol_res
maxLineSearches = 0
nonlinearSolverConvergenceTest      = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

periodicDirichletConditions=None

conservativeFlux=None
auxiliaryVariables=[ct.pressure_gauges]
