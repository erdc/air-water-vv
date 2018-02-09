from proteus.default_n import *
from pressureInitial_p import *

triangleOptions = ct.triangleOptions

femSpaces = {0:ct.pbasis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

stepController=FixedStep

matrix = SparseMatrix

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
nonlinearSmoother = None
linearSmoother    = None

numericalFluxType = ConstantAdvection_exterior

linear_solver_options_prefix = 'pinit_'

multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 0.01*ct.phi_nl_atol_res
tolFac = 0.0
nl_atol_res = ct.phi_nl_atol_res
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
periodicDirichletConditions=None
conservativeFlux=None
nLevels=ct.nLevels
