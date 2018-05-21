from proteus.default_n import *
from pressureincrement_p import *

triangleOptions = ct.triangleOptions

femSpaces = {0:ct.pbasis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

stepController=FixedStep

numericalFluxType = PresInc.NumericalFlux
matrix = SparseMatrix

parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver = KSP_petsc4py

if ct.openTop:
    nonlinearSmoother = None
    linearSmoother    = None
else:
    nonlinearSmoother = None
    linearSmoother    = NavierStokesPressureCorrection # pure neumann laplacian solver

linear_solver_options_prefix = 'phi_'

multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = ct.phi_nl_atol_res
tolFac = 0.0
nl_atol_res = ct.phi_nl_atol_res
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
periodicDirichletConditions=None

#conservativeFlux = {0:'point-eval'} #'point-eval','pwl-bdm-opt'
conservativeFlux=None
nLevels=ct.nLevels
