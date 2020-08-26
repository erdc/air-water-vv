from proteus import *
from proteus.default_n import *
from pressureincrement_p import *

mesh = domain.MeshOptions
triangleOptions = triangleOptions
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel

femSpaces = {0:pbasis}

stepController=FixedStep

numericalFluxType = PresInc.NumericalFlux
#numericalFluxType = NumericalFlux.ConstantAdvection_exterior
matrix = LinearAlgebraTools.SparseMatrix
conservativeFlux = {0:'point-eval'} #'point-eval','pwl-bdm-opt'

#numericalFluxType = None
#matrix = SparseMatrix
#conservativeFlux=None

linear_solver_options_prefix = 'phi_'

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
    parallelPartitioningType    = parallelPartitioningType
    nLayersOfOverlapForParallel = nLayersOfOverlapForParallel
    nonlinearSmoother = None
    linearSmoother    = None


multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 0.01*phi_nl_atol_res
tolFac = 0.0
nl_atol_res = phi_nl_atol_res
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
