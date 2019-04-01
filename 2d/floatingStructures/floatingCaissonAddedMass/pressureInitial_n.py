from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
try:
    from .pressureInitial_p import *
except:
    from pressureInitial_p import *

timeIntegration=NoIntegration

triangleOptions = triangleOptions

femSpaces = {0:pbasis}

stepController=FixedStep

#numericalFluxType = NumericalFlux.ConstantAdvection_Diffusion_SIPG_exterior #weak boundary conditions (upwind ?)
matrix = LinearAlgebraTools.SparseMatrix

if useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
    parallelPartitioningType = parallelPartitioningType
    nLayersOfOverlapForParallel = nLayersOfOverlapForParallel
    nonlinearSmoother = None
    linearSmoother    = None

numericalFluxType = NumericalFlux.ConstantAdvection_exterior

linear_solver_options_prefix = 'pinit_'

multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 1.0e-10
tolFac = 0.0
nl_atol_res = 1.0e-10
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
periodicDirichletConditions=None

conservativeFlux=None
