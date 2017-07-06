from proteus.default_n import *
from proteus import (FemTools,
                     Quadrature,
                     TimeIntegration,
                     NumericalFlux,
                     NonlinearSolvers,
                     LinearSolvers)
import moveMesh_p as physics
from proteus import Context
ct = Context.get()
nLevels = ct.nLevels
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = ct.restrictFineSolutionToAllMeshes
triangleOptions = ct.triangleOptions

timeIntegration = TimeIntegration.NoIntegration

femSpaces = {0:ct.basis,
             1:ct.basis}

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature


subgridError = None

massLumping = False

numericalFluxType = NumericalFlux.Stress_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother = None

fullNewtonFlag = True


matrix = SparseMatrix

if True:#cek hack, profiling, ct.opts.parallel:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py
else:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'mesh_'
linearSmoother = None
linearSolverConvergenceTest = 'r-true'
tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ct.mesh_nl_atol_res
nl_atol_res = ct.mesh_nl_atol_res
maxNonlinearIts = 4#should be linear
maxLineSearches = 0

conservativeFlux = None

auxiliaryVariables=[physics.fo]
