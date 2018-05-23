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
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

# time stepping
timeIntegration = TimeIntegration.NoIntegration

# mesh options
nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
triangleOptions = mesh.triangleOptions
nn = ct.nn

from petsc4py import PETSc
OptDB = PETSc.Options()
OptDB.setValue("ksp_constant_null_space", 1)
OptDB.setValue("info", 1)

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

femSpaces = {0: ct.basis}

massLumping       = False
numericalFluxType = NumericalFlux.Diffusion_IIPG_exterior
conservativeFlux  = None

subgridError = None
shockCapturing = None

fullNewtonFlag = True
multilevelNonlinearSolver  = NonlinearSolvers.Newton
levelNonlinearSolver       = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother = None

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'mesh_'
linearSmoother = None
linearSolverConvergenceTest = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 1.0e-8
nl_atol_res = 1.0e-8
maxNonlinearIts = 4#should be linear
maxLineSearches = 0
