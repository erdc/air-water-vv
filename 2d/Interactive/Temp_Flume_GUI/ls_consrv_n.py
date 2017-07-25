from proteus.default_n import *
from proteus import (#Context,
                     LinearAlgebraTools,
                     LinearSolvers,
                     NonlinearSolvers,
                     StepControl,
                     TimeIntegration
                     )
import ls_consrv_p as physics
import Context
ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

#time stepping
runCFL = ct.runCFL
timeIntegrator  = TimeIntegration.ForwardIntegrator
timeIntegration = TimeIntegration.NoIntegration

#mesh options
nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
triangleOptions = mesh.triangleOptions

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature
femSpaces = {0:ct.basis}

subgridError      = None
massLumping       = False
numericalFluxType = DoNothing
conservativeFlux  = None
shockCapturing    = None

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = LinearAlgebraTools.SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'mcorr_'
nonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest  = 'r-true'
levelNonlinearSolverConvergenceTest  = 'r'

tolFac = 0.0
nl_atol_res = ct.mcorr_nl_atol_res

linTolFac = 0.0
l_atol_res = 0.001*ct.mcorr_nl_atol_res
useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0

auxiliaryVariables = ct.domain.auxiliaryVariables['ls_consrv']
