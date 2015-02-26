from proteus import *
ct = Context.get()
from ls_consrv_p import *
from proteus.default_n import *   
triangleOptions = ct.triangleOptions
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
parallelPartitioningType=ct.parallelPartitioningType
nLevels = ct.nLevels
elementQuadrature=ct.elementQuadrature
elementBoundaryQuadrature=ct.elementBoundaryQuadrature
runCFL=ct.runCFL

timeIntegrator  = ForwardIntegrator
timeIntegration = NoIntegration

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

matrix = SparseMatrix

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
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest  = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.mcorr_nl_atol_res
nl_atol_res = ct.mcorr_nl_atol_res
useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
