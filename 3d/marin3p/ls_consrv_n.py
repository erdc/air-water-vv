from proteus.default_n import *
from ls_consrv_p import *

timeIntegration = NoIntegration

femSpaces = {0:ct.basis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

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

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'mcorr_'
linearSolverConvergenceTest  = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.mcorr_nl_atol_res
nl_atol_res = ct.mcorr_nl_atol_res

maxNonlinearIts = 50
maxLineSearches = 0
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
nLevels=ct.nLevels
