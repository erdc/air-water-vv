from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.default_n import *
from proteus import Context
from twp_navier_stokes_p import *

ct = Context.get()
nd = ct.domain.nd
runCFL = ct.runCFL
nLevels = ct.nLevels
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = ct.restrictFineSolutionToAllMeshes
triangleOptions = ct.triangleOptions
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_controller


if nd == 2:
    femSpaces = {0:ct.basis,
	         1:ct.basis,
	         2:ct.basis}
if nd == 3:
    femSpaces = {0:ct.basis,
                 1:ct.basis,
                 2:ct.basis,
                 3:ct.basis}
    
massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,
                                   ct.domain.nd,
                                   lag = ct.ns_lag_subgridError,
                                   hFactor = ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,
                                       ct.domain.nd,
                                       ct.ns_shockCapturingFactor,
                                       lag = ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
if nd == 2:
    linearSmoother = LinearSolvers.SimpleNavierStokes2D
elif nd == 3:
    linearSmoother = LinearSolvers.SimpleNavierStokes2D
    
matrix = LinearAlgebraTools.SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = LinearSolvers.PETSc
    levelLinearSolver      = LinearSolvers.PETSc
else:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest = 'r-true'

tolFac = 0.0
linTolFac = 0.00001
l_atol_res = 0.00001*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}

#auxiliaryVariables = [ct.pointGauges]
