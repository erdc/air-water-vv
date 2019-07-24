from proteus.default_n import *
from proteus import (Context,
                     LinearAlgebraTools,
                     LinearSolvers,
                     NonlinearSolvers,
                     StepControl,
                     TimeIntegration
                     )
from proteus.mprans import RANS2P
import twp_navier_stokes_p as physics

#from broad_crested_weir import *
ct = Context.get()
mesh = ct.domain.MeshOptions

runCFL = ct.runCFL
if ct.timeDiscretization == 'vbdf':
    timeIntegration = TimeIntegration.VBDF
    timeOrder = 2
    stepController = StepControl.Min_dt_cfl_controller
elif ct.timeDiscretization == 'flcbdf':
    timeIntegration = TimeIntegration.FLCBDF
    #stepController = FLCBDF_controller_sys
    stepController = StepControl.Min_dt_cfl_controller
    time_tol = 10.0 * ct.ns_nl_atol_res
    atol_u = {1: time_tol, 2: time_tol}
    rtol_u = {1: time_tol, 2: time_tol}
else:
    timeIntegration = TimeIntegration.BackwardEuler_cfl
    stepController = StepControl.Min_dt_cfl_controller

# mesh options
nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
triangleOptions = mesh.triangleOptions

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

femSpaces = {0: ct.basis,
             1: ct.basis,
             2: ct.basis}

massLumping = False

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients=physics.coefficients,
                                   nd=ct.domain.nd,
                                   lag=ct.ns_lag_subgridError,
                                   hFactor=ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients=physics.coefficients,
                                       nd=ct.domain.nd,
                                       shockCapturingFactor=ct.ns_shockCapturingFactor,
                                       lag=ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

nonlinearSmoother = None

linearSmoother = LinearSolvers.SimpleNavierStokes2D

matrix = LinearAlgebraTools.SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = LinearSolvers.PETSc
    levelLinearSolver = LinearSolvers.PETSc
else:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver = LinearSolvers.KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU

linear_solver_options_prefix = 'rans2p_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest = 'rits-true'

tolFac = 0.0
nl_atol_res = ct.ns_nl_atol_res

linTolFac = 0.01
l_atol_res = 0.01 * ct.ns_nl_atol_res
useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = None#{0: 'pwl-bdm-opt'}

auxiliaryVariables = ct.domain.auxiliaryVariables['twp']
