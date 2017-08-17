from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.default_n import *
import twp_navier_stokes_p as physics
from proteus.mprans import RANS2P
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions
adaptMesh = ct.adaptMesh
adaptMesh_nSteps = ct.adaptMesh_nSteps
adaptMesh_numIter = ct.adaptMesh_numIter
MeshAdaptMesh=ct.MeshAdaptMesh
useModel=ct.useModel
if ct.useHex or ct.structured:
    nnx = ct.nnx
    nny = ct.nny

    if ct.useHex:
        quad = True


#time stepping
runCFL = ct.runCFL
if ct.timeDiscretization=='vbdf':
    timeIntegration = TimeIntegration.VBDF
    timeOrder=2
    stepController  = StepControl.Min_dt_cfl_controller
elif ct.timeDiscretization=='flcbdf':
    timeIntegration = TimeIntegration.FLCBDF
    #stepController = FLCBDF_controller_sys
    stepController  = StepControl.Min_dt_cfl_controller
    time_tol = 10.0*ct.ns_nl_atol_res
    atol_u = {1:time_tol,2:time_tol}
    rtol_u = {1:time_tol,2:time_tol}
else:
    timeIntegration = TimeIntegration.BackwardEuler_cfl
    stepController  = StepControl.Min_dt_cfl_controller

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

massLumping       = False

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients=physics.coefficients,
                                   nd=nd,
                                   lag=ct.ns_lag_subgridError,
                                   hFactor=ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients=physics.coefficients,
                                       nd=nd,
                                       shockCapturingFactor=ct.ns_shockCapturingFactor,
                                       lag=ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = LinearSolvers.SimpleNavierStokes2D

matrix = LinearAlgebraTools.SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = LinearSolvers.PETSc
    levelLinearSolver      = LinearSolvers.PETSc
else:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
    schur_solver = ct.opts.schur_solver
    if schur_solver == 'Qp':
        linearSmoother=NavierStokes3D_Qp
    elif schur_solver == 'petsc_ASM':
        linearSmoother = petsc_ASM
    elif schur_solver == 'two_phase_Qp':
        linearSmoother=NavierStokes_TwoPhaseQp
    elif schur_solver == 'two_phase_LSC':
        linearSmoother=NavierStokes_TwoPhaseLSC
    elif schur_solver == 'two_phase_PCD':
        linearSmoother=NavierStokes_TwoPhasePCD
    elif schur_solver == 'LSC':
        linearSmoother=NavierStokes3D_LSC
    elif schur_solver == 'pcd':
        linearSmoother=NavierStokes3D_PCD
    elif schur_solver == 'selfp_proteus':
        linearSmoother = Schur_Sp
    elif schur_solver == 'selfp_petsc':
        linearSmoother = SimpleNavierStokes2D
    elif schur_solver == 'petsc_LU':
        linearSmoother=petsc_LU
    else:
        raise Exception, 'invalid solver type'
    linearSolverConvergenceTest = 'r-true'

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'rans2p_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
useEisenstatWalker = False#True
maxNonlinearIts = 50
maxLineSearches = 0
if ct.useHex:
    conservativeFlux = None
else:
    if ct.adaptMesh and ct.opts.parallel:
        conservativeFlux = None
    else:
        conservativeFlux = {0: 'pwl-bdm-opt'}

auxiliaryVariables = ct.domain.auxiliaryVariables['twp']
