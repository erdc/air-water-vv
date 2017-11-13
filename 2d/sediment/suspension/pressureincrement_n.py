from proteus import *
from proteus.default_n import *
from pressureincrement_p import *

mesh = domain.MeshOptions
triangleOptions = triangleOptions
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel

femSpaces = {0:ct.pbasis}

stepController=FixedStep

numericalFluxType = PresInc.NumericalFlux
matrix = LinearAlgebraTools.SparseMatrix

#numericalFluxType = PresInc.NumericalFlux
#subgridError = PresInc.SubgridError(coefficients,nd,lag=ct.ns_sed_lag_subgridError,hFactor=hFactor)
#shockCapturing = PresInc.ShockCapturing(coefficients,nd,ct.ns_sed_shockCapturingFactor,lag=ct.ns_sed_lag_shockCapturing)


#if openTop:
#    linearSmoother    = None
#    multilevelLinearSolver = LinearSolvers.LU
#    levelLinearSolver      = LinearSolvers.LU
#else:
#    

linear_solver_options_prefix = 'phi_'

if ct.useSuperlu:
    if ct.openTop is False:
        linearSmoother         = LinearSolvers.NavierStokesPressureCorrection # pure neumann laplacian solver
        multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        levelLinearSolver      = LinearSolvers.KSP_petsc4py
    else:
        multilevelLinearSolver = LinearSolvers.LU
        levelLinearSolver      = LinearSolvers.LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
    parallelPartitioningType    = parallelPartitioningType
    nLayersOfOverlapForParallel = nLayersOfOverlapForParallel
    nonlinearSmoother = None
    linearSmoother    = None
    if ct.openTop is False:
        linearSmoother         = LinearSolvers.NavierStokesPressureCorrection # pure neumann laplacian solver


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
periodicDirichletConditions=None

#conservativeFlux = {0:'pwl-bdm-opt'} #'point-eval','pwl-bdm-opt'
conservativeFlux=None
