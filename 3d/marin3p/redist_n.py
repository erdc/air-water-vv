from proteus.default_n import *
from redist_p import *

nl_atol_res = ct.rd_nl_atol_res
tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.rd_nl_atol_res
useEisenstatWalker = False#True

if ct.redist_Newton:
    timeIntegration = NoIntegration
    stepController = Newton_controller
    maxNonlinearIts = 25
    maxLineSearches = 10
    nonlinearSolverConvergenceTest = 'r'
    levelNonlinearSolverConvergenceTest = 'r'
    linearSolverConvergenceTest = 'r-true'
else:
    timeIntegration = BackwardEuler_cfl
    stepController = RDLS3P.PsiTC
    runCFL=0.5
    psitc['nStepsForce']=6
    psitc['nStepsMax']=25
    psitc['reduceRatio']=2.0
    psitc['startRatio']=1.0
    rtol_res[0] = 0.0
    atol_res[0] = ct.rd_nl_atol_res
    maxNonlinearIts = 1
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
    linearSolverConvergenceTest = 'r-true'

femSpaces = {0:ct.basis}
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature
       
massLumping       = False
numericalFluxType = NumericalFlux.DoNothing    
conservativeFlux  = None
subgridError      = RDLS3P.SubgridError(coefficients,
                                      ct.nd)
shockCapturing    = RDLS3P.ShockCapturing(coefficients,
                                        ct.nd,
                                        shockCapturingFactor=ct.rd_shockCapturingFactor,
                                        lag=ct.rd_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver  = NonlinearSolvers.Newton
levelNonlinearSolver       = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'rdls_'

parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
