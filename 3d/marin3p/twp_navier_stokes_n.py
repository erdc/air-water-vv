from proteus.default_n import *
from twp_navier_stokes_p import *

triangleOptions = ct.triangleOptions
timeIntegration = VBDF
timeOrder = ct.timeOrder
stepController  = Min_dt_cfl_controller
runCFL = ct.runCFL
femSpaces = {0:ct.basis,
	     1:ct.basis,
	     2:ct.basis}

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS3PF.NumericalFlux
subgridError = RANS3PF.SubgridError(coefficients,
                                    nd=ct.nd,
                                    lag=ct.ns_lag_subgridError,
                                    hFactor=ct.hFactor)
shockCapturing = RANS3PF.ShockCapturing(coefficients,
                                        ct.nd,
                                        ct.ns_shockCapturingFactor,
                                        ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.001*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = None#{0:'pwl-bdm-opt'}
nLevels=ct.nLevels
parallelPartitioningType = ct.parallelPartitioningType
nLayersOfOverlapForParallel = ct.nLayersOfOverlapForParallel
