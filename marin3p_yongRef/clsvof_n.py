from proteus import *
from proteus.default_n import *
from clsvof_p import *
from marin import *

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = CLSVOFNewton
fullNewtonFlag = True
updateJacobian = True
timeIntegration = BackwardEuler_cfl

eps_tolerance_clsvof=False
if eps_tolerance_clsvof:
    nl_atol_res = 1E-12
else:
    nl_atol_res=ct.clsvof_nl_atol_res
#
l_atol_res = nl_atol_res
tolFac=0.
maxNonlinearIts = 100
stepController = Min_dt_controller

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature
femSpaces = {0:ct.pbasis}

#numericalFluxType = CLSVOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

matrix = SparseMatrix
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'clsvof_'
linearSolverConvergenceTest = 'r-true'
