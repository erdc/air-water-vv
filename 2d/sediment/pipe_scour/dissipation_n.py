from proteus import StepControl
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
import dissipation_p as physics
from proteus import Context
from proteus.mprans import Dissipation

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

triangleOptions = mesh.triangleOptions
elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature
timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_cfl_controller

femSpaces = {0:ct.basis}

massLumping       = False
numericalFluxType = Dissipation.NumericalFlux
conservativeFlux  = None
subgridError      = Dissipation.SubgridError(coefficients=physics.coefficients,nd=ct.nd)
shockCapturing    = Dissipation.ShockCapturing(physics.coefficients,ct.nd,shockCapturingFactor=ct.dissipation_shockCapturingFactor,
                                         lag=ct.dissipation_lag_shockCapturing)

fullNewtonFlag  = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None
#printNonlinearSolverInfo = True
matrix = SparseMatrix
if not ct.useOldPETSc and not ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py
else:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'dissipation_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac = 0.0
linTolFac =0.0
l_atol_res = 0.001*ct.dissipation_nl_atol_res
nl_atol_res = ct.dissipation_nl_atol_res
useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0

