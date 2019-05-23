from proteus import *
from marin import *
from proteus.default_n import *
from twp_navier_stokes_p import *

from proteus import Context
ct = Context.get()

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis,
	     3:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NewtonNS
levelNonlinearSolver      = NewtonNS

nonlinearSmoother = None
linearSmoother    = SimpleNavierStokes3D

matrix = SparseMatrix

if usePETSc:    
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
    #schur_solver = 'two_phase_PCD' #'two_phase_PCD' #'selfp_petsc'
    schur_solver = ct.opts.schur_solver
    A_block_amg = ct.opts.A_block_amg
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
        linearSmootherOptions = (False,
                                 False,
                                 ct.opts.pcd_settings_lumped,
                                 0,
                                 laplace_null_space,
                                 A_block_amg) #(density_scaling, numerical_viscosity, lumped, chebyshev_its, velocity_block_preconditioner)
    elif schur_solver == 'LSC':
        linearSmoother=NavierStokes3D_LSC
    elif schur_solver == 'pcd':
        linearSmoother=NavierStokes3D_PCD
    elif schur_solver == 'selfp_proteus':
        linearSmoother = Schur_Sp
        linearSmootherOptions = (A_block_amg,) # (velocity_block_preconditioner)
    elif schur_solver == 'selfp_petsc':
        linearSmoother = SimpleNavierStokes3D
        linearSmootherOptions = (A_block_amg,) # (velocity_block_preconditioner)
    elif schur_solver == 'petsc_LU':
        linearSmoother=petsc_LU
    else:
        raise Exception, 'invalid solver type'
    linearSolverConvergenceTest = 'r-true'
    parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
    nLayersOfOverlapForParallel = 0
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*vof_nl_atol_res
nl_atol_res = ns_nl_atol_res
useEisenstatWalker = False#True
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}
