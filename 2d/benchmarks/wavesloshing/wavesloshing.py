"""
Wavesloshing Problem
"""
import numpy as np
from math import cos
from proteus import (Domain, Context,
                     FemTools as ft,
                     #SpatialTools as st,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D

opts=Context.Options([
    # predefined test cases
    ("water_depth_fraction", 0.5, "Depth of the water relative to the tank height"),
    ("amplitude_fraction", 0.02, "Amplitude of the sloshing relative to tank height"),
    # tank
    ("tank_dim", (3.141 , 3.141), "Dimensions of the tank"),
    #gravity
    ("g",(0,-9.81,0), "Gravity vector"),
    # probe dx
    ("dxProbe",0.25, "Probe spacing"),
    # refinement
    ("refinement", 16 ,"Refinement level"),
    ("cfl", 0.33 ,"Target cfl"),
    # run time
    ("T", 10.0 ,"Simulation time"),
    ("dt_fixed", 0.01, "Fixed time step"),
    ("dt_init", 0.001 ,"Maximum initial time step"),
    # run details
    ("structured", False, "Use a structured mesh"),
    ("gen_mesh", True ,"Generate new mesh"),
    ("parallel", True ,"Run in parallel")])

# ----- CONTEXT ------ #

# tank
tank_dim = opts.tank_dim

# water
water_depth = opts.water_depth_fraction * tank_dim[0]
water_amplitude = opts.amplitude_fraction * tank_dim[0]

##########################################
#     Discretization Input Options       #
##########################################

#[temp] temporary location
backgroundDiffusionFactor = 0.01 #[temp]

refinement = opts.refinement
genMesh = opts.gen_mesh
movingDomain = False
checkMass = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega

# ----- INPUT CHECKS ----- #
if spaceOrder not in [1,2]:
    raise ValueError("INVALID: spaceOrder(" + str(spaceOrder) + ")")

if useRBLES not in [0.0, 1.0]:
    raise ValueError("INVALID: useRBLES(" + str(useRBLES) + ")")

if useMetrics not in [0.0, 1.0]:
    raise ValueError("INVALID: useMetrics(" + str(useMetrics) + ")")

# ----- DISCRETIZATION ----- #
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
        basis=ft.C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd,2)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd-1,2)
    else:
    	 basis=ft.C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = ft.SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd-1,3)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
        basis=ft.C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = ft.CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = ft.CubeGaussQuadrature(nd-1,4)
    else:
        basis=ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = ft.SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd-1,4)

##########################################
# Numerical Options and Other Parameters #
##########################################

weak_bc_penalty_constant = 100.0
nLevels = 1

# ----- PHYSICAL PROPERTIES ----- #

# Water
rho_0 = 998.2
nu_0 = 1.004e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

# Surface Tension
sigma_01 = 0.0

# Gravity
g = opts.g

# ----- TIME STEPPING & VELOCITY----- #

T = opts.T
dt_fixed = opts.dt_fixed
dt_init = min(0.1 * dt_fixed, opts.dt_init)
runCFL = opts.cfl
nDTout = int(round(T / dt_fixed))

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# if useHex:
#     nnx=4*Refinement+1
#     nny=2*Refinement+1
#     hex=True
#     domain = Domain.RectangularDomain(tank_dim)
# else:
#     if structured:
#         nnx=4*Refinement
#         nny=2*Refinement
#         domain = Domain.RectangularDomain(tank_dim)
#         boundaryTags = domain.boundaryTags
#[temp] This is the same structured setup as in the dambreak problems

# ----- TANK ----- #

tank = Tank2D(domain, tank_dim)

# ----- GAUGES ----- #

#[temp] pressure (twp) and so on at these points:
# pointGauges = PointGauges(gaugeLocations={'pointGauge_pressure':(L[0]/2.0,L[0]/2.0,0.0)})
# lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825':((L[0]/2.0,0.0,0.0),(L[0]/2.0,L[1],0.0))},linePoints=20)
# lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)
# ----- EXTRA BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

# ----- MESH CONSTRUCTION ----- #

he = tank_dim[0] / float(4 * refinement - 1)
domain.MeshOptions.he = he
st.assembleDomain(domain)

# ----- STRONG DIRICHLET ----- #

ns_forceStrongDirichlet = False

# ----- NUMERICAL PARAMETERS ----- #
if useMetrics:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_viscosity = epsFact_curvature \
        = epsFact_vof = ecH \
        = epsFact_consrv_dirac = epsFact_density \
        = 1.5
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = epsFact_viscosity = epsFact_curvature \
        = epsFact_vof = ecH \
        = epsFact_consrv_dirac = epsFact_density \
        = 1.5
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

# ----- NUMERICS: TOLERANCES ----- #
ns_nl_atol_res = max(1.0e-8,0.01*he**2)
vof_nl_atol_res = max(1.0e-8,0.01*he**2)
ls_nl_atol_res = max(1.0e-8,0.01*he**2)
rd_nl_atol_res = max(1.0e-8,0.01*he)
mcorr_nl_atol_res = max(1.0e-8,0.01*he**2)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2)

# ----- TURBULENCE MODELS ----- #
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure = 4

##########################################
#            Signed Distance             #
##########################################

def signedDistance(x):
   d=abs(x[1] - water_depth  - water_amplitude * cos(x[0]))
   if x[1] < (water_depth + water_amplitude * cos(x[0])):
       return - d
   else: 
       return d
     




 
