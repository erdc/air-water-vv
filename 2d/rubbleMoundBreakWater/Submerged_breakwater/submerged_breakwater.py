"""
Submerged Breakwater
"""
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
import math
import numpy as np

opts=Context.Options([
    # predefined test cases
    ("water_level", 0.315, "Height of free surface above seabed"),
    # tank
    ("h", 0.75, "Height of domain in meters"),
    ("Ls", 2., "Distance of the front toe of the structure end from generation zone in wavelengths"),
    ("Lend", 3.5, "Distance of the back toe of the structure end from absorption zone in meters"),
    ("Lgen", 1., "Length of generation zone in wavelengths"),
    ("Labs", 2., "Length of absorption zone in wavelengths"),
    ("free_slip", True, "Should tank walls have free slip conditions "
                        "(otherwise, no slip conditions will be applied)."),
    # Geometry of a trapezoidal breakwater
    ("b", 0.25, "Width of breakwater at the top"),
    ("hs", 0.25, "Height of the breakwater"),
    ("slope1", 1./2., "Slope1 of the breakwater"),
    ("slope2", 2./3., "Slope2 of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', 0.058, "Mean diameter of the medium"),
    # waves
    ("wave_period", 1., "Period of the waves"),
    ("wave_height", 0.12, "Height of the waves"),
    ("wave_dir", np.array([1.,0.,0.]), "Direction of the waves (from left boundary)"),
    ("wave_wavelength",1.489, "Direction of the waves (from left boundary)"),
    ("wave_type", 'Fenton', "type of wave"),
    ('Ycoeff', np.array([0.23585651, 0.05332571, 0.01505139, 0.00508848,
                         0.00191607, 0.00078766, 0.00037795, 0.00027376]), 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', np.array([0.24551497, 0.01803292, 0.00053150, -0.00000162,
                         0.00001059, 0.00000239, 0.00000023, 0.00000002]), 'Bcoeff only if Fenton is activated'),
    ("fast", True, "switch for fast cosh calculations in WaveTools"),
    # mesh refinement
    ("refinement", False, "Gradual refinement"),
    ("he", 0.04, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_max_water", 10, "Set maximum characteristic in water phase"),
    ("refinement_freesurface", 0.1,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", False, "True/False"),
    ("T", 30.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.5 , "Target cfl"),
    ("nsave",  5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("useVF", 0., "Smoothing function at the water/air interface"),
    ("conservativeFlux", False, "Switches ON/OFF postprocessing on velocity field")])

# ----- Physical constants ----- #

rho_0 = 998.2
nu_0 = 1.004e-6
rho_1 = 1.205
nu_1 = 1.500e-5
sigma_01 = 0.0
g = np.array([0., -9.81, 0.])
gAbs = math.sqrt(sum(g**2))

# ----- CONTEXT ------ #

# waves
period = opts.wave_period
omega = 2*np.pi/opts.wave_period
height = opts.wave_height
mwl = depth = opts.water_level
direction = opts.wave_dir
wave = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                             g=g, waveDir=direction,
                             wavelength=opts.wave_wavelength,
                             waveType=opts.wave_type,
                             Ycoeff=np.array(opts.Ycoeff),
                             Bcoeff=np.array(opts.Bcoeff),
                             Nf=len(opts.Bcoeff),
                             fast=opts.fast)
wavelength = wave.wavelength
waterLevel = opts.water_level

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# refinement
he = opts.he
smoothing = he*3.

# ----- TANK ------ #

hs = opts.hs
b = opts.b
slope1 = opts.slope1
slope2 = opts.slope2

x0 = opts.Ls*wavelength
x1 = x0 + hs/slope1 
x2 = x1 + b
x3 = x2 + hs/slope2

obs = [[[x0, 0.],
        [x1, hs],
        [x2, hs],
        [x3, 0.]],]

obstacle_regions=[[x1+(x2-x1)/2., hs/2.]]

tank_dim = [x3+opts.Lend*wavelength, opts.h]

tank = st.TankWithObstacles2D(domain=domain, dim=tank_dim, obstacles=obs, hole=False,
                              obstacle_regions=obstacle_regions)


# ---- POROUS MEDIA ---- #

porosity = opts.porosity
voidFrac = 1.0-porosity
d50 = opts.d50
d15 = d50/1.2

term1 = 3.12*(10**-3.)
term2 = (gAbs/(nu_0**2.))**(2./3.)
term3 = (d15**2.)
Alpha1 = 1684+term1*term2*term3

term1 = -5.10*(10**-3.)
term2 = (gAbs/(nu_0**2.))**(1./3.)
term3 = (d15)
Beta1 = 1.72+1.57*math.exp(term1*term2*term3)
Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2))
Beta=Beta1*voidFrac/((porosity**3)*d15)

#Proteus scale in viscosity, so i need to divide alpha and beta by nu_0
dragAlpha_P=(porosity**2)*Alpha/nu_0
dragBeta=0.0 #(porosity**3)*Beta/nu_0

# ----- GENERATION / ABSORPTION LAYERS ----- #

tank.setSponge(x_n=opts.Lgen*wavelength, x_p=opts.Labs*wavelength)
dragAlpha = 10.*omega/nu_0
tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing = smoothing)
tank.setAbsorptionZones(x_p=True, dragAlpha = dragAlpha)
tank.setPorousZones(flags=tank.regionFlags[tank.regionIndice['obstacle1']], dragAlpha=dragAlpha_P, dragBeta=dragBeta, porosity=porosity)

# ----- BOUNDARY CONDITIONS ----- #

# Waves
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

# open top
tank.BC['y+'].setAtmosphere()

if opts.free_slip:
    tank.BC['y-'].setFreeSlip()
    tank.BC['x+'].setFreeSlip()

else:  # no slip
    tank.BC['y-'].setNoSlip()
    tank.BC['x+'].setNoSlip()

# sponge
tank.BC['sponge'].setNonMaterial()
tank.BC['obstacle1'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()

# ----- GAUGES ----- #

# Layout 1
probes1 = []
probes1.append(((x0-0.47, 0., 0.), (x0-0.47, tank_dim[1], 0.)),)
probes1.append(((x0-0.35, 0., 0.), (x0-0.35, tank_dim[1], 0.)),)
for i in range(20):
    probes1.append(((x0+1.125+0.85+i*0.1, 0.0, 0.0), (x0+1.125+0.85+i*0.1, tank_dim[1], 0.0)),)

# Layout 2
probes2 = []
probes2.append(((x0+1.125+1.50, 0.0, 0.0),                (x0+1.125+1.50, tank_dim[1], 0.0)),)   
probes2.append(((x0+1.125+1.50+0.56, 0.0, 0.0),           (x0+1.125+1.50+0.56, tank_dim[1], 0.0)),)  
probes2.append(((x0+1.125+1.50+0.56+0.30, 0.0, 0.0),      (x0+1.125+1.50+0.56+0.30, tank_dim[1], 0.0)),)  
probes2.append(((x0+1.125+1.50+0.56+0.30+0.16, 0.0, 0.0), (x0+1.125+1.50+0.56+0.30+0.16, tank_dim[1], 0.0)),)  

# GenZone
genProbes = []
dx = opts.Lgen*wavelength/10.
for j in range(11):
    genProbes.append(((-opts.Lgen*wavelength+j*dx, 0., 0.), (-opts.Lgen*wavelength+j*dx, tank_dim[1], 0.)),)

columnLines1 = tuple(map(tuple,probes1))
columnLines2 = tuple(map(tuple,probes2))
columnLinesG = tuple(map(tuple,genProbes))

tank.attachLineIntegralGauges('vof', gauges=((('vof',),columnLines1),), fileName='line_integral_gauges_1.csv')
tank.attachLineIntegralGauges('vof', gauges=((('vof',),columnLines2),), fileName='line_integral_gauges_2.csv')
tank.attachLineIntegralGauges('vof', gauges=((('vof',),columnLinesG),), fileName='line_integral_gauges_G.csv')

# ----- ASSEMBLE DOMAIN ----- #

tank.facets = np.array([[[0, 1, 2, 3]]]+[[[i for i in range(8)]]]+[[[7, 8, 9, 6]]]+[[[4, 10, 11, 5]]])

domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.gen_mesh
domain.MeshOptions.he = opts.he
domain.use_gmsh = opts.use_gmsh
st.assembleDomain(domain)


# ----- REFINEMENT OPTIONS ----- #

import py2gmsh 
from MeshRefinement import geometry_to_gmsh 
mesh = geometry_to_gmsh(domain)

field_list = []
box_size = 0.07

box = py2gmsh.Fields.Box(mesh=mesh)
box.VIn = 0.02 #he-he/4.
box.VOut = he
box.XMin = -opts.Lgen*wavelength
box.XMax = tank_dim[0]+opts.Labs*wavelength
box.YMin = waterLevel-box_size
box.YMax = waterLevel+box_size
field_list += [box]

p0 = py2gmsh.Entity.Point([-opts.Lgen*wavelength, waterLevel+box_size, 0.], mesh=mesh)
p1 = py2gmsh.Entity.Point([tank_dim[0]+opts.Labs*wavelength, waterLevel+box_size, 0.], mesh=mesh)
p2 = py2gmsh.Entity.Point([-opts.Lgen*wavelength, waterLevel-box_size, 0.], mesh=mesh)
p3 = py2gmsh.Entity.Point([tank_dim[0]+opts.Labs*wavelength, waterLevel-box_size, 0.], mesh=mesh)
l1 = py2gmsh.Entity.Line([p0, p1], mesh=mesh)
l2 = py2gmsh.Entity.Line([p2, p3], mesh=mesh)

grading = 1.05
bl = py2gmsh.Fields.BoundaryLayer(mesh=mesh)
bl.hwall_n = 0.02 #he-he/4.
bl.ratio = grading
bl.EdgesList = [l1, l2]
field_list += [bl]

fmin = py2gmsh.Fields.Min(mesh=mesh)
fmin.FieldsList = field_list
mesh.setBackgroundField(fmin)

mesh.Options.Mesh.CharacteristicLengthMax = he

domain.MeshOptions.genMesh = opts.gen_mesh
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.use_gmsh = opts.use_gmsh

geofile = 'mesh'
mesh.writeGeo(geofile+'.geo')
domain.geofile = geofile


##########################################
# Numerical Options and other parameters #
##########################################

from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = opts.useVF #1.0
useOnlyVF = False
useRANS = opts.useRANS # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters
sc = 0.5 # default: 0.5. Test: 0.25
sc_beta = 1.5 # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 1. # default: 1.0. Test: 0.1
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
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
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = 1e-6#max(1.0e-6,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = 1e-6#max(1.0e-6,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = 1e-6#max(1.0e-6,0.001*domain.MeshOptions.he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
