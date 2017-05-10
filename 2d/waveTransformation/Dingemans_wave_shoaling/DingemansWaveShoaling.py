"""
Dingemans Wave Shoaling
"""
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.86, "Height of free surface above seabed"),
    # tank
    ("tank_dim", (58., 1.26), "Dimensions of the tank"),
    ("tank_sponge", (5., 5.), "Length of relaxation zones zones (left, right)"),
    ("tank_BC", 'freeslip', "Boundary condition at the bottom of the tank (freeslip/noslip)"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 2.02, "Period of the waves"),
    ("wave_height", 0.02, "Height of the waves"),
    ("wave_depth", 0.86, "Wave depth"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("wave_wavelength", 5.037, "Direction of the waves (from left boundary)"), #calculated by FFT
    ("wave_type", 'Fenton', "type of wave"),
    ("Bcoeff", np.array([0.01402408, 0.00008097, 0.00000013, 0.00000000, 0.00000000,
                          0.00000000, 0.00000000, 0.00000000]), "Bcoeffs"),
    ("Ycoeff", np.array([0.01246994, 0.00018698, 0.00000300, 0.00000006, 0.00000000,
                          0.00000000, 0.00000000, 0.00000000]), "Ycoeffs"),
    ("fast", True, "switch for fast cosh calculations in WaveTools"),
    # mesh refinement
    ("refinement", False, "Gradual refinement"),
    ("he", 0.04, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_max_water", 10, "Set maximum characteristic in water phase"),
    ("refinement_freesurface", 0.1,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    #RZ Darcy corrections#  
    ("alpha_value",0.5,"alphaDarcy value"),
    ("fscaling",1,"use fscaling!=0, to switch on frequnecy scaling"),
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
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# waves
omega = 1.
if opts.waves is True:
    period = opts.wave_period
    omega = 2*np.pi/opts.wave_period
    height = opts.wave_height
    mwl = opts.water_level
    depth = opts.wave_depth
    direction = opts.wave_dir
    waves = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                                 g=np.array([0., -9.81, 0.]), waveDir=direction,
                                 wavelength=opts.wave_wavelength,
                                 waveType=opts.wave_type,
                                 Ycoeff=np.array(opts.Ycoeff),
                                 Bcoeff=np.array(opts.Bcoeff),
                                 Nf=len(opts.Bcoeff),
                                 fast=opts.fast)
    wavelength = waves.wavelength

# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

sloped_shore = [[[9.22, 0.],
                 [9.64, 0.06],
                 [15.01, 0.06],
                 [27.04, 0.66],
                 [31.04, 0.66],
                 [37.07, 0.06],
                 [45.39, 0.06],
                 [45.81, 0.]],]

tank = st.TankWithObstacles2D(domain=domain,
                              dim=tank_dim,
                              obstacles=sloped_shore)

tank_sponge = opts.tank_sponge
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
tank.setGenerationZones(x_n=True, waves=waves)
tank.setAbsorptionZones(x_p=True)

he = opts.he
smoothing = he*3.

aa = opts.alpha_value
scalef = omega
if(opts.fscaling ==0):
	scalef=1

# ----- BOUNDARY CONDITIONS ----- #

# Waves
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(waves, smoothing=smoothing, vert_axis=1)

tank.BC['y+'].setAtmosphere()
if opts.tank_BC == 'noslip':
    tank.BC['y-'].setNoSlip()
if opts.tank_BC == 'freeslip':
    tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()


# ----- GAUGES ----- #

gauge_x = [6.26, 10.26, 12.66, 23.26, 27.26, 29.26, 31.26, 33.66, 36.86, 40.26, 44.26]
gauge_y = []
column_gauge_locations = []

for i in range(len(gauge_x)):
    
    if 9.22 < gauge_x[i] < 9.64:
        gauge_y.append( (gauge_x[i]-9.22)*0.06/(9.64-9.22) )
    elif 9.64 <= gauge_x[i] <= 15.01:
        gauge_y.append(0.06)
    elif 15.01 < gauge_x[i] < 27.04:
        gauge_y.append( 0.06+(gauge_x[i]-15.01)*(0.66-0.06)/(27.04-15.01) )
    elif 27.04 <= gauge_x[i] <= 31.04:
        gauge_y.append(0.66)
    elif 31.04 < gauge_x[i] < 37.07:
        gauge_y.append( 0.66+(gauge_x[i]-31.04)*(0.06-0.66)/(37.07-31.04) )
    elif 37.07 <= gauge_x[i] <= 45.39:
        gauge_y.append(0.06)
    elif 45.39 < gauge_x[i] < 45.81:
        gauge_y.append( 0.06+(gauge_x[i]-45.39)*(0.-0.06)/(45.81-45.39) )
    else:
        gauge_y.append(0.)
        
    column_gauge_locations.append(((gauge_x[i], gauge_y[i], 0.), (gauge_x[i], tank_dim[1], 0.)))

tank.attachLineIntegralGauges('vof', gauges=((('vof',),column_gauge_locations),), fileName='column_gauges.csv')

tank.facets = np.array([[[i for i in range(12)]]]+[[[11, 12, 13, 10]]]+[[[8, 14, 15, 9]]])

# ----- ASSEMBLE DOMAIN ----- #

domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.gen_mesh
domain.MeshOptions.he = he
domain.use_gmsh = opts.use_gmsh
st.assembleDomain(domain)

# ----- REFINEMENT OPTIONS ----- #

import py2gmsh 
from MeshRefinement import geometry_to_gmsh 
mesh = geometry_to_gmsh(domain) 

field_list = []
box = 0.1001

box1 = py2gmsh.Fields.Box(mesh=mesh) 
box1.VIn = he/3.
box1.VOut = he 
box1.XMin = -tank_sponge[0] 
box1.XMax = tank_dim[0]+tank_sponge[1] 
box1.YMin = waterLevel-box 
box1.YMax = waterLevel+box 
field_list += [box1]

p0 = py2gmsh.Entity.Point([-tank_sponge[0], waterLevel+box, 0.], mesh=mesh)
p1 = py2gmsh.Entity.Point([tank_dim[0]+tank_sponge[1], waterLevel+box, 0.], mesh=mesh) 
p2 = py2gmsh.Entity.Point([-tank_sponge[0], waterLevel-box, 0.], mesh=mesh) 
p3 = py2gmsh.Entity.Point([tank_dim[0]+tank_sponge[1], waterLevel-box, 0.], mesh=mesh) 
l1 = py2gmsh.Entity.Line([p0, p1], mesh=mesh) 
l2 = py2gmsh.Entity.Line([p2, p3], mesh=mesh)

grading = 1.05
bl2 = py2gmsh.Fields.BoundaryLayer(mesh=mesh) 
bl2.hwall_n = he/3. 
bl2.ratio = grading 
bl2.EdgesList = [l1, l2] 
field_list += [bl2] 


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

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




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
useVF = 1.0
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
mesh.writeGeo(geofile+'.geo') 

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
