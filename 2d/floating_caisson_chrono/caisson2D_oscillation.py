from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
import ChRigidBody as crb
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 0.9, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (2.3, 1.2,), "Dimensions of the tank"),
    ("tank_sponge", (2., 2.), "Length of absorption zones (front/back, left/right)"),
    # waves
    ("waves", False, "Generate waves (True/False)"),
    ("wave_period", 0.8, "Period of the waves"),
    ("wave_height", 0.029, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # caisson
    ("caisson_dim", (0.3, 0.1), "Dimensions of the caisson"),
    ("caisson_coords", (1.15, 0.9), "Dimensions of the caisson"),
    ("caisson_width", 0.9, "Width of the caisson"),
    ("caisson_corner_r", 0., "radius of the corners of the caisson"),
    ("free_x", (0.0, 0.0, 0.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 1.0), "Rotational DOFs"),
    ("VCG", None, "vertical position of the barycenter of the caisson"),
    ("draft", 0.425, "Draft of the caisson"),
    ("caisson_mass", 15., "Mass of the caisson"),
    ("caisson_inertia", 0.236, "Inertia of the caisson"),
    ("rotation_angle", np.pi/12., "Initial rotation angle (in radians)"),
    # mooring
    ("mooring", False, "add moorings"),
    ("mooring_type", 'spring', "type of moorings"),
    ("mooring_anchor", (0.,0.,0.), "anchor coordinates (absolute coorinates)"),
    ("mooring_fairlead", (0.,0.,0.), "fairlead cooridnates (relative coordinates from barycenter)"),
    ("mooring_K", 197.58, "mooring (spring) stiffness"),
    ("mooring_R", 19.8, "mooring (spring) damping"),
    ("mooring_restlength", 0., "mooring (spring) rest length"),
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("refinement_level", 0 ,"Set maximum element diameter to he/2**refinement_level"),
    ("T", 10.0 ,"Simulation time"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.33 ,"Target cfl"),
    ("nsave",  20,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

# waves
if opts.waves is True:
    period = opts.wave_period
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.MonochromaticWaves(period, height, mwl, depth,
                                 np.array([0., -9.81, 0.]), direction)
    wavelength = wave.wavelength
    # tank options
    tank_dim = 8*wavelength
    tank_sponge = (2*wavelength, 2*wavelength)

else:
    tank_dim = opts.tank_dim
    tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
domain2 = Domain.PlanarStraightLineGraphDomain()
# caisson options
dim = opts.caisson_dim
VCG = opts.VCG
if VCG is None:
    VCG = dim[1]/2.
draft = opts.draft
free_x = opts.free_x
free_r = opts.free_r
rotation = opts.rotation_angle
if opts.caisson_coords is None:
    coords = [tank_dim[0]/2., waterLevel]
else:
    coords = opts.caisson_coords
# barycenter = (coords[0], coords[1]-dim[1]/2.+VCG)
# inertia = opts.inertia
width = opts.caisson_width



def quarter_circle(center, radius, p_nb, angle, angle0=0., v_start=0.):
    # p_nb = int(np.ceil(2*np.pi*radius/dx))  # number of points on segment
    # p_nb = refinement
    vertices = []
    segments = []
    for i in range(p_nb):
        x = radius*np.sin(angle0+angle*float(i)/(p_nb-1))
        y = radius*np.cos(angle0+angle*float(i)/(p_nb-1))
        vertices += [[center[0]+x, center[1]+y]]
        if i > 0:
            segments += [[v_start+(i-1), v_start+i]]
        elif i == p_nb-1:
            segments += [[v_start+i, v_start]]
    return vertices, segments

radius = opts.caisson_corner_r
if radius != 0:
    vertices = []
    vertexFlags = []
    segments = []
    segmentFlags = []
    dim = [0.5, 0.5]
    angle0 = [np.pi/2., 0., 3*np.pi/2, np.pi]
    angle1 = [-np.pi/2., -np.pi/2., -np.pi/2., -np.pi/2.]
    centers = [[dim[0]/2.-radius, dim[1]/2.-radius], [-dim[0]/2.+radius, dim[1]/2.-radius],
               [-dim[0]/2.+radius, -dim[1]/2.+radius], [dim[0]/2.-radius, -dim[1]/2.+radius]]
    p_nb = 10
    center = [0., 0.]
    flag = 1
    v_start = 0
    for i in range(len(angle0)):
        v_start = len(vertices)
        v, s = quarter_circle(center=centers[i], radius=radius, p_nb=p_nb,
                                      angle=angle1[i], angle0=angle0[i],
                                      v_start=v_start)
        vertices += v
        vertexFlags += [1]*len(v)
        segments += s+[[len(vertices)-1, len(vertices)]]
        segmentFlags += [1]*len(s)+[1]
    segments[-1][1] = 0  # last segment links to vertex 0
    boundaryTags = {'caisson': 1}
    caisson = st.CustomShape(domain, barycenter=(0.,-0.115, 0.),
                             vertices=vertices, vertexFlags=vertexFlags,
                             segments=segments, segmentFlags=segmentFlags,
                             boundaryTags=boundaryTags)
else:
    caisson = st.Rectangle(domain, dim=opts.caisson_dim)

caisson.setHoles([[0., 0.]])
caisson.holes_ind = np.array([0])
caisson.translate([opts.caisson_coords[0], opts.caisson_coords[1]])
caisson.rotate(np.pi/12., pivot=caisson.barycenter)
# system = crb.System(np.array([0., -9.81, 0.]))
ang = np.pi/12.

system = crb.System(np.array([0., -9.81, 0.]))
body = crb.RigidBody(shape=caisson,
                     system=system,
                     center=caisson.barycenter[:2],
                     rot=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.]),
                     mass = opts.caisson_mass,
                     inertia = np.array([0., 0., opts.caisson_inertia]),
                     free_x = np.array([0., 0., 0.]),
                     free_r = np.array([1., 1., 1.]))

# body.rotation_init=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
body.setRecordValues(all_values=True)
if opts.mooring is True:
    if opts.mooring_type == 'spring':
        body.addSpring(stiffness=opts.mooring_K, damping=opts.mooring_R,
                       fairlead=opts.mooring_barycenter,
                       anchor=opts.mooring_anchor,
                       rest_length=caisson.barycenter[0])



# ----- SHAPES ----- #

tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if left:
    if opts.waves is True:
        tank.setGenerationZones(x_n=left, waves=wave)
        tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave, vert_axis=1)
    else:
        tank.setAbsorptionZones(x_n=left)
if right:
    tank.setAbsorptionZones(x_p=right)
tank.setChildShape(caisson, 0)

# ----- BOUNDARY CONDITIONS ----- #

for bc in caisson.BC_list:
    bc.setNoSlip()

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()





##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]


refinement_level = opts.refinement_level
he = opts.caisson_dim[1]/10.0*(0.5**refinement_level)
domain.MeshOptions.he = he #coarse grid


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

st.assembleDomain(domain)
domain.writeGeo('mesh')
domain.MeshOptions.use_gmsh = True

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=True
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
dt_out =  (T-dt_init)/nDTout
runCFL = opts.cfl

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
useRANS = 0 # 0 -- None
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
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
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

ns_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = max(1.0e-12,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)

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
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))

# tank.BC['y+'].p_dirichlet = twpflowPressure_init
