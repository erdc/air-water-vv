from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
#from proteus.mbd import ChRigidBody as crb
#from proteus.mbd import pyChronoCore as pych
import numpy as np

from proteus import Comm
comm=Comm.init()

opts=Context.Options([
    ("water_level", 0.9, "Height of free surface above bottom"),
    # tank
    ("tank_wavelength_scale", True, "if True: tank_x=value*wavelength, tank_y=value*wavelength"),
    ("tank_x", 2., "Length of tank"),
    ("tank_y", 1., "Width of tank"),
    ("tank_z", 1.8, "Height of tank"),
    # chrono options
    ("sampleRate", 0., "sampling rate for chrono. 0 for every timestep"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 1., "Period of the waves"),
    ("wave_height", 0.1, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # mesh refinement
    ("he", 0.05, "Set characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "use_gmsh"),
    ("refinement", True, "ref"),
    ("refinement_freesurface", 0.05, "ref"),
    ("refinement_grading", 1.2, "ref"),
    ("movingDomain", True, "True/False"),
    ("T", 1.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", 1., "Fixed (maximum) time step"),
    ("chrono_dt", 1e-5, "time step in chrono"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4, "Target cfl"),
    ("nsave", 0, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ])


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = np.array([0., 0., -9.81])

# ----- CONTEXT ------ #

# general options
water_level = opts.water_level

# tank options
tank_dim = [opts.tank_x, opts.tank_y, opts.tank_z]

wavelength=1.
# general options

if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = np.array(opts.wave_dir)
    period = opts.wave_period
    BCoeffs = np.zeros(3)
    YCoeffs = np.zeros(3)
    wave = wt.MonochromaticWaves(period=period,
                                 waveHeight=height,
                                 mwl=mwl,
                                 depth=depth,
                                 g=g,
                                 waveDir=direction,
                                 wavelength=wavelength,
                                 waveType='Linear',
                                 Ycoeff=YCoeffs,
                                 Bcoeff=BCoeffs,
                                 Nf=len(BCoeffs),
                                 fast=False)
    wavelength = wave.wavelength

# tank options
if opts.waves is True:
    if opts.tank_wavelength_scale:
        tank_dim = (opts.tank_x*wavelength, opts.tank_y*wavelength, water_level*2)
    tank_sponge = (1*wavelength, 2*wavelength, 1*wavelength, 1*wavelength)
else:
    tank_sponge=(0.,0.,0.,0.)


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

boundaryTags = bt = {'notfixed': 1,
'fixed': 2}
vertices = [[0.,0.], [1.,0.], [1., 1.], [0., 1.]]
vertexFlags = [bt['fixed'], bt['notfixed'], bt['notfixed'], bt['notfixed']]
segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
segmentFlags = [bt['notfixed'], bt['notfixed'], bt['notfixed'], bt['notfixed']]
regions = [[0.5, 0.5]]
regionFlags = [1]

rect = st.CustomShape(domain,
vertices=vertices, vertexFlags=vertexFlags,
segments=segments, segmentFlags=segmentFlags,
regions=regions, regionFlags=regionFlags,
boundaryTags=boundaryTags)

rect.BC['fixed'].setFixedNodes()



# ----- GENERATION / ABSORPTION LAYERS ----- #

class AuxVar():
    def __init__(self):
        self.dt = 0.
        self.dt = 0.
        self.t = 0.
        self.t_last = 0.

    def attachModel(self, model, ar):
        self.model = model
        self.ar = ar
        return self

    def attachAuxiliaryVariables(self,avDict):
        pass

    def calculate_init(self):
        pass

    def calculate(self):
        self.t_last = self.t
        t_last = self.model.stepController.t_model_last
        self.t = t_last


aux = AuxVar()
aux.dt = 0.01

he = opts.he

mesh_fileprefix = 'mesh'+str(int(he*1000))
domain.MeshOptions.he = he
domain.MeshOptions.setTriangleOptions()
domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)

st.assembleDomain(domain)  # must be called after defining shapes

domain.auxiliaryVariables['moveMesh'] = [aux]


if opts.use_gmsh and opts.refinement is True:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he

    mesh.writeGeo(mesh_fileprefix+'.geo')




















##########################################
# Numerical Options and other parameters #
##########################################

from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10./nu_0#Re
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
nd = domain.nd
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
sc = 0.5
sc_beta = 1.5
epsFact_consrv_diffusion = 10.0
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
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
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

ns_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-6,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-6,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-6,opts.mesh_tol*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - water_level
    phi = x[nd-1] - water_level
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))

#isosurface = None
from proteus.Isosurface import Isosurface
#isosurface = Isosurface(isosurfaces=(('phi', (0.,)),), domain=domain, format='h5', sampleRate=opts.sampleRate)
isosurface = None


def load_simulation_globals():
    """Put some variables we need in engine namespace.

    These can then be retrieved by clients for inspection, visualization, etc.
    """
    nodes = isosurface.nodes_array
    triangles = isosurface.elements_array
    x = nodes[:,0]
    y = nodes[:,1]
    z = nodes[:,2]
    vertices = nodes
    nn = len(x)
