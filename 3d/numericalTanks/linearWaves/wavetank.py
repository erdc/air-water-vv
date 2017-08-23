from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.mbd import ChRigidBody as crb
from proteus.mbd import pyChronoCore as pych
import numpy as np

from proteus import Comm
comm=Comm.init()

opts=Context.Options([
    ("water_level", 0.9, "Height of free surface above bottom"),
    # tank
    ("tank_x", 2., "Length of tank"),
    ("tank_y", 2., "Width of tank"),
    ("tank_z", 2., "Height of tank"),
    ("auto_sponge", True, "if True: sets gen/abs according to wavelength"),
    ("sponge_xn", 2., "sponge layer at x-"),
    ("sponge_xp", 2., "sponge layer at x+"),
    ("sponge_yp", 2., "sponge layer at y-"),
    ("sponge_yn", 2., "sponge layer at y+"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 1., "Period of the waves"),
    ("wave_height", 0.04, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # mesh refinement
    ("he", 0.2, "Set characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "use_gmsh"),
    ("refinement", True, "ref"),
    ("refinement_freesurface", 0.1, "ref"),
    ("refinement_grading", 1.2, "ref"),
    ("T", 50.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("cfl", 0.4, "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ])

#sampleRate = 0.05

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = np.array([0., 0., -9.81])

# ----- CONTEXT ------ #

# general options
water_level = opts.water_level

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
                                 waveType='Linear',
                                 Ycoeff=YCoeffs,
                                 Bcoeff=BCoeffs,
                                 Nf=len(BCoeffs),
                                 fast=False)
    wavelength = wave.wavelength

# tank options
tank_dim = [opts.tank_x, opts.tank_y, opts.tank_z]
if opts.auto_sponge is True:
    tank_sponge = [1*wavelength, 2*wavelength, 0., 0.]
else:
    tank_sponge = [opts.tank_sponge_xn, opts.tank_sponge_xp, opts.tank_sponge_yn, opts.tank_sponge_yp]
# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()

# ----- SHAPES ----- #

tank = st.Tank3D(domain, tank_dim)

# ----- BOUNDARY CONDITIONS ----- #

# to change individual BC functions, example:
# tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 3.*t
tank.BC['z+'].setAtmosphere()
tank.BC['z-'].setFreeSlip()
tank.BC['x-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['y+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

if opts.waves is True:
    tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1], y_n=tank_sponge[2], y_p=tank_sponge[3])
    left = right = False
    if tank_sponge[0]: left = True
    if tank_sponge[1]: right = True
    if opts.waves is True:
        dragAlpha = 5*(2*np.pi/period)/nu_0
    else:
        dragAlpha = 0.5/nu_0
    if left:
        if opts.waves is True:
            smoothing = opts.he*3.
            tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=dragAlpha)
            tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=2)
        else:
            tank.setAbsorptionZones(x_n=left, dragAlpha=dragAlpha)
            tank.BC['x-'].setNoSlip()
    else:
        tank.BC['x-'].setNoSlip()
    if right:
        tank.setAbsorptionZones(x_p=right, dragAlpha=dragAlpha)




domain.MeshOptions.he = opts.he
domain.MeshOptions.setTriangleOptions()
st.assembleDomain(domain)  # must be called after defining shapes

he = opts.he
domain.writePoly('mesh')

grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
geofile = 'mesh'+str(int(tank_dim[0]*1000))+str(int(tank_sponge[0]*1000))+str(int(he*1000))
if opts.use_gmsh is True and opts.refinement is True:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    he = opts.he
    he_max = 10.
    he_max_water = 10.
    ecH = 3.
    if opts.refinement_freesurface > 0:
        box = opts.refinement_freesurface
    else:
        box = ecH*he
    field_list = []

    # refinement free surface
    box1 = py2gmsh.Fields.Box(mesh=mesh)
    box1.VIn = he
    box1.VOut = he_max
    box1.XMin = -tank_sponge[0]
    box1.XMax = tank_dim[0]+tank_sponge[1]
    box1.YMin = 0
    box1.YMax = tank_dim[1]
    box1.ZMin = water_level-box
    box1.ZMax = water_level+box
    field_list += [box1]

    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)

    math1 = py2gmsh.Fields.MathEval(mesh=mesh)
    math1.F = mesh_grading(start='sqrt((z-{zmax})^2)'.format(zmax=water_level+box), he=he, grading=grading)
    field_list += [math1]

    math2 = py2gmsh.Fields.MathEval(mesh=mesh)
    math2.F = mesh_grading(start='sqrt((z-{zmin})^2)'.format(zmin=water_level+box), he=he, grading=grading)
    field_list += [math2]

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')

domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.geofile = geofile
domain.MeshOptions.use_gmsh = opts.use_gmsh




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
movingDomain=False
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

ns_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-12,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-12,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-12,opts.mesh_tol*he**2)

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
