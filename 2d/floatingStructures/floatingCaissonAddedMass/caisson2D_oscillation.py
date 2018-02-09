from proteus import Domain, Context, Comm
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.Profiling import logEvent
#import ChRigidBody as crb
from proteus.mbd import ChRigidBody as crb
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 1.5, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (3., 3.,), "Dimensions of the tank"),
    ("tank_sponge", (0., 0.), "Length of absorption zones (front/back, left/right)"),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", False, "Places Gauges in tank"),
    ("gauge_fixed", False, "Places Gauges in tank"),
    # caisson
    ("addedMass", True, "added mass"),
    ("caisson", True, "caisson"),
    ("caisson_dim", (0.5, 0.2), "Dimensions of the caisson"),
    ("caisson_coords", (1.5, 1.5), "Dimensions of the caisson"),
    ("caisson_width", 1., "Width of the caisson"),
    ("caisson_BC", 'freeslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (1.0, 1.0, 1.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 0.0), "Rotational DOFs"),
    ("VCG", 0.05, "vertical position of the barycenter of the caisson"),
    ("caisson_mass", 25., "Mass of the caisson"),
    ("caisson_inertia", 4.05, "Inertia of the caisson"),
    ("rotation_angle", 0., "Initial rotation angle (in degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # mesh refinement
    ("refinement", True, "Gradual refinement"),
    ("he", 0.03, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_max_water", 10, "Set maximum characteristic in water"),
    ("refinement_freesurface", 0.25,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.75,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("sc", 0.25, "shockCapturing factor"),
    ("weak_factor", 10., "weak bc penalty factor"),
    ("strong_dir", False, "strong dirichlet (True/False)"),
    ])



# ----- CONTEXT ------ #

wavelength=1.
# general options
waterLevel = water_level = opts.water_level
rotation_angle = np.radians(opts.rotation_angle)


# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options
if opts.caisson is True:
    free_x = opts.free_x
    free_r = opts.free_r
    rotation = np.radians(opts.rotation_angle)
    width = opts.caisson_width
    inertia = opts.caisson_inertia/width

    caisson_dim = opts.caisson_dim
    caisson = st.Rectangle(domain, dim=opts.caisson_dim, coords=[0.,0.], barycenter=np.zeros(3))
    ang = rotation_angle
    caisson.setHoles([[0., 0.]])
    caisson.holes_ind = np.array([0])
    print(caisson.regions+np.ones(2))
    print(caisson.regionFlags)

    trans = np.array([opts.caisson_coords[0], opts.caisson_coords[1]])
    print(trans+np.ones(2))
    caisson.translate(trans)
    # system = crb.System(np.array([0., -9.81, 0.]))
    # rotation = np.array([1, 0., 0., 0.])
    rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)
    system = crb.ProtChSystem(np.array([0., -9.81, 0.]))
    system.setTimeStep(opts.chrono_dt)
    system.step_start = 10
    body = crb.ProtChBody(shape=caisson,
                          system=system)
    chbod = body.ChBody
    from proteus.mbd import pyChronoCore as pych
    x, y, z = caisson.barycenter
    pos = pych.ChVector(x, y, z)
    e0, e1, e2, e3 = rotation_init
    rot = pych.ChQuaternion(e0, e1, e2, e3)
    inertia = pych.ChVector(1., 1., inertia)
    chbod.SetPos(pos)
    chbod.SetRot(rot)
    chbod.SetMass(opts.caisson_mass)
    chbod.SetInertiaXX(inertia)
    body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))
    system.setCouplingScheme("CSS", prediction="backwardEuler")

    # body.setInitialRot(rotation_init)
    # body.rotation_init=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    body.setRecordValues(all_values=True)

    for bc in caisson.BC_list:
        if opts.caisson_BC == 'noslip':
            bc.setNoSlip()
        if opts.caisson_BC == 'freeslip':
            bc.setFreeSlip()

    
    def prescribed_motion(t):
        new_x = np.array(caisson_coords)
        new_x[1] = caisson_coords[1]+0.01*cos(2*np.pi*(t/4)+np.pi/2)
        return new_x

    #body.setPrescribedMotion(prescribed_motion)

# ----- SHAPES ----- #
tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)


# ----- BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes


# ----- GAUGES ----- #



domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.use_gmsh = opts.use_gmsh
geofile='mesh'+str(opts.he)
domain.geofile=geofile


# MESH REFINEMENT

if opts.use_gmsh:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = opts.refinement_grading
    he = opts.he
    he_max = opts.he_max
    he_max_water = opts.he_max_water
    ecH = 3.
    if opts.refinement_freesurface > 0:
        box = opts.refinement_freesurface
    else:
        box = ecH*he
    field_list = []

    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)
    box1 = py2gmsh.Fields.Box(mesh=mesh)
    box1.VIn = he
    box1.VOut = he_max
    box1.XMin = -100
    box1.XMax = 100
    box1.YMin = -100
    box1.YMax = 100
    box1.ZMin = water_level-box
    box1.ZMax = water_level+box
    field_list += [box1]

    me01 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist = '(Sqrt(({zcoord}-y)*({zcoord}-y)))'.format(zcoord=water_level+box)
    me01.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me01]
    me02 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist = '(Sqrt(({zcoord}-y)*({zcoord}-y)))'.format(zcoord=water_level-box)
    me02.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me02]

    me3 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist_z = '(abs(abs({z_p}-y)+abs(y-{z_n})-({z_p}-{z_n}))/2.)'.format(z_p=max(caisson.vertices[:,1]), z_n=min(caisson.vertices[:,1]))
    dist_x = '(abs(abs({z_p}-x)+abs(x-{z_n})-({z_p}-{z_n}))/2.)'.format(z_p=max(caisson.vertices[:,0]), z_n=min(caisson.vertices[:,0]))
    #dist_y = 'abs(({y_center}-y)-{radius})/2.)'.format(y_center=caisson.barycenter[1], radius=caisson.radius)
    me3.F = '{he}*{grading}^(Sqrt({dist_x}^2+{dist_z}^2)/{he})'.format(he=he, grading=grading, dist_x=dist_x, dist_z=dist_z)
    me3.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me3]

    field_list += [me02]
    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')



# passed in added_mass_p.py coefficients
max_flag = 0
max_flag = max(domain.vertexFlags)
max_flag = max(domain.segmentFlags+[max_flag])
max_flag = max(domain.facetFlags+[max_flag])
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
for s in system.subcomponents:
    if type(s) is crb.ProtChBody:
        for i in range(s.i_start, s.i_end):
            flags_rigidbody[i] = 1


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
weak_bc_penalty_constant = opts.weak_factor/nu_0#Re
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
if opts.sc == 0.25:
    sc = 0.25 # default: 0.5. Test: 0.25
    sc_beta = 1. # default: 1.5. Test: 1.
    epsFact_consrv_diffusion = 0.1 # default: 1.0. Test: 0.1. Safe: 10.
elif opts.sc == 0.5:
    sc = 0.5
    sc_beta = 1.5
    epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
else:
    import sys
    sys.quit()
ns_forceStrongDirichlet = opts.strong_dir
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

tolfac = 0.001
mesh_tol = 0.001
ns_nl_atol_res = max(1.0e-8,tolfac*he**2)
vof_nl_atol_res = max(1.0e-8,tolfac*he**2)
ls_nl_atol_res = max(1.0e-8,tolfac*he**2)
mcorr_nl_atol_res = max(1.0e-8,0.1*tolfac*he**2)
rd_nl_atol_res = max(1.0e-8,tolfac*he)
kappa_nl_atol_res = max(1.0e-8,tolfac*he**2)
dissipation_nl_atol_res = max(1.0e-8,tolfac*he**2)
mesh_nl_atol_res = max(1.0e-8,mesh_tol*he**2)
am_nl_atol_res = max(1.0e-8,mesh_tol*he**2)

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
