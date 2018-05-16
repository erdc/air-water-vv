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
    ("water_level", 0.4, "Height of free surface above bottom"),
    # tank
    ('tank_as_experiment', False, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_wavelength_scale', True, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_x', 2., 'Length of tank'),
    ('tank_y', 0.8, 'Height of tank'),
    ('tank_sponge_wavelength_scale', True, 'True: sponge length relative to wavelength'),
    ("tank_sponge_xn", 1., "Length of absorption zones (front/back, left/right)"),
    ("tank_sponge_xp", 2., "Length of absorption zones (front/back, left/right)"),
    ('tank_sponge_gen', 'x-', 'sponge layers to be gen zones (if waves)'),
    ('tank_sponge_abs', 'x+', 'sponge layers to be abs zones'),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", True, "Places Gauges in tank"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 1., "Period of the waves"),
    ("wave_height", 0.062, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wave_type", 'Fenton', "type of wave"),
    # caisson
    ("addedMass", True, "added mass"),
    ("caisson", True, "caisson"),
    ("caisson_BC", 'freeslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (0.0, 1.0, 0.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 1.0), "Rotational DOFs"),
    ("rotation_angle", 0., "Initial rotation angle (in degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # mesh refinement
    ("refinement", True, "Gradual refinement"),
    ("he", 0.01, "Set characteristic element size"),
    ("refinement_freesurface", 0.062,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("epsFact_density", 1.5, "epsFact_density"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("ns_closure", 0, "ns closure"),
    ("ELLIPTIC_REDISTANCING_TYPE", 0, "Elliptic redistancing type for redist"),
    ])



# ----- CONTEXT ------ #

wavelength=1.
# general options
waterLevel = water_level = opts.water_level

# waves
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    period = opts.wave_period
    BCoeffs = np.zeros(3)
    YCoeffs = np.zeros(3)
    if opts.wave_type in ['Linear', 'Fenton']:
        wave = wt.MonochromaticWaves(period=period,
                                     waveHeight=height,
                                     mwl=mwl,
                                     depth=depth,
                                     g=np.array([0., -9.81, 0.]),
                                     waveDir=direction,
                                     waveType=opts.wave_type,
                                     Nf=len(BCoeffs),
                                     fast=False)
        wavelength = wave.wavelength
    elif opts.wave_type == 'Focused':
        N = 100
        Hs = 0.0194486
        fp = 1.
        xf = 7.
        tf = 20.
        bandFactor = 1.6
        fmax = bandFactor*fp
        fmin = fp/(bandFactor)
        frange = np.linspace(fmin, fmax, N)
        kk = wt.dispersion(2*np.pi*frange, 0.4)
        phi = -kk*xf+(2*np.pi*frange)*tf
        depth = water_level
        wave = wt.RandomWaves(Tp=1./fp,
                              Hs=Hs,
                              mwl=water_level,
                              depth=depth,
                              waveDir=np.array([1.,0.,0.]),
                              g=np.array([0.,-9.81,0.]),
                              N=N,
                              bandFactor=bandFactor,
                              spectName='JONSWAP',
                              spectral_params={'gamma': 3.3,
                                               'TMA': False,
                                               'depth': depth},
                              phi=phi,
                              fast=False)
        wavelength = np.max(2*np.pi/kk)
        period = 1./fp


# tank options
tank_dim = np.array([opts.tank_x, opts.tank_y])
tank_sponge = np.array([opts.tank_sponge_xn, opts.tank_sponge_xp])
if opts.tank_wavelength_scale and opts.waves:
    tank_dim[0] *= wavelength
if opts.tank_as_experiment:
    tank_dim = [18., opts.tank_y]
sponges = {'x-': opts.tank_sponge_xn,
           'x+': opts.tank_sponge_xp}
if opts.tank_sponge_wavelength_scale and opts.waves:
    for key in sponges:
        sponges[key] *= wavelength
logEvent("TANK SPONGE: "+str(tank_sponge))
logEvent("TANK DIM: "+str(tank_dim))

# ----- DOMAIN ----- #

# domain
domain = Domain.PlanarStraightLineGraphDomain()

# caisson options
if opts.caisson is True:
    free_x = opts.free_x
    free_r = opts.free_r
    width = 0.29
    #inertia = 0.34165
    mass = 14.5+0.276
    inertia = (mass-0.276)*(0.1535**2+(0.1-0.0796)**2)
    vertices = np.array([[0.,0.], [0.5,0.], [0.5,0.123], [0.35,0.123], [0.35,0.373],
                         [0.15,0.373], [0.15,0.123], [0., 0.123]])
    vertexFlags = [1 for i in range(len(vertices))]
    segments = np.array([[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,0]])
    segmentFlags = [1 for i in range(len(segments))]
    facets = [[[i for i in range(len(vertices))]]]
    facetFlags = [1]
    regions = np.array([0.25,0.1])
    boundaryTags = {'caisson': 1}
    barycenter = np.array([0.25,0.1,0.])
    caisson = st.CustomShape(domain, barycenter=barycenter,
                            vertices=vertices, vertexFlags=vertexFlags,
                            segments=segments, segmentFlags=segmentFlags,
                            facets=facets, facetFlags=facetFlags,
                            boundaryTags=boundaryTags)
    #caisson = st.Rectangle(domain, dim=(0.5,0.2), barycenter=[0.25,opts.VCG,0.])
    caisson_dim = [0.5, 0.378]
    caisson_coords = [tank_dim[0]/2.-0.25, 0.4-0.1+0.378/2.]
    if opts.tank_as_experiment:
        caisson_coords = [7.-0.25, 0.4-0.1]

    caisson.facetFlags = np.array([1])
    caisson.regionFlags = np.array([1])
    caisson.setHoles([[0.1, 0.1]])
    caisson.holes_ind = np.array([0])
    if opts.wave_type == 'Focused':
        #caisson.translate([7.-sponges['x-'], 0.4-0.1])
        caisson_coords = [7.-0.25, 0.4-0.1]
        caisson.translate([7.-0.25, 0.4-0.1])
    else:
        if opts.tank_as_experiment:
            caisson.translate([7.-0.25, 0.4-0.1])
        else:
            caisson.translate([caisson_coords[0], 0.4-0.1])

    # CHRONO 
    # system
    system = crb.ProtChSystem(np.array([0., -9.81, 0.]))
    system.setTimeStep(opts.chrono_dt)
    system.step_start = 10
    # floating body
    body = crb.ProtChBody(system=system)
    body.attachShape(caisson)
    body.setWidth2D(width)
    chbod = body.ChBody
    from proteus.mbd import pyChronoCore as pych
    x, y, z = caisson.barycenter
    pos = pych.ChVector(x, y, z)
    inertia = pych.ChVector(1., 1., inertia)
    chbod.SetPos(pos)
    chbod.SetMass(mass)
    chbod.SetInertiaXX(inertia)
    body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))
    system.setCouplingScheme("CSS", prediction="backwardEuler")

    body.setRecordValues(all_values=True)

    for bc in caisson.BC_list:
        if opts.caisson_BC == 'noslip':
            bc.setNoSlip()
        if opts.caisson_BC == 'freeslip':
            bc.setFreeSlip()

# tank
tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=sponges['x-'], x_p=sponges['x+'])
left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if opts.waves is True:
    smoothing = opts.he*1.5
    dragAlpha = 5*2*np.pi/period/(1.004e-6)
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
else:
    dragAlpha = 0.5*1./(1.004e-6)
if left:
    if opts.waves is True:
        tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=dragAlpha)
    else:
        tank.setAbsorptionZones(x_n=left, dragAlpha=dragAlpha)
        tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right, dragAlpha=dragAlpha)
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)


# ----- BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes


# ----- GAUGES ----- #

if opts.gauge_output:
    if left or right:
        gauge_dx = tank_sponge[0]/10.
    else:
        gauge_dx = tank_dim[0]/10.
    gauge_dx = 0.01
    probes=np.linspace(-sponges['x-']+0.0001, tank_dim[0]+sponges['x+'], (sponges['x-']+tank_dim[0]+sponges['x+'])/gauge_dx+1)
    PG=[]
    PG2=[]
    LIG = []
    zProbes=waterLevel*0.5
    for i in probes:
        PG.append((i, zProbes, 0.),)
        PG2.append((i, waterLevel, 0.),)
        if i == probes[0]:
            LIG.append(((i, 0.+0.0001, 0.),(i, tank_dim[1]-0.0001,0.)),)
        elif i != probes[0]:
            if opts.caisson:
                if not caisson_coords[0]-caisson_dim[0]*2. < i < caisson_coords[0]+caisson_dim[0]*2.:
                    LIG.append(((i-0.0001, 0.+0.0001, 0.),(i-0.0001, tank_dim[1]-0.0001,0.)),)
            else:
                LIG.append(((i-0.0001, 0.+0.0001, 0.),(i-0.0001, tank_dim[1]-0.0001,0.)),)
    tank.attachPointGauges(
        'twp',
        gauges = ((('p',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_pressure.csv'
    )
    tank.attachPointGauges(
        'ls',
        gauges = ((('phi',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_levelset.csv'
    )

    tank.attachLineIntegralGauges(
        'vof',
        gauges=((('vof',), LIG),),
        activeTime = (0., opts.T),
        sampleRate = 0,
        fileName = 'lineGauge.csv'
    )


domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.use_gmsh = opts.use_gmsh
geofile='mesh'+str(opts.he)+'_'+opts.wave_type
domain.geofile=geofile


# MESH REFINEMENT

if opts.use_gmsh:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = opts.refinement_grading
    he = opts.he
    he_max = 10.
    ecH = 3.
    if opts.refinement_freesurface > 0:
        box = opts.refinement_freesurface
    else:
        box = ecH*he
    field_list = []

    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)

    def dist_plane(xn, xp, plane='x'):
        x_range = abs(xp-xn)
        dist = '0.5*(abs({plane}-({xn}))+abs({plane}-({xp}))-{x_range})'.format(xn=xn, xp=xp, x_range=x_range, plane=plane)
        return dist

    me1 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist_z = dist_plane(xn=water_level-box, xp=water_level+box, plane='y')
    dist_x = dist_plane(xn=-sponges['x-'], xp=tank_dim[0]+sponges['x+'], plane='x')
    dist = 'sqrt(({dist_x})^2+({dist_z})^2)'.format(dist_x=dist_x, dist_z=dist_z)
    #dist = 'sqrt(({dist_z})^2)'.format(dist_z=dist_z)
    me1.F = mesh_grading(start=dist, he=he, grading=grading)
    #me1.F = '{he}*{grading}^({dist}/{he})'.format(dist=dist, he=he, grading=grading)
    field_list += [me1]


    if opts.caisson is True:
        me1 = py2gmsh.Fields.MathEval(mesh=mesh)
        maxv = max(caisson.vertices[:,1])
        minv = min(caisson.vertices[:,1])
        dist_z = dist_plane(xn=minv, xp=maxv, plane='y')
        maxv = max(caisson.vertices[:,0])
        minv = min(caisson.vertices[:,0])
        dist_x = dist_plane(xn=minv, xp=maxv, plane='x')
        dist = 'sqrt(({dist_x})^2+({dist_z})^2)'.format(dist_x=dist_x, dist_z=dist_z)
        #dist = 'sqrt(({dist_z})^2)'.format(dist_z=dist_z)
        me1.F = mesh_grading(start=dist, he=he, grading=grading)
        #me1.F = '{he}*{grading}^({dist}/{he})'.format(dist=dist, he=he, grading=grading)
        field_list += [me1]

        #for s in caisson.segments:
        #    v1 = caisson.vertices[s[0]]
        #    v2 = caisson.vertices[s[1]]
        #    vv = v2-v1
        #    print("VV", vv)
        #    dist = '((({vx})*x-({vy})*y+{v2x}*{v1y}-{v2y}*{v1x})/sqrt(({vx})^2+({vy})^2))'.format(vx=vv[1], vy=vv[0], v1x=v1[0], v1y=v1[1], v2x=v2[0], v2y=v2[1])
        #    me = py2gmsh.Fields.MathEval(mesh=mesh)
        #    #dist_z = dist_plane(xn=water_level-box, xp=water_level+box, plane='y')
        #    #dist_x = dist_plane(xn=-sponges['x-'], xp=tank_dim[0]+sponges['x+'], plane='x')
        #    #dist = 'sqrt(({dist_x})^2+({dist_z})^2)'.format(dist_x=dist_x, dist_z=dist_z)
        #    #dist = 'sqrt(({dist_z})^2)'.format(dist_z=dist_z)
        #    me.F = mesh_grading(start=dist, he=he, grading=grading)
        #    #me1.F = '{he}*{grading}^({dist}/{he})'.format(dist=dist, he=he, grading=grading)
        #    field_list += [me]

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')



if opts.addedMass is True:
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
sc = 0.5
sc_beta = 1.5
epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
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
    epsFact_density    = opts.epsFact_density
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
ns_closure=opts.ns_closure #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

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

if opts.caisson:
    system.calculate_init()
