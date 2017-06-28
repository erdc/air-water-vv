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
    # cylinder
    ("cylinder", True, "cylinder"),
    ("cylinder_radius", 0.515/2., "radius of cylinder"),
    ("cylinder_height", 0.401, "radius of cylinder"),
    ("cylinder_draft", 0.172, "radius of cylinder"),
    ("cylinder_coords", (1., 1., 0.929), "coordinates of cylinder"),
    ("cylinder_mass", 35.85, "mass of cylinder"),
    ("free_x", (1., 1., 1.), "Translational DOFs"),
    ("free_r", (1., 1., 1.), "Rotational DOFs"),
    ("VCG", 0.0758, "VCG"),
    ("moorings", True, "moorings"),
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
    ("movingDomain", True, "True/False"),
    ("T", 50.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("chrono_dt", 1e-5, "time step in chrono"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4, "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("sc", 0.25, "shockCapturing factor"),
    ])

sampleRate = 0.05

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
cylinder_radius = opts.cylinder_radius
cylinder_height = opts.cylinder_height
free_x = np.array(opts.free_x)
free_r = np.array(opts.free_r)
chrono_dt = opts.chrono_dt

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

print("WAVELENGTH: ", wavelength)

# tank options
if opts.waves is True:
    tank_dim = (2*wavelength, opts.tank_y, water_level*2)
    tank_sponge = (1*wavelength, 2*wavelength, 0., 0.)
else:
    tank_sponge = [0.,0., 0., 0.]


#Palm = True
#if Palm is True:
#    tank_dim = (6., 5., water_level*2.)
#    tank_sponge = (1*wavelength, 2*wavelength, 0., 0.)
    


# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()


# ----- SHAPES ----- #

tank = st.Tank3D(domain, tank_dim)
#tank = st.Cuboid(domain, dim=tank_dim, coords=np.array(tank_dim)/2.)
#tank = st.Tank3D(domain, dim=(5., 5., 1.8))

if opts.cylinder is True:
    nPoints = int(2*np.pi*cylinder_radius/opts.he)
    cylinder_coords = [tank_dim[0]/2., tank_dim[1]/2., water_level+cylinder_height/2.-opts.cylinder_draft]
    barycenter = cylinder_coords-np.array([0.,0.,cylinder_height/2.])+np.array([0.,0.,opts.VCG])
    cylinder = st.Cylinder(domain, radius=cylinder_radius, height=cylinder_height, coords=cylinder_coords, barycenter=barycenter, nPoints=nPoints)
    #cylinder.translate(cylinder_coords)
    cylinder.setHoles([cylinder_coords])
    cylinder.holes_ind = np.array([0])
    cylinder.setBarycenter(barycenter)
    tank.setChildShape(cylinder, 0)



# ----- CHRONO ----- #

g = np.array([0., 0., -9.81])
system = crb.ProtChSystem(g, sampleRate=sampleRate)
#system.chrono_processor = 0
#system.parallel_mode = False
system.setTimeStep(chrono_dt)
system.build_kdtree = True
timestepper = "Euler"
if timestepper == "HHT":
    system.setTimestepperType("HHT")

if opts.cylinder is True:
    body = crb.ProtChBody(system, shape=cylinder)
    body.setConstraints(free_x=free_x, free_r=free_r)
    body.ChBody.SetMass(opts.cylinder_mass)
    body.setRecordValues(all_values=True)
    body.ChBody.SetBodyFixed(False)
    Ixx = Iyy = 0.9
    Izz = cylinder_radius**2/2.*opts.cylinder_mass
    from proteus.mbd import pyChronoCore as pcc
    inert = pcc.ChVector(Ixx, Iyy, Izz)
    body.ChBody.SetInertiaXX(inert)

if opts.moorings is True:
    mesh = crb.Mesh(system)
    mooring_X = 6.660

    anchor = np.array([mooring_X, 0., 0.9])
    dist_from_center = cylinder_radius+0.015
    fairlead1_offset = np.array([dist_from_center*np.cos(np.radians(120)), dist_from_center*np.sin(np.radians(120)), 0.])
    fairlead2_offset = np.array([dist_from_center,0.,0.])
    fairlead3_offset = np.array([dist_from_center*np.cos(np.radians(-120)), dist_from_center*np.sin(np.radians(-120)), 0.])

    fairlead_center = np.array([cylinder_coords[0], cylinder_coords[1], water_level])

    fairlead1 = fairlead_center + fairlead1_offset
    fairlead2 = fairlead_center + fairlead2_offset
    fairlead3 = fairlead_center + fairlead3_offset

    anchor1 = fairlead1-[0.,0.,fairlead1[2]] + np.array([mooring_X*np.cos(np.radians(120)), mooring_X*np.sin(np.radians(120)), 0.])
    anchor2 = fairlead2-[0.,0.,fairlead2[2]] + np.array([mooring_X, 0., 0.])
    anchor3 = fairlead3-[0.,0.,fairlead3[2]] + np.array([mooring_X*np.cos(np.radians(-120)), mooring_X*np.sin(np.radians(-120)), 0.])
        
    from catenary import MooringLine

    L = 6.950
    w = 1.243
    EA = 1e20
    d = 4.786e-3
    A0 = (np.pi*d**2/4)
    nb_elems = 100
    dens = 0.1447/A0
    E = 1e9
    E = 1.6e6/A0
    l1 = MooringLine(L=L, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1, tol=1e-8)
    l2 = MooringLine(L=L, w=w, EA=EA, anchor=anchor2, fairlead=fairlead2, tol=1e-8)
    l3 = MooringLine(L=L, w=w, EA=EA, anchor=anchor3, fairlead=fairlead3, tol=1e-8)

    l1.setVariables()
    l2.setVariables()
    l3.setVariables()

    L = np.array([L])
    d = np.array([d])
    nb_elems = np.array([nb_elems])
    E = np.array([E])
    dens = np.array([dens])
    cable_type = "CableANCF"
    #cable_type = "BeamEuler"
    m1 = crb.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E, beam_type=cable_type)
    m2 = crb.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E, beam_type=cable_type)
    m3 = crb.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E, beam_type=cable_type)

    m1.setNodesPositionFunction(l1.s)
    m2.setNodesPositionFunction(l2.s)
    m3.setNodesPositionFunction(l3.s)

    moorings = [m1, m2, m3]

    body2 = crb.ProtChBody(system) # for anchor
    body2.ChBody.SetBodyFixed(True)

    for m in moorings:
        m.setDragCoefficients(tangential=0.5, normal=2.5, segment_nb=0)
        m.setAddedMassCoefficients(tangential=0., normal=3.8, segment_nb=0)
        if opts.waves is True:
            m.setFluidVelocityFunction(wave.u) # acts only out of domain
        m.external_forces_from_ns = True
        m.external_forces_manual = True
        m.setNodesPosition()
        m.buildNodes()
        m.setFluidDensityAtNodes(np.array([rho_0 for i in range(m.nodes_nb)]))

        if opts.cylinder is True:
            m.attachBackNodeToBody(body)
            m.attachFrontNodeToBody(body)

    # FLOOR
    pos1 = m1.getNodesPosition()
    pos2 = m2.getNodesPosition()
    pos3 = m3.getNodesPosition()

    Efloor = 2e4
    material = pych.ChMaterialSurfaceSMC()
    material.SetYoungModulus(Efloor)
    material.SetFriction(0.3)
    material.SetRestitution(0.2)
    material.SetAdhesion(0)

    m1.setContactMaterial(material)
    m2.setContactMaterial(material)
    m3.setContactMaterial(material)

    m1.setName('mooring1')
    m2.setName('mooring2')
    m3.setName('mooring3')

    vec = pych.ChVector(0., 0., -0.11)
    box_dim = [20.,20.,0.2]
    box = pych.ChBodyEasyBox(box_dim[0], box_dim[1], box_dim[2], 1000, True)
    box.SetPos(vec)
    box.SetMaterialSurface(material)
    box.SetBodyFixed(True)
    system.addBodyEasyBox(box)


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
for key, bc in tank.BC.items():
    # fix the nodes on the wall of tank
    # in case of moving domain
    bc.setFixedNodes()

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




if opts.cylinder is True:
    for key, bc in cylinder.BC.items():
        bc.setNoSlip()
# moving mesh BC were created automatically when
# making a chrono rigid body for the cylinder

domain.MeshOptions.he = opts.he
domain.MeshOptions.setTriangleOptions()
st.assembleDomain(domain)  # must be called after defining shapes

he = opts.he
domain.writePoly('mesh')








tank_sponge = [0,0]




grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
geofile = 'mesh'+str(int(tank_dim[0]*1000))+str(int(tank_sponge[0]*1000))+str(int(he*1000))
if opts.refinement is True:
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

#isosurface = None
from proteus.Isosurface import Isosurface
isosurface = Isosurface(isosurfaces=(('phi', (0.,)),), domain=domain, format='h5', sampleRate=sampleRate)


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
