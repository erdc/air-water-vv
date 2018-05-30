from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.mbd import ChRigidBody as crb
from proteus.mbd import pyChronoCore as pych
import numpy as np
from proteus.Profiling import logEvent

from proteus import Comm
comm=Comm.init()

opts=Context.Options([
    ("water_level", 4., "Height of free surface above bottom"),
    # tank
    ('tank_wavelength_scale', True, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_x', 2., 'Length of tank'),
    ('tank_y', 2., 'Width of tank'),
    ('tank_z', 6., 'Height of tank'),
    ('tank_sponge_lambda_scale', True, 'True: sponge length relative to wavelength'),
    ('tank_sponge_xn', 1., 'length of sponge layer x-'),
    ('tank_sponge_xp', 2., 'length of sponge layer x+'),
    ('tank_sponge_yn', 0., 'length of sponge layer y-'),
    ('tank_sponge_yp', 0., 'length of sponge layer y+'),
    ('tank_sponge_gen', 'x-y-y+', 'sponge layers to be gen zones (if waves)'),
    ('tank_sponge_abs', 'x+', 'sponge layers to be abs zones'),
    ('IC', 'Perturbed', 'Initial Conditions: Perturbed or AtRest'),
    # chrono options
    ("sampleRate", 0., "sampling rate for chrono. 0 for every timestep"),
    # cylinder
    ("scale", 50., "scale used to reduce the dimensions of the structure"),
    ("cylinder", True, "cylinder"),
    ("turbine", True, "takes turbine mass into account (different COG)"),
    ("cylinder_draft", 20., "draft of the cylinder without scaling"),
    ("free_x", (1., 1., 1.), "Translational DOFs"),
    ("free_r", (1., 1., 1.), "Rotational DOFs"),
    # moorings
    ("moorings", True, "moorings"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("waves_case", None, "Takes predefined wave conditions [1 to 7]"),
    ("wave_period", 1.2, "Period of the waves"),
    ("wave_height",0.10, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # mesh refinement
    ("he", 0.1, "Set characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "use_gmsh"),
    ("refinement", True, "ref"),
    ("refinement_freesurface", 0.05, "ref"),
    ("refinement_grading", 1.2, "ref"),
    ("movingDomain", True, "True/False"),
    ("T", 5.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("chrono_dt", 1e-4, "time step in chrono"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4, "Target cfl"),
    ("nsave", 20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ])


wave_heights = [1.92, 7.578, 7.136, 7.574, 10.304, 10.74, 11.122]
wave_periods = [7.5, 12.1, 14.3, 20., 12.1, 14.3, 20.]
scale = opts.scale


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
tank_dim = np.array([opts.tank_x, opts.tank_y, opts.tank_z])
#cylinder_radius = opts.cylinder_radius
#cylinder_height = opts.cylinder_height
free_x = np.array(opts.free_x)
free_r = np.array(opts.free_r)
chrono_dt = opts.chrono_dt

wavelength=1.
# general options

wave_to_initial_conditions = True

if opts.waves is True:
    mwl = depth = opts.water_level
    direction = np.array(opts.wave_dir)
    height = opts.wave_height
    period = opts.wave_period
    if opts.waves_case is not None:
        assert opts.waves_case < 7, 'wave case must be between 0 and 6'
        heigth = wave_heights[opts.waves_case]/scale
        period = wave_periods[opts.waves_case]/np.sqrt(scale)
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


# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()
domain2 = Domain.PiecewiseLinearComplexDomain()


# ----- SHAPES ----- #

# TANK
if opts.tank_wavelength_scale and opts.waves:
    tank_dim[0:2] *= wavelength
tank = st.Tank3D(domain, tank_dim)
sponges = {'x-': opts.tank_sponge_xn,
           'x+': opts.tank_sponge_xp,
           'y-': opts.tank_sponge_yn,
           'y+': opts.tank_sponge_yp}
if opts.tank_sponge_lambda_scale and opts.waves:
    for key in sponges:
        sponges[key] *= wave.wavelength
tank.setSponge(x_n=sponges['x-'], x_p=sponges['x+'], y_n=sponges['y-'], y_p=sponges['y+'])
# to change individual BC functions, example:
# tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 3.*t
for key, bc in tank.BC.items():
    # fix the nodes on the wall of tank
    # in case of moving domain
    bc.setFixedNodes()
tank.BC['z+'].setAtmosphere()
tank.BC['z-'].setFreeSlip()
if opts.waves:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
else:
    tank.BC['x-'].setFreeSlip()
#tank.BC['x+'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
tank.BC['x+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
#tank.BC['y-'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
tank.BC['y+'].setFreeSlip()
#tank.BC['y+'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
tank.BC['wall'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()
dragAlpha = 0.5/nu_0
gen = opts.tank_sponge_gen
abso = opts.tank_sponge_abs
if opts.waves is True:
    dragAlpha = 5*(2*np.pi/period)/nu_0
    smoothing = opts.he*3
    tank.setGenerationZones(x_n=('x-' in gen and sponges['x-'] != 0),
                            x_p=('x+' in gen and sponges['x+'] != 0),
                            y_n=('y-' in gen and sponges['y-'] != 0),
                            y_p=('y+' in gen and sponges['y+'] != 0),
                            waves=wave,
                            smoothing=smoothing,
                            dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_n=('x-' in abso and sponges['x-'] != 0),
                        x_p=('x+' in abso and sponges['x+'] != 0),
                        y_n=('y-' in abso and sponges['y-'] != 0),
                        y_p=('y+' in abso and sponges['y+'] != 0),
                        dragAlpha=dragAlpha)


# CYLINDER
if opts.cylinder is True:
    #nPoints = int(2*np.pi*cylinder_radius/opts.he)
    #cylinder_coords = [tank_dim[0]/2., tank_dim[1]/2., water_level+cylinder_height/2. - opts.cylinder_draft]
    #barycenter = [tank_dim[0]/2., tank_dim[1]/2.,  
    cylinder = st.ShapeSTL(domain, 'semi-sub0_50.stl')
    #cylinder_draft = opts.cylinder_draft
    #scale = scale
    #draft = cylinder_draft/scale
    #cylinder = st.Cylinder(domain, radius=cylinder_radius, height=cylinder_height, coords=cylinder_coords, barycenter=barycenter, nPoints=nPoints )
    location = [tank_dim[0]/2., tank_dim[1]/2., opts.water_level - opts.cylinder_draft/scale]
    mx = (min(cylinder.vertices[:,0])+max(cylinder.vertices[:,0]))/2.
    my = (min(cylinder.vertices[:,1])+max(cylinder.vertices[:,1]))/2.
    mz = (min(cylinder.vertices[:,2])+max(cylinder.vertices[:,2]))/2.
    cylinder.setRegions([[mx, my, mz]], [1])
    cylinder.setHoles([[mx, my, mz]])
    if opts.turbine is True:
        KG = 10.11/scale  # from Maine report
    else:
        KG = 5.6/scale
    cylinder.setBarycenter(np.array([0.5, 0.28876, min(cylinder.vertices[:,2])+KG])) # X and Y coords ASD taken from CFX, scaled and translated - z coords from paper
    cylinder.rotate(-np.pi/6., axis=np.array([0.,0.,10.]))
    cylinder.translate(location-np.array([cylinder.barycenter[0], cylinder.barycenter[1], 0.]))
    cylinder.holes_ind = np.array([0])
    tank.setChildShape(cylinder, 0)
    for key, bc in cylinder.BC.items():
        bc.setNoSlip()


# ----- CHRONO ----- #

g = np.array([0., 0., -9.81])
system = crb.ProtChSystem(gravity=g, sampleRate=opts.sampleRate)
system.update_substeps = True  # update drag and added mass forces on cable during substeps
#system.chrono_processor = 0
#system.parallel_mode = False
system.setTimeStep(chrono_dt)
system.build_kdtree = True
timestepper = "Euler"
if timestepper == "HHT":
    system.setTimestepperType("HHT")


if opts.cylinder is True:
    body = crb.ProtChBody(system)
    body.attachShape(cylinder)
    body.setConstraints(free_x=free_x, free_r=free_r)
    mass = 13444000/scale**3
    Ixx = Iyy = ((23.91+24.90)/scale/2.)**2*mass
    Izz = (32.17/scale)**2*mass
    if opts.turbine is True:
        mass = 14040000/scale**3
        Ixx = Iyy = ((31.61+32.34)/scale/2.)**2*mass
        Izz = (32.17/scale)**2*mass
    body.ChBody.SetMass(mass)
    from proteus.mbd import pyChronoCore as pcc
    inert = pcc.ChVector(Ixx, Iyy, Izz)
    body.ChBody.SetInertiaXX(inert)
    body.setRecordValues(all_values=True)

if opts.moorings is True:
    # variables
    #scale = 50.
    L = 835.5/scale
    d = 0.13376/scale  # equivalent diameter (chain -> cylinder)
    A0 = (np.pi*d**2/4)
    w = 108.63/(scale**3)  # kg/m  # 116.6 vs 108.63
    nb_elems =  50
    dens = w/A0+rho_0
    E = (753.6e6/scale**3)/A0
    fairlead_radius = 40.868/scale
    anchor_radius = 837.6/scale
    fairlead_depth = 14./scale  # only for 20m draft!
    anchor_depth = 200./scale  # only for 20m draft!
    anchor_depth -= 0.01
    mooring_X = anchor_radius-fairlead_radius
    fairlead_height = water_level-fairlead_depth
    anchor_height = water_level-anchor_depth
    # prescribed motion of body
    # start with fully stretched line and go to initial position before simulation starts
    prescribed_init = False
    if prescribed_init:
        fairlead_height = np.sqrt(L**2-mooring_X**2)
        offset_body = fairlead_height-water_level
        coords = body.ChBody.GetPos()+np.array([0.,0.,offset_body])
        vec = pych.ChVector(coords[0], coords[1], coords[2])
        body.ChBody.SetPos(vec)
        time_init = 10.
        time_init_dt = 1e-3
        prescribed_t = np.arange(0., time_init, 1e-1)
        prescribed_z = np.linspace(0., water_level-fairlead_height, len(prescribed_t))
        body.setPrescribedMotionCustom(t=prescribed_t, z=prescribed_z, t_max=time_init)
    # fairleads
    dist_from_center = -fairlead_radius
    fairlead1_offset = np.array([dist_from_center*np.cos(np.radians(120)), dist_from_center*np.sin(np.radians(120)), 0.])
    fairlead2_offset = np.array([dist_from_center,0.,0.])
    fairlead3_offset = np.array([dist_from_center*np.cos(np.radians(-120)), dist_from_center*np.sin(np.radians(-120)), 0.])
    fairlead_center = np.array([cylinder.barycenter[0], cylinder.barycenter[1], fairlead_height])
    fairlead1 = fairlead_center + fairlead1_offset
    fairlead2 = fairlead_center + fairlead2_offset
    fairlead3 = fairlead_center + fairlead3_offset
    # anchors
    anchor1 = fairlead1-[0.,0.,fairlead1[2]] + np.array([-mooring_X*np.cos(np.radians(120)), -mooring_X*np.sin(np.radians(120)), 0.]) + np.array([0.,0.,anchor_height])
    anchor2 = fairlead2-[0.,0.,fairlead2[2]] + np.array([-mooring_X, 0., 0.]) + np.array([0.,0.,anchor_height])
    anchor3 = fairlead3-[0.,0.,fairlead3[2]] + np.array([-mooring_X*np.cos(np.radians(-120)), -mooring_X*np.sin(np.radians(-120)), 0.]) + np.array([0.,0.,anchor_height])
    # quasi-statics for laying out cable
    from catenary import MooringLine
    EA = 1e20  # strong EA so there is no stretching
    l1 = MooringLine(L=L, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1, tol=1e-8)
    l2 = MooringLine(L=L, w=w, EA=EA, anchor=anchor2, fairlead=fairlead2, tol=1e-8)
    l3 = MooringLine(L=L, w=w, EA=EA, anchor=anchor3, fairlead=fairlead3, tol=1e-8)
    if prescribed_init is False:
        l1.setVariables()
        l2.setVariables()
        l3.setVariables()
        # functions to use for laying out the cables
        m1_s = l1.s
        m2_s = l2.s
        m3_s = l3.s
        m1_ds = l1.ds
        m2_ds = l2.ds
        m3_ds = l3.ds
    else:
        m1_s = lambda s: l1.anchor+(l1.fairlead-l1.anchor)*s/L
        m2_s = lambda s: l2.anchor+(l2.fairlead-l2.anchor)*s/L
        m3_s = lambda s: l3.anchor+(l3.fairlead-l3.anchor)*s/L
        m1_ds = lambda s: l1.fairlead-l1.anchor
        m2_ds = lambda s: l2.fairlead-l2.anchor
        m3_ds = lambda s: l3.fairlead-l3.anchor

    # make chrono cables
    mesh = crb.Mesh(system)
    #mesh.setAutomaticGravity(True)
    # make variables arrays
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
    m1.setName('mooring1')
    m2.setName('mooring2')
    m3.setName('mooring3')
    # set functions to lay out nodes
    m1.setNodesPositionFunction(m1_s, m1_ds)
    m2.setNodesPositionFunction(m2_s, m2_ds)
    m3.setNodesPositionFunction(m3_s, m3_ds)
    # for anchor
    body2 = crb.ProtChBody(system)
    body2.ChBody.SetBodyFixed(True)
    body2.barycenter0 = np.zeros(3)
    moorings = [m1, m2, m3]
    for m in moorings:
        #if opts.waves is True:
        #    m.setFluidVelocityFunction(wave.u) # acts only out of domain
        m.external_forces_from_ns = True
        m.external_forces_manual = True
        m.setNodesPosition()
        # set NodesPosition must be calle dbefore buildNodes!
        m.buildNodes()
        m.external_forces_manual = True
        m.external_forces_from_ns = True
        m.setApplyDrag(True)
        m.setApplyBuoyancy(True)
        m.setApplyAddedMass(True)
        m.setFluidDensityAtNodes(np.array([rho_0 for i in range(m.nodes_nb)]))
        m.setDragCoefficients(tangential=1.15, normal=0.213, segment_nb=0)
        m.setAddedMassCoefficients(tangential=0.269, normal=0.865, segment_nb=0)
        # small Iyy for bending
        m.setIyy(1e-20, 0)
        if opts.cylinder is True:
            m.attachBackNodeToBody(body)
        m.attachFrontNodeToBody(body2)
        etas = np.array([-0.5, 0., 0.5, 1.])  # etas to record strain
        m.recordStrainEta(etas)

    # seabed contact
    pos1 = m1.getNodesPosition()
    pos2 = m2.getNodesPosition()
    pos3 = m3.getNodesPosition()
    Efloor = 2e4
    material = pych.ChMaterialSurfaceSMC()
    #material.SetYoungModulus(Efloor)
    material.SetKn(300e6)  # normal stiffness
    material.SetGn(1.)  # normal damping coefficient
    material.SetFriction(0.3)
    material.SetRestitution(0.2)
    material.SetAdhesion(0)
    m1.setContactMaterial(material)
    m2.setContactMaterial(material)
    m3.setContactMaterial(material)
    # build floor
    vec = pych.ChVector(0., 0., -0.1)
    box_dim = [100.,100.,0.2]
    box = pych.ChBodyEasyBox(box_dim[0], box_dim[1], box_dim[2], 1000, True)
    box.SetPos(vec)
    box.SetMaterialSurface(material)
    box.SetBodyFixed(True)
    system.addBodyEasyBox(box)



he = opts.he

mesh_fileprefix = 'mesh'+str(int(he*1000))
domain.MeshOptions.he = he
domain.MeshOptions.setTriangleOptions()
domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)

st.assembleDomain(domain)  # must be called after defining shapes

prescribed_init = False

if prescribed_init:
    logEvent('Calculating chrono prescribed motion before starting simulation with dt='+str(time_init_dt)+' for '+str(time_init)+' seconds (this migth take some time)')
    system.calculate_init()
    system.setTimeStep(time_init_dt)
    system.calculate(time_init)  # run chrono before fluid sim for intended time to executed prescribed motion
    for i in range(int(time_init/1e-3)):
        system.calculate(1e-3)  # run chrono before fluid sim for intended time to executed prescribed motion
    logEvent('finished prescribed motion with body at position '+str(body.ChBody.GetPos()))
    system.setTimeStep(opts.chrono_dt)



if opts.use_gmsh and opts.refinement is True:
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    he = opts.he
    he_max = 100.
    he_max_water = 100.
    field_list = []
    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)
    # refinement free surface
    def dist_plane(xn, xp, plane='x'):
        x_range = abs(xp-xn)
        dist = '0.5*(abs({plane}-({xn}))+abs({plane}-({xp}))-{x_range})'.format(xn=xn, xp=xp, x_range=x_range, plane=plane)
        return dist

    box = opts.wave_height/2.
    me1 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist_z = dist_plane(xn=water_level-box, xp=water_level+box, plane='z')
    #dist = 'sqrt(({dist_x})^2+({dist_y})^2+({dist_z})^2)'.format(dist_x=dist_x, dist_y=dist_y, dist_z=dist_z)
    dist = dist_z
    me1.F = mesh_grading(start=dist, he=he, grading=grading)
    field_list += [me1]

    me3 = py2gmsh.Fields.MathEval(mesh=mesh)
    radius = 1.5
    mz = (min(cylinder.vertices[:,2])+max(cylinder.vertices[:,2]))/2.
    center = cylinder.barycenter-np.array([0.,0.,cylinder.barycenter[2]])+np.array([0.,0.,mz])
    dist_x = '(abs((Sqrt(({x_center}-x)^2+({y_center}-y)^2)-{radius})))'.format(x_center=center[0], radius=radius, y_center=center[1])
    dist_z = '(abs(abs({z_p}-z)+abs(z-{z_n})-({z_p}-{z_n}))/2.)'.format(z_p=max(cylinder.vertices[:,2]), z_n=min(cylinder.vertices[:,2]))
    dist = 'Sqrt(((abs({dist_x}-{radius})+{dist_x})/2.)^2+{dist_z}^2)'.format(dist_x=dist_x, dist_z=dist_z, radius=radius)
    me3.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me3]

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(mesh_fileprefix+'.geo')




system.calculate_init()














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

### initial conditions

waterLine_x = 2 * tank_dim[0]
waterLine_y = 2 * tank_dim[1]
#waterLine_z = wave.eta(x,0.) + wave.mwl


def signedDistance(x):
    water_level = (wave.eta(x,0)+wave.mwl)
    phi_z = x[2]-water_level
    return water_level,phi_z

def vel_u(x,t):
    return wave.u(x,0)

def eta_IC(x,t):
    return wave.eta(x,0)

