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
    ("tank_dim", (3.137*2, 2.5, 3.,), "Dimensions of the tank"),
    ("tank_sponge", (3.137, 3.137*2), "Length of absorption zones (front/back, left/right)"),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", False, "Places Gauges in tank"),
    ("gauge_fixed", False, "Places Gauges in tank"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 1.4185, "Period of the waves"),
    ("wave_height", 0.07, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wave_wavelength", 1., "Direction of the waves (from left boundary)"),
    ("wave_type", 'Fenton', "type of wave"),
    ("eps", 0.5, "eps"),
    ("w", 0., "frequency"),
    ("Bcoeff", np.array([6.9992068e-002,4.8922522e-005,-6.7217371e-007,1.6872617e-008,-1.7685194e-010,1.8800231e-012,0.0,0.0]), "BCoeffs"),
    ("Ycoeff", np.array([6.9865422e-002,2.5069224e-003,1.3384697e-004,8.4828079e-006,5.9092756e-007,4.3715891e-008,3.3715231e-009,2.6808955e-010]), "YCoeffs"),
    # caisson
    ("caisson", True, "caisson"),
    ("caisson_dim", (0.5, 2.48, 0.32), "Dimensions of the caisson"),
    ("caisson_coords", (3.137, 1.25, 1.406), "Dimensions of the caisson"),
    ("caisson_width", 1., "Width of the caisson"),
    ("caisson_corner_r", 0.064, "radius of the corners of the caisson"),
    ("caisson_corner_side", 'bottom', "corners placement"),
    ("caisson_BC", 'freeslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (1.0, 1.0, 1.0), "Translational DOFs"),
    ("free_r", (1.0, 1.0, 1.0), "Rotational DOFs"),
    ("VCG", 0.135, "vertical position of the barycenter of the caisson"),
    ("caisson_mass", 125., "Mass of the caisson"),
    ("caisson_inertia", 4.05, "Inertia of the caisson"),
    ("rotation_angle", 0., "Initial rotation angle (in degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # mooring
    ("mooring", True, "add moorings"),
    ("mooring_type", 'prismatic', "type of moorings"),
    ("mooring_anchor", (2./2.,2.,0.), "anchor coordinates (absolute coorinates)"),
    ("mooring_fairlead", (0.,0.,0.), "fairlead cooridnates (relative coordinates from barycenter)"),
    ("mooring_K", 197.58, "mooring (spring) stiffness"),
    ("mooring_R", 19.8, "mooring (spring) damping"),
    ("mooring_restlength", 0., "mooring (spring) rest length"),
    # mesh refinement
    ("refinement", True, "Gradual refinement"),
    ("he", 0.01, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_caisson", 0.01, "Set maximum characteristic element size on caisson boundary"),
    ("he_max_water", 10, "Set maximum characteristic in water"),
    ("refinement_freesurface", 0.25,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.75,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", 1.2, "Grading of refinement/coarsening (default: 20% volume)"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "True: use Gmsh. False: use Triangle/Tetgen"),
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
    ("parallel", True ,"Run in parallel")])


rho_0 = 998.2
nu_0 = 1.004e-6
rho_1 = 1.205
nu_1 = 1.500e-5
sigma_01 = 0.0
g = np.array([0., 0., -9.81])


# ----- CONTEXT ------ #

wavelength=1.
# general options
waterLevel = opts.water_level
rotation_angle = np.radians(opts.rotation_angle)

# waves
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    if opts.wave_type == 'Linear':
        period = opts.wave_period
        BCoeffs = opts.Bcoeff
        YCoeffs = opts.Ycoeff
    if opts.wave_type == 'Fenton':
        import Fenton
        if opts.w:
            period = 2*np.pi/opts.w
        else:
            period = 2*np.pi/np.sqrt(opts.eps*2*9.81/0.5)
        comm = Comm.get()
        if comm.isMaster():
            Fenton.writeInput(height, depth, period)
            Fenton.runFourier()
        comm.barrier()
        BCoeffs, YCoeffs = Fenton.getBYCoeffs()
        logEvent("BCOEFFS: "+str(BCoeffs))
        logEvent("YCOEFFS: "+str(YCoeffs))
        with open('Solution.res', 'r') as f:
            for line in f:
                if 'Wave period' in line:
                    words = line.split()
                    period = float(words[5])/np.sqrt(9.81/depth)
                if 'Wave length' in line:
                    words = line.split()
                    wavelength = float(words[5])*depth
                    logEvent("PERIOD: "+str(period))
                    logEvent("WAVELENGTH: "+str(wavelength))
        #getFFT.copyFiles()
        wave = wt.MonochromaticWaves(period=period,
                                     waveHeight=height,
                                     mwl=mwl,
                                     depth=depth,
                                     g=g,
                                     waveDir=direction,
                                     wavelength=wavelength,
                                     waveType=opts.wave_type,
                                     Ycoeff=YCoeffs,
                                     Bcoeff=BCoeffs,
                                     Nf=len(BCoeffs),
                                     fast=False)
        wavelength = wave.wavelength


# tank options
if opts.waves is True:
    tank_dim = (2*wavelength, 2.5, opts.water_level*2)
    tank_sponge = (1*wavelength, 2*wavelength)
else:
    tank_dim = opts.tank_dim
    tank_sponge = opts.tank_sponge
logEvent("TANK SPONGE: "+str(tank_sponge))
logEvent("TANK DIM: "+str(tank_dim))




# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()
# caisson options
if opts.caisson is True:
    dim = opts.caisson_dim
    VCG = opts.VCG
    if VCG is None:
        VCG = dim[1]/2.
    free_x = opts.free_x
    free_r = opts.free_r
    rotation = np.radians(opts.rotation_angle)
    if opts.caisson_coords is None:
        coords = [tank_dim[0]/2., opts.water_level]
    else:
        coords = opts.caisson_coords
    barycenter = (0., 0., -dim[1]/2.+VCG)
    width = opts.caisson_width
    inertia = opts.caisson_inertia/width


    caisson_coords = coords
    caisson_dim = opts.caisson_dim

    if opts.he_caisson:
        he_caisson = opts.he_caisson
    else:
        he_caisson = opts.he




    def quarter_circle(center, radius, p_nb, angle, angle0=0., v_start=0.):
        # p_nb = int(np.ceil(2*np.pi*radius/dx))  # number of points on segment
        # p_nb = refinement
        vertices = []
        segments = []
        for i in range(p_nb):
            x = radius*np.sin(angle0+angle*float(i)/(p_nb-1))
            y = radius*np.cos(angle0+angle*float(i)/(p_nb-1))
            vertices += [[center[0]+x, center[1]-dim[1]/2., center[2]+y]]
            if i > 0:
                segments += [[v_start+(i-1), v_start+i]]
            elif i == p_nb-1:
                segments += [[v_start+i, v_start]]
        return vertices, segments

    radius = opts.caisson_corner_r

    nb = int((np.pi*2*radius/4.)/(2*he_caisson))
    #nb = int(np.pi/2/(np.arcsin(he_caisson/2./radius)*2))
    if radius != 0:
        vertices = []
        vertexFlags = []
        segments = []
        segmentFlags = []
        dim = opts.caisson_dim
        if opts.caisson_corner_side == 'bottom':
            angle0 = [np.pi/2., 0., 3*np.pi/2, np.pi]
            angle1 = [-np.pi/2., -np.pi/2., -np.pi/2., -np.pi/2.]
            centers = [[dim[0]/2., 0., dim[2]/2.], [-dim[0]/2., 0., dim[2]/2.],
                    [-dim[0]/2.+radius, 0., -dim[2]/2.+radius], [dim[0]/2.-radius, 0., -dim[2]/2.+radius]]
            p_nb = [0, 0, nb, nb]
        else:
            angle0 = [np.pi/2., 0., 3*np.pi/2, np.pi]
            angle1 = [-np.pi/2., -np.pi/2., -np.pi/2., -np.pi/2.]
            centers = [[dim[0]/2.-radius, 0., dim[2]/2.-radius], [-dim[0]/2.+radius, 0., dim[2]/2.-radius],
                    [-dim[0]/2.+radius, 0., -dim[2]/2.+radius], [dim[0]/2.-radius, 0., -dim[2]/2.+radius]]
            p_nb = [nb, nb, nb, nb]
        center = [0., 0., 0.]
        flag = 1
        v_start = 0
        for i in range(len(angle0)):
            v_start = len(vertices)
            if p_nb[i] != 0:
                v, s = quarter_circle(center=centers[i], radius=radius, p_nb=p_nb[i],
                                            angle=angle1[i], angle0=angle0[i],
                                            v_start=v_start)
            else:
                v = [(np.array(centers[i])-[0., dim[1]/2., 0.]).tolist()]
            vertices += v
        facets = [[[i for i in range(len(vertices))]]]

        vertices_side2 = np.array(vertices)+[0., dim[1], 0.]
        facets += [[[i+len(vertices) for i in range(len(vertices_side2))]]]

        for i in range(len(vertices)-1):
            facets += [[[i, i+1, i+1+len(vertices), i+len(vertices)]]]
        facets += [[[0, len(vertices)-1, 2*len(vertices)-1, len(vertices)]]]
        facetFlags = [1 for i in range(len(facets))]

        vertices += vertices_side2.tolist()
        vertexFlags = [1 for i in range(len(vertices))]

            #segments_side2 = np.array(copy.deepcopy(segments([[len(vertices)-1, len(vertices)]]
        #segments[-1][1] = 0  # last segment links to vertex 0
        boundaryTags = {'caisson': 1}
        caisson = st.CustomShape(domain, barycenter=barycenter,
                                vertices=vertices, vertexFlags=vertexFlags,
                                #segments=segments, segmentFlags=segmentFlags,
                                facets=facets, facetFlags=facetFlags,
                                boundaryTags=boundaryTags)
        caisson.regionFlags = np.array([1])
    else:
        caisson = st.Rectangle(domain, dim=opts.caisson_dim)
    caisson.volumes = [[[i for i in range(len(facets))]]]
    caisson.setHoles([[0., 0., 0.]])
    caisson.holes_ind = np.array([0])
    caisson.translate([caisson_coords[0], caisson_coords[1], caisson_coords[2]])


    ang = rotation_angle
    rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)
    system = crb.ProtChSystem(g)
    system.setTimeStep(opts.chrono_dt)
    body = crb.ProtChBody(shape=caisson, system=system)
    from proteus.mbd import pyChronoCore as pych
    x, y, z = caisson.barycenter
    pos = pych.ChVector(x, y, z)
    e0, e1, e2, e3 = rotation_init
    rot = pych.ChQuaternion(e0, e1, e2, e3)
    inertia = pych.ChVector(1., 1., inertia)
    body.ChBody.SetPos(pos)
    body.ChBody.SetRot(rot)
    body.ChBody.SetMass(opts.caisson_mass)
    body.ChBody.SetInertiaXX(inertia)
    body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))

    # body.setInitialRot(rotation_init)
    # body.rotation_init=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    body.setRecordValues(all_values=True)
    if opts.mooring is True:
        if opts.mooring_type == 'spring':
            body.addSpring(stiffness=opts.mooring_K, damping=opts.mooring_R,
                           fairlead=np.array(opts.mooring_fairlead),
                           anchor=np.array(opts.mooring_anchor),
                           rest_length=opts.mooring_restlength)
        elif opts.mooring_type == 'prismatic':
            body.addPrismaticLinksWithSpring(stiffness=opts.mooring_K, damping=opts.mooring_R,
                           pris1=np.array([0., 0., caisson.barycenter[2]]),
                           pris2=np.array([0.,0.,0.]),
                           rest_length=caisson.barycenter[0])



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
tank = st.Tank3D(domain, tank_dim)
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)

# remove segments
tank.segments = None
tank.segmentFlags = None

# ----- BOUNDARY CONDITIONS ----- #

tank.BC['z+'].setAtmosphere()
if opts.tank_BC == 'noslip':
    tank.BC['z-'].setNoSlip()
    tank.BC['y-'].setNoSlip()
    tank.BC['y+'].setNoSlip()
if opts.tank_BC == 'freeslip':
    tank.BC['z-'].setFreeSlip()
    tank.BC['y-'].setFreeSlip()
    tank.BC['y+'].setFreeSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes
tank.BC['z+'].setTank()  # sliding mesh nodes
tank.BC['z-'].setTank()  #sliding mesh nodes

# absorption/generation zones
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
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
        tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
    else:
        tank.setAbsorptionZones(x_n=left, dragAlpha=dragAlpha)
        tank.BC['x-'].setNoSlip()
else:
    tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right, dragAlpha=dragAlpha)


# ----- GAUGES ----- #

if opts.gauge_output:
    if left or right:
        gauge_dx = tank_sponge[0]/10.
    else:
        gauge_dx = tank_dim[0]/10.
    probes=np.linspace(-tank_sponge[0], tank_dim[0]+tank_sponge[1], (tank_sponge[0]+tank_dim[0]+tank_sponge[1])/gauge_dx+1)
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
    if opts.gauge_fixed:
        PGF = []
        for i in range(4):
            PGF.append((caisson_coords[0]-0.15+0.1*i, waterLevel-0.28, 0.), )
        tank.attachPointGauges(
            'twp',
            gauges = ((('p', 'u', 'v'), PGF),),
            activeTime=(0, opts.T),
            sampleRate=0,
            fileName='pointGauge_fixed.csv'
        )

    #he = opts.caisson_dim[1]/10.0*(0.5**opts.refinement_level)


domain.MeshOptions.he = opts.he
fileprefix = 'meshgeo_'+str(wavelength)+'a'
domain.geofile = fileprefix
domain.MeshOptions.setOutputFiles(name=fileprefix)
#domain.polyfile = fileprefix
st.assembleDomain(domain)
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.use_gmsh = opts.use_gmsh
#geofile='meshgeo_he'+str(opts.he)+'max'+str(opts.he_max)+str(opts.he_max_water)+'c'+str(opts.he_caisson)+'_T'+str(tank_dim)+str(tank_sponge)
#if opts.caisson:
#    geofile += 'C'+str(caisson_dim)+str(caisson_coords)
#if opts.refinement:
#    geofile += 'g'+str(round(opts.refinement_grading, 3))+str(opts.refinement_caisson)+str(opts.refinement_freesurface)

# MESH REFINEMENT

grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
if opts.refinement is True:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    he = opts.he
    he_max = opts.he_max
    he_max_water = opts.he_max_water
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
    box1.ZMin = waterLevel-box
    box1.ZMax = waterLevel+box
    field_list += [box1]

    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)

    math1 = py2gmsh.Fields.MathEval(mesh=mesh)
    math1.F = mesh_grading(start='sqrt((z-{zmax})^2)'.format(zmax=waterLevel+box), he=he, grading=grading)
    field_list += [math1]

    math2 = py2gmsh.Fields.MathEval(mesh=mesh)
    math2.F = mesh_grading(start='sqrt((z-{zmin})^2)'.format(zmin=waterLevel+box), he=he, grading=grading)
    field_list += [math2]

    #cylinder1 = py2gmsh.Fields.Cylinder(mesh=mesh)
    #cylinder1.Radius = 
    #cylinder1.VIn = he
    #cylinder1.VOut = he_max
    #cylinder1.XAxis = 
    #cylinder1.Yaxis = 
    #cylinder1.ZAxis =
    #cylinder1.XCenter = caisson_coords[0]
    #cylinder1.YCenter = caisson_coords[1]
    #cylinder1.ZCenter = caisson_coords[2]


#p0 = py2gmsh.Entity.Point([-tank_sponge[0], waterLevel+box, 0.], mesh=mesh)
#p1 = py2gmsh.Entity.Point([tank_dim[0]+tank_sponge[1], waterLevel+box, 0.], mesh=mesh)
#p2 = py2gmsh.Entity.Point([-tank_sponge[0], waterLevel-box, 0.], mesh=mesh)
#p3 = py2gmsh.Entity.Point([tank_dim[0]+tank_sponge[1], waterLevel-box, 0.], mesh=mesh)
#l1 = py2gmsh.Entity.Line([p0, p1], mesh=mesh)
#l2 = py2gmsh.Entity.Line([p2, p3], mesh=mesh)
#
#bl2 = py2gmsh.Fields.BoundaryLayer(mesh=mesh)
#bl2.hwall_n = he
#bl2.ratio = grading
#bl2.EdgesList = [l1, l2]
#field_list += [bl2]

# max element size in water phase
#box2 = py2gmsh.Fields.Box(mesh=mesh)
#box2.VIn = he_max_water
#box2.VOut = he_max
#box2.XMin = -tank_sponge[0]
#box2.XMax = tank_dim[0]+tank_sponge[1]
#box2.YMin = 0
#box2.YMax = waterLevel
#box2.ZMin = waterLevel+0.1
#box2.ZMax = waterLevel-0.1
#field_list += [box2]

#if opts.caisson:
    # boundary layer on caisson
#bl1 = py2gmsh.Fields.BoundaryLayer()
#bl1.hwall_n = he_caisson
#bl1.ratio = grading
#bl1.EdgesList = mesh.getLinesFromIndex([i+1 for i in range(len(caisson.segments))])
#mesh.addField(bl1)
#field_list += [bl1]

#if opts.refinement_caisson:
#    # create circle (non-physical) around caisson
#    refinement_caisson = opts.refinement_caisson
#    p0 = py2gmsh.Entity.Point([caisson_coords[0], caisson_coords[1], 0.], mesh=mesh)
#    p1 = py2gmsh.Entity.Point([caisson_coords[0]-refinement_caisson, caisson_coords[1], 0.], mesh=mesh)
#    p2 = py2gmsh.Entity.Point([caisson_coords[0]+refinement_caisson, caisson_coords[1], 0.], mesh=mesh)
#    p3 = py2gmsh.Entity.Point([caisson_coords[0]-refinement_caisson+0.00001, caisson_coords[1], 0.], mesh=mesh)
#    c1 = py2gmsh.Entity.Circle(p1, p0, p2, mesh=mesh)
#    c2 = py2gmsh.Entity.Circle(p2, p0, p3, mesh=mesh)
#
#    # refined circle around caisson
#    b1 = py2gmsh.Fields.Ball(mesh=mesh)
#    b1.VIn = he
#    b1.VOut = he_max
#    b1.Radius = refinement_caisson
#    b1.XCenter = caisson_coords[0]
#    b1.YCenter = caisson_coords[1]
#    b1.ZCenter = 0.
#    field_list += [b1]

    # boundary layer on circle around caisson
#    bl3 = py2gmsh.Fields.BoundaryLayer(mesh=mesh)
#    bl3.hwall_n = he
#    bl3.ratio = grading
#    bl3.EdgesList = [c1, c2]
#    field_list += [bl3]

# background field
fmin = py2gmsh.Fields.Min(mesh=mesh)
fmin.FieldsList = field_list
mesh.setBackgroundField(fmin)

# max element size
mesh.Options.Mesh.CharacteristicLengthMax = he_max

mesh.writeGeo(fileprefix+'.geo')





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
