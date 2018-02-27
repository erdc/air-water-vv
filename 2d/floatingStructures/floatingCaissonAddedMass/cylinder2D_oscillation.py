from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.mbd import ChRigidBody as crb
from math import *
import numpy as np



opts=Context.Options([
    ("addedMass", False, "addedMass"),
    # predefined test cases
    ("water_level", 1.2192, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (1.5, 1.2192*2), "Dimensions of the tank"),
    ("tank_sponge", (1.5, 1.5), "Length of absorption zones (front/back, left/right)"),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", False, "Places Gauges in tank"),
    ("gauge_fixed", False, "Places Gauges in tank"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 0.8, "Period of the waves"),
    ("wave_height", 0.029, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wave_wavelength", 1.596, "Direction of the waves (from left boundary)"),
    ("wave_type", 'Linear', "type of wave"),
    ("Bcoeff", np.zeros(2,), "BCoeffs"),
    ("Ycoeff", np.zeros(2,), "YCoeffs"),
    # caisson
    ("caisson", True, "caisson"),
    ("caisson_dim", 3*0.0254, "Dimensions of the caisson"),
    ("caisson_coords", (0.75, 1.2192-0.0254), "Dimensions of the caisson"),
    ("caisson_width", 1., "Width of the caisson"),
    ("caisson_corner_r", 0.064, "radius of the corners of the caisson"),
    ("caisson_corner_side", 'bottom', "corners placement"),
    ("caisson_BC", 'freeslip', "BC on caisson ('noslip'/'freeslip')"),
    ("free_x", (0, 1., 0), "Translational DOFs"),
    ("free_r", (0, 0, 0), "Rotational DOFs"),
    ("caisson_mass", np.pi*(3*0.0254)**2*998.2/2., "Mass of the caisson"),
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
    ("he", 0.003048, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_max_water", 10, "Set maximum characteristic in water"),
    ("refinement_freesurface", 0.0254*4,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", True, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave",  20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# general options
waterLevel = water_level = opts.water_level
rotation_angle = np.radians(opts.rotation_angle)

# waves
if opts.waves is True:
    period = opts.wave_period
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                                 g=np.array([0., -9.81, 0.]), waveDir=direction,
                                 wavelength=opts.wave_wavelength,
                                 waveType=opts.wave_type,
                                 Ycoeff=np.array(opts.Ycoeff),
                                 Bcoeff=np.array(opts.Bcoeff),
                                 Nf=len(opts.Bcoeff),
                                 fast=False)
    wavelength = wave.wavelength

# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options
if opts.caisson is True:
    dim = opts.caisson_dim
    free_x = opts.free_x
    free_r = opts.free_r
    rotation = np.radians(opts.rotation_angle)
    if opts.caisson_coords is None:
        coords = [tank_dim[0]/2., waterLevel]
    else:
        coords = opts.caisson_coords
    barycenter = (0., 0., 0.)
    width = opts.caisson_width
    inertia = opts.caisson_inertia/width


    caisson_dim = opts.caisson_dim
    caisson_coords = opts.caisson_coords

    def quarter_circle(center, radius, p_nb, angle, angle0=0., v_start=0.):
        # p_nb = int(np.ceil(2*np.pi*radius/dx))  # number of points on segment
        # p_nb = refinement
        vertices = []
        segments = []
        for i in range(p_nb+1):
            x = radius*np.sin(angle0+angle*float(i)/(p_nb))
            y = radius*np.cos(angle0+angle*float(i)/(p_nb))
            vertices += [[center[0]+x, center[1]+y]]
            if i > 0:
                segments += [[v_start+(i-1), v_start+i]]
            elif i == p_nb-1:
                segments += [[v_start+i, v_start]]
        return vertices, segments

    radius = opts.caisson_dim
    if radius != 0:
        p_nb = int((np.pi*2*radius)/(opts.he))
        v, s = quarter_circle(center=[0.,0.], radius=caisson_dim, p_nb=p_nb,
                       angle=2*np.pi, angle0=0., v_start=0.)
        vertices = []
        vertexFlags = []
        segments = []
        segmentFlags = []
        dim = opts.caisson_dim
        center = [0., 0.]
        flag = 1
        v_start = 0
        vertices += v[:-1]
        vertexFlags += [1]*len(vertices)
        segments += s[:-1]+[[len(vertices)-1, 0]]
        segmentFlags += [1]*len(segments)
        segments[-1][1] = 0  # last segment links to vertex 0
        boundaryTags = {'caisson': 1}
        caisson = st.CustomShape(domain, barycenter=barycenter,
                                vertices=vertices, vertexFlags=vertexFlags,
                                segments=segments, segmentFlags=segmentFlags,
                                boundaryTags=boundaryTags)
        facet = []
        for i, vert in enumerate(caisson.vertices):
            facet += [i]
        caisson.facets = np.array([[facet]])
        caisson.facetFlags = np.array([1])
        caisson.regionFlags = np.array([1])
    ang = rotation_angle
    caisson.setHoles([[0., 0.]])
    caisson.holes_ind = np.array([0])
    caisson.translate([caisson_coords[0], caisson_coords[1]])
    # system = crb.System(np.array([0., -9.81, 0.]))
    # rotation = np.array([1, 0., 0., 0.])
    rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)

    use_chrono = True
    if use_chrono:
        system = crb.ProtChSystem(np.array([0., -9.81, 0.]))
        system.setTimeStep(opts.chrono_dt)
        body = crb.ProtChBody(shape=caisson,
                            system=system)
        from proteus.mbd import pyChronoCore as pych
        x, y, z = caisson.barycenter
        pos = pych.ChVector(x, y, z)
        e0, e1, e2, e3 = rotation_init
        rot = pych.ChQuaternion(e0, e1, e2, e3)
        inertia = pych.ChVector(1., 1., inertia)
        # chrono functions
        body.ChBody.SetPos(pos)
        body.ChBody.SetRot(rot)
        body.ChBody.SetMass(opts.caisson_mass)
        body.ChBody.SetInertiaXX(inertia)
        # customised functions
        body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))
        body.setRecordValues(all_values=True)

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
                               pris1=np.array([0.,caisson.barycenter[1],0.]),
                               pris2=np.array([0.,0.,0.]),
                               rest_length=caisson.barycenter[0])
    else:
        A_factor = 1
        from proteus import AuxiliaryVariables
        class BodyAddedMass(AuxiliaryVariables.AV_base):
            def __init__(self):
                AuxiliaryVariables.AV_base.__init__(self)

            def attachModel(self, model, ar):
                self.model=model
                return self
            
            def calculate(self):
                pass
        from proteus.mprans import BodyDynamics as bd
        mass = 25.
        Ixx = [1.,1.,1.]
        body = system = bd.RigidBody(shape=caisson)
        bodyAddedMass = BodyAddedMass()
        body.ProtChAddedMass=body.bodyAddedMass=bodyAddedMass
        system = body  # so it is added in twp_navier_stokes_p.py
        body.setMass(mass)
        body.setConstraints(free_x=free_x,
                            free_r=free_r)
        body.It = np.array([[Ixx[0], 0., 0.],
                            [0., Ixx[1], 0.],
                            [0., 0., Ixx[2]]])
        body.setRecordValues(filename='record_rectangle'+str(A_factor), all_values=True)
        body.coords_system = body.Shape.coords_system  # hack
        def step(dt):
            print("STEPPING")
            nd = body.nd
            body.h[:] = np.zeros(3)
            body.ang_disp[:] = np.zeros(3)
            # calculate displacements (can be changed here) forces: body.F; moments: body.M
            A = body.bodyAddedMass.model.levelModelList[-1].Aij[1].copy()
            A *= A_factor
            acc = (body.F[1])/(body.mass + A[1,1])
            body.acceleration[:] = np.array([0.,acc,0.])
            hz = 0
            velz = body.last_velocity[1]
            substeps = 20
            dt_sub = dt/substeps
            for i in range(substeps):
                velz += acc*dt_sub
                hz += velz*dt_sub
            body.velocity[1] = velz
            body.h[1] = hz
            #body.h[:] = body.getDisplacement(dt)
            #body.ang_disp[:] = body.getAngularDisplacement(dt)
            # translate
            body.Shape.translate(body.h[:nd])
            # rotate
            body.ang = np.linalg.norm(body.ang_disp[:])
            if nd == 2 and body.ang_vel[2] < 0:
                body.ang = -body.ang
            if body.ang != 0.:
                body.Shape.rotate(body.ang, body.ang_vel, body.Shape.barycenter)
                body.rotation[:nd, :nd] = body.Shape.coords_system
                body.rotation_matrix[:] = np.dot(np.linalg.inv(body.last_rotation),
                body.rotation)
                body.rotation_euler[:] = bd.getEulerAngles(body.rotation)
            else:
                body.rotation_matrix[:] = np.eye(3)
                body.barycenter[:] = body.Shape.barycenter
                body.position[:] = body.Shape.barycenter

        body.step = step
        #body.scheme = 'Forward_Euler'



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
left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if left:
    if opts.waves is True:
        smoothing = opts.he*3.
        tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
        tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
    else:
        tank.setAbsorptionZones(x_n=left, dragAlpha=0.5/1.004e-6)
        tank.BC['x-'].setNoSlip()
else: 
    tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right, dragAlpha=0.5/1.004e-6)
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)

# ----- BOUNDARY CONDITIONS ----- #

st_circle = False
if st_circle:
    radius=caisson_dim*3.
    circle_inner1 = st.Circle(domain, radius=radius, coords=caisson_coords, barycenter=barycenter, nPoints=20)
    circle_inner1.setChildShape(caisson, 0)
    tank.setChildShape(circle_inner1, 0)
    circle_inner1.BC['circle'].setNonMaterial()
    #circle.BC['circle'].setFixedNodes()
    circle_inner1.BC['circle'].hx_dirichlet = caisson.BC['caisson'].hx_dirichlet
    circle_inner1.BC['circle'].hy_dirichlet = caisson.BC['caisson'].hy_dirichlet
    circle_inner1.BC['circle'].hz_dirichlet = caisson.BC['caisson'].hz_dirichlet
    circle_inner1.BC['circle'].u_stress = caisson.BC['caisson'].u_stress
    circle_inner1.BC['circle'].v_stress = caisson.BC['caisson'].v_stress
    circle_inner1.BC['circle'].w_stress = caisson.BC['caisson'].w_stress

tank.BC['y+'].setAtmosphere()
if opts.tank_BC == 'noslip':
    tank.BC['y-'].setNoSlip()
if opts.tank_BC == 'freeslip':
    tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()


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



#meshfile='T'+str(tank_dim[0])+str(tank_dim[1])
from proteus import Comm
comm = Comm.get()
comm.barrier()
geofile='mesh_he'+str(opts.he)+'max'+str(opts.he_max)+str(opts.he_max_water)+'_T'+str(tank_dim)+str(tank_sponge)
if opts.caisson:
    geofile += 'C'+str(caisson_dim)+str(caisson_coords)
if opts.refinement:
    geofile += 'g'+str(round(opts.refinement_grading, 3))+str(opts.refinement_caisson)+str(opts.refinement_freesurface)
geofile = geofile.replace(' ', '')
geofile = geofile.replace('(', '')
geofile = geofile.replace(')', '')
geofile = geofile.replace('[', '')
geofile = geofile.replace(']', '')
geofile = geofile.replace('.', '-')
geofile = geofile.replace(',', '_')
domain.MeshOptions.setOutputFiles(name=geofile)
domain.MeshOptions.he = opts.he
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.use_gmsh = opts.use_gmsh

st.assembleDomain(domain)



he = opts.he
if opts.refinement:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = np.sqrt(opts.refinement_grading*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3))
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

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')

#f.write('Point(204) = {{{0}}};\n'.format(str([-tank_sponge[0], waterLevel-box, 0])[1:-1]))
#f.write('Point(205) = {{{0}}};\n'.format(str([tank_dim[0]+tank_sponge[1], waterLevel-box, 0])[1:-1]))
#f.write('Point(206) = {{{0}}};\n'.format(str([-tank_sponge[0], waterLevel+box, 0])[1:-1]))
#f.write('Point(207) = {{{0}}};\n'.format(str([tank_dim[0]+tank_sponge[1], waterLevel+box, 0])[1:-1]))
#f.write('Line(204) = {{{0}}};\n'.format(str([204, 205])[1:-1]))
#f.write('Line(205) = {{{0}}};\n'.format(str([206, 207])[1:-1]))


# passed in added_mass_p.py coefficients
# passed in added_mass_p.py coefficients
max_flag = 0
max_flag = max(domain.vertexFlags)
max_flag = max(domain.segmentFlags+[max_flag])
max_flag = max(domain.facetFlags+[max_flag])
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
flags_rigidbody[1] = 1


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
sc = 0.25 # default: 0.5. Test: 0.25
sc_beta = 1. # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 0.1 # default: 1.0. Test: 0.1
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

tolfac = 0.001
ns_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
am_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mcorr_nl_atol_res = 1e-7 #max(1.0e-12,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he)
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
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
