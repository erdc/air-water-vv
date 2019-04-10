
import numpy as np
from math import cos, ceil, pi
from proteus import (Domain,
                     Context,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
from proteus.Profiling import logEvent
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
#from proteus.mbd import ChRigidBody as crb 				#chrono update
from proteus.mbd import CouplingFSI as cfsi
#from proteus.mbd import pyChronoCore as pych
import pychrono as pych
from proteus.Profiling import logEvent

from proteus import Comm
comm=Comm.init()

#   ____            _            _      ___        _   _
#  / ___|___  _ __ | |_ _____  _| |_   / _ \ _ __ | |_(_) ___  _ __  ___
# | |   / _ \| '_ \| __/ _ \ \/ / __| | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |__| (_) | | | | ||  __/>  <| |_  | |_| | |_) | |_| | (_) | | | \__ \
#  \____\___/|_| |_|\__\___/_/\_\\__|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                           |_|
# Context options
# used in command line directly with option -C
# e.g.: parun [...] -C "g=(0.,-9.81,0.) rho_0=998.2 genMesh=False"
#
# only change/add context options in the "other options section"
# other sections have variables used in _p and _n files

context_options = []
# physical constants
context_options += [
#    ("rho_0", 998.2, "Water density"),
    ("rho_0", 1037.3, "Water (bottom fluid) density"),
	("nu_0", 1.004e-6, "Water kinematic viscosity m/sec^2"),
#    ("rho_1", 1.205, "Air Densiy"),
    ("rho_1", 1014.1, "Air (top fluid) Densiy"),
	("nu_1", 1.5e-5, "Air kinematic viscosity m/sec^2"),
    ("sigma_01", 0., "Surface tension"),
    ("g", (0, -9.81, 0), "Gravitational acceleration vector"),
    ]
# run time options
context_options += [
    ("T", 0.5 ,"Simulation time in s"),
    ("dt_init", 0.001 ,"Value of initial time step"),
    ("dt_fixed", None, "Value of maximum time step"),
    ("archiveAllSteps", False, "archive every steps"),
    ("dt_output", 0.05, "number of saves per second"),
    ("runCFL", 0.5 ,"Target CFL value"),
    ("cfl", 0.5 ,"Target CFL value"),
    ]

context_options += [
    ("water_level", 3.5, "Height of free surface above bottom"),
    ("ecH", 1.5, "how many elements around phi=0"),
    # tank
    ('tank_wavelength_scale', True, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_x', 5., 'Length of tank'),
    ('tank_z', 6., 'Height of tank'),
    ('tank_sponge_lambda_scale', True, 'True: sponge length relative to wavelength'),
    ('tank_sponge_xn', 0., 'length of sponge layer x-'),
    ('tank_sponge_xp', 0., 'length of sponge layer x+'),
    ('tank_sponge_gen', 'x-', 'sponge layers to be gen zones (if waves)'),
    ('tank_sponge_abs', 'x+', 'sponge layers to be abs zones'),
    ('IC', 'AtRest', 'Initial Conditions: Perturbed or AtRest'),
    # chrono options
    ("sampleRate", 0., "sampling rate for chrono. 0 for every timestep"),
    # sphere
    ("fixStructure", False, "fix structure in place"),
    ("free_x", (1., 1., 0.), "Translational DOFs"),
    ("free_r", (0., 0., 1.), "Rotational DOFs"),
    # waves
    ("waves", False, "Generate waves (True/False)"),
    ("wave_period", 1.2, "Period of the waves"),
    ("wave_height",0.10, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # mesh refinement
#    ("he", 0.05, "Set characteristic element size"),
    ("he", 0.5, "Set characteristic element size"),
	# numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "use_gmsh"),
    ("refinement", True, "ref"),
    ("refinement_freesurface", 0.05, "ref"),
    ("refinement_grading", 1.2, "ref"),
    ("movingDomain", True, "True/False"),
    ("addedMass", False, "True/False"),
    ("chrono_dt", 1e-4, "time step in chrono"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("useRANS", 0, "RANS model"),
    ("useSuperlu", True, "RANS model"),
    ]
# instantiate context options
opts=Context.Options(context_options)


rho_0 = opts.rho_0
nu_0 = opts.nu_0
rho_1 = opts.rho_1
nu_1 = opts.nu_1
sigma_01= opts.sigma_01
#g = np.array(opts.g)
pyg = np.array(opts.g)
he = opts.he

# ----- CONTEXT ------ #

# general options
water_level = opts.water_level

# tank options
tank_dim = np.array([opts.tank_x, opts.tank_z])
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


#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

#domain = Domain.PiecewiseLinearComplexDomain()
#domain2 = Domain.PiecewiseLinearComplexDomain()
domain = Domain.PlanarStraightLineGraphDomain()
domain2 = Domain.PlanarStraightLineGraphDomain()

nd=2

# ----- SHAPES ----- #

# TANK
if opts.tank_wavelength_scale and opts.waves:
    tank_dim[0:1] *= wavelength
tank = st.Tank2D(domain, tank_dim)
sponges = {'x-': opts.tank_sponge_xn,
           'x+': opts.tank_sponge_xp}
if opts.tank_sponge_lambda_scale and opts.waves:
    for key in sponges:
        sponges[key] *= wave.wavelength
tank.setSponge(x_n=sponges['x-'], x_p=sponges['x+'])
# to change individual BC functions, example:
# tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 3.*t


# SPHERE
eps=2.0;
coords = np.array([tank_dim[0]/2., tank_dim[1]/2.+eps])
barycenter = np.array([tank_dim[0]/2., tank_dim[1]/2.+eps])
sphere = st.Circle(domain,
                   radius=0.5,
                   coords=coords,
                   barycenter=barycenter,
                   nPoints=int(ceil(2.*pi*tank_dim[0]/opts.he)))
sphere.setHoles([sphere.coords])
# for gmsh:
sphere.holes_ind = np.array([0])
tank.setChildShape(sphere, 0)

domain.MeshOptions.he = 0.5
st.assembleDomain(domain)  # must be called after defining shapes
#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

for key, bc in sphere.BC.items():
    bc.setNoSlip()

for key, bc in tank.BC.items():
    # fix the nodes on the wall of tank
    # in case of moving domain
    bc.setFixedNodes()
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
if opts.waves:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
else:
    tank.BC['x-'].setFreeSlip()
#tank.BC['x+'].setUnsteadyTwoPhaseVelocityInlet(wave = wave, vert_axis = 2, smoothing= 3.0*opts.he)
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()
dragAlpha = 0.5/nu_0
gen = opts.tank_sponge_gen
abso = opts.tank_sponge_abs
if opts.waves is True:
    dragAlpha = 5*(2*np.pi/period)/nu_0
    smoothing = opts.he*3
    tank.setGenerationZones(x_n=('x-' in gen and sponges['x-'] != 0),
                            x_p=('x+' in gen and sponges['x+'] != 0),
                            waves=wave,
                            smoothing=smoothing,
                            dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_n=('x-' in abso and sponges['x-'] != 0),
                        x_p=('x+' in abso and sponges['x+'] != 0),
                        dragAlpha=dragAlpha)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

smoothing = opts.ecH * opts.he
if opts.IC == 'Perturbed':
    k = 0.18
    nd = domain.nd
    def signedDistance(x):
        water_level = (wave.eta(x,0)+wave.mwl)
        phi_z = x[2]-water_level
        return water_level,phi_z

    def eta_IC(x,t):
        return wave.eta(x,0)

    def vel_u(x,t):
        return wave.u(x,0)

    class P_IC:
        #def __init__(self,waterLevel):
        # self.waterLevel = opts.waterLine_z
        def uOfXT(self,x,t):
            if signedDistance(x)[1] < 0:
                H = 0
                dyn = -(eta_IC(x,t)*((cosh(k*(water_level-(water_level-x[2]))))/(cosh(water_level*k)))*rho_0*g[2])
                hyd = -(tank_dim[2] - signedDistance(x)[0])*rho_1*g[2]-(water_level - x[2])*rho_0*g[2]
                pressure = dyn + hyd
                p_air = 0.0
            elif 0 < signedDistance(x)[1] <= smoothing:
                H = smoothedHeaviside(smoothing/2. , signedDistance(x)[1] - smoothing / 2.)
                x_max = x[:]
                x_max[0]=x[0]
                x_max[1]=x[1]
                x_max[2]=x[2] - signedDistance(x)[1]
                dyn = -(eta_IC(x,t)*((cosh(k*(water_level-(-eta_IC(x,t)))))/(cosh(water_level*k)))*rho_0*g[2])
                hyd = -(tank_dim[2] - signedDistance(x_max)[0])*rho_1*g[2]-(water_level - x_max[2])*rho_0*g[2]
                pressure = dyn + hyd
                p_air =  -(tank_dim[2] - signedDistance(x)[0])*rho_1*g[2]
            else:
                H = 1
                p_air =  -(tank_dim[2] - signedDistance(x)[0])*rho_1*g[2]
                pressure = 0.0
            p = H * p_air + (1-H) * pressure
            return p
    class U_IC:
        def uOfXT(self,x,t):
            return tank.BC['x-'].u_dirichlet.uOfXT(x, t)
    class V_IC:
        def uOfXT(self,x,t):
            return tank.BC['x-'].v_dirichlet.uOfXT(x, t)
    class W_IC:
        def uOfXT(self,x,t):
            return tank.BC['x-'].w_dirichlet.uOfXT(x, t)
    class VF_IC:
        def uOfXT(self, x, t):
            return tank.BC['x-'].vof_dirichlet.uOfXT(x, t)
    class PHI_IC:
        def uOfXT(self, x, t):
            return x[nd-1] - signedDistance(x)[0]

elif opts.IC == 'AtRest':
    class P_IC:
        def uOfXT(self, x, t):
            p_L = 0.0
            phi_L = tank_dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            return p_L -pyg[nd-1]*(rho_0*(phi_L - phi)
                                 +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                  -smoothedHeaviside_integral(smoothing,phi)))
    class U_IC:
        def uOfXT(self, x, t):
            return 0.0
    class V_IC:
        def uOfXT(self, x, t):
            return 0.0
    class W_IC:
        def uOfXT(self, x, t):
            return 0.0
    class VF_IC:
        def uOfXT(self, x, t):
            return smoothedHeaviside(smoothing,x[nd-1]-water_level)
    class PHI_IC:
        def uOfXT(self, x, t):
            return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

g = pych.ChVectorD(0.,-9.81,0.0)
system = cfsi.ProtChSystem(sampleRate=opts.sampleRate)
system.ChSystem.Set_G_acc(g)
system.update_substeps = True  # update drag and added mass forces on cable during substeps
#system.chrono_processor = 0
#system.parallel_mode = False
system.setTimeStep(chrono_dt)
system.build_kdtree = True
timestepper = "Euler"
if timestepper == "HHT":
    system.setTimestepperType("HHT")


body = cfsi.ProtChBody(system)
body.attachShape(sphere)
body.setConstraints(free_x=free_x, free_r=free_r)
#mass = 4*np.pi/3*sphere.radius**3*(opts.rho_0 + 1)
mass = 4*np.pi/3*sphere.radius*1040.0
Ixx = Iyy = Izz = 2./5.*mass*sphere.radius**2
body.ChBody.SetMass(mass)
inert = pych.ChVectorD(Ixx, Iyy, Izz)
body.ChBody.SetInertiaXX(inert)
body.setRecordValues(all_values=True)
if opts.fixStructure:
    body.ChBody.SetBodyFixed(True)

#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

outputStepping = TpFlow.OutputStepping(
    final_time=opts.T,
    dt_init=opts.dt_init,
    # cfl=opts.cfl,
    dt_output=opts.dt_output,
    nDTout=None,
    dt_fixed=opts.dt_fixed,
)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
    ns_model=None,
    ls_model=None,
    nd=domain.nd,
    cfl=opts.cfl,
    outputStepping=outputStepping,
    structured=False,
    he=he,
    nnx=None,
    nny=None,
    nnz=None,
    domain=domain,
    initialConditions=initialConditions,
    boundaryConditions=None, # set with SpatialTools,
    useSuperlu=opts.useSuperlu,
)

myTpFlowProblem.movingDomain = opts.movingDomain

params = myTpFlowProblem.Parameters

# line below needed for relaxation zones
# (!) hack
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
myTpFlowProblem.archiveAllSteps = True

# MESH PARAMETERS
params.mesh.genMesh = opts.genMesh
params.mesh.he = opts.he

# PHYSICAL PARAMETERS
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.gravity = np.array(opts.g)
params.physical.surf_tension_coeff = opts.sigma_01

# MODEL PARAMETERS
ind = -1
if opts.movingDomain:
    params.Models.moveMeshElastic.index = ind+1
    ind += 1
params.Models.rans2p.index = ind+1
ind += 1
params.Models.vof.index = ind+1
ind += 1
params.Models.ncls.index = ind+1
ind += 1
params.Models.rdls.index = ind+1
ind += 1
params.Models.mcorr.index = ind+1
ind += 1
if opts.addedMass is True:
    params.Models.addedMass.index = ind+1
    ind += 1

# auxiliary variables
params.Models.rans2p.auxiliaryVariables += [system]
#params.Models.rans2p.weak_bc_penalty_constant = 10./nu_0#Re

if opts.addedMass is True:
    # passed in added_mass_p.py coefficients
    params.Models.addedMass.auxiliaryVariables += [system.ProtChAddedMass]
    max_flag = 0
    max_flag = max(domain.vertexFlags)
    max_flag = max(domain.segmentFlags+[max_flag])
    max_flag = max(domain.facetFlags+[max_flag])
    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
    for s in system.subcomponents:
        if type(s) is cfsi.ProtChBody:
            for i in range(s.i_start, s.i_end):
                flags_rigidbody[i] = 1
    params.Models.addedMass.flags_rigidbody = flags_rigidbody

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|

he = opts.he

mesh_fileprefix = 'mesh'+str(int(he*1000))
domain.MeshOptions.he = he
domain.MeshOptions.setTriangleOptions()
domain.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)

st.assembleDomain(domain)  # must be called after defining shapes

# prescribed_init = False

# if prescribed_init:
#     logEvent('Calculating chrono prescribed motion before starting simulation with dt='+str(time_init_dt)+' for '+str(time_init)+' seconds (this migth take some time)')
#     system.calculate_init()
#     system.setTimeStep(time_init_dt)
#     system.calculate(time_init)  # run chrono before fluid sim for intended time to executed prescribed motion
#     for i in range(int(time_init/1e-3)):
#         system.calculate(1e-3)  # run chrono before fluid sim for intended time to executed prescribed motion
#     logEvent('finished prescribed motion with body at position '+str(body.ChBody.GetPos()))
#     system.setTimeStep(opts.chrono_dt)

if opts.use_gmsh and opts.refinement is True:
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    import py2gmsh
    from py2gmsh import geometry2mesh
    mesh = geometry2mesh(domain)
    grading = np.cbrt(opts.refinement_grading*12/np.sqrt(2))/np.cbrt(1.*12/np.sqrt(2))  # convert change of volume to change of element size
    he = opts.he
    he_max = 100.
    he_max_water = 100.
    field_list = []
    def mesh_grading(start, he, grading):
        return '(abs({start})*({grading}-1)+{he})/{grading}'.format(he=he, start=start, grading=grading)
    def dist_plane(xn, xp, plane='x'):
        x_range = abs(xp-xn)
        dist = '0.5*(abs({plane}-({xn}))+abs({plane}-({xp}))-{x_range})'.format(xn=xn, xp=xp, x_range=x_range, plane=plane)
        return dist

    box = opts.wave_height/2.
    me1 = py2gmsh.Field.MathEval(mesh=mesh)
    dist_z = dist_plane(xn=water_level-box, xp=water_level+box, plane='z')
    #dist = 'sqrt(({dist_x})^2+({dist_y})^2+({dist_z})^2)'.format(dist_x=dist_x, dist_y=dist_y, dist_z=dist_z)
    dist = dist_z
    me1.F = mesh_grading(start=dist, he=he, grading=grading)
    field_list += [me1]

    me3 = py2gmsh.Field.MathEval(mesh=mesh)
    dist = 'Sqrt(({x_center}-x)^2+({y_center}-y)^2+({z_center}-z)^2)-radius'.format(x_center=sphere.coords[0],
                                                                                    y_center=sphere.coords[1],
                                                                                    z_center=sphere.coords[2])
    me3.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me3]

    # background field
    fmin = py2gmsh.Field.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max
    mesh.Coherence = True

    mesh.writeGeo(mesh_fileprefix+'.geo')

# system.calculate_init()
# for i in range(100000):
#     print(i, system.GetChTime(), body.ChBody.GetPos())
#     system.calculate(0.01)
#     print(m1.getTensionBack(), m2.getTensionBack(), m3.getTensionBack())
