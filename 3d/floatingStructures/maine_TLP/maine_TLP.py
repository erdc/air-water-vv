import numpy as np
from proteus import Domain, Context, Comm
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt
from proteus.mprans import BoundaryConditions as bc
from proteus.Profiling import logEvent
from proteus import MeshTools as mt
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.mbd import CouplingFSI as fsi
import pychrono as pc


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
   ("rho_0", 998.2, "Water density"),
   ("nu_0", 1.004e-6, "Water kinematic viscosity m/sec^2"),
   ("rho_1", 1.205, "Air Densiy"),
   ("nu_1", 1.5e-5, "Air kinematic viscosity m/sec^2"),
   ("sigma_01", 0., "Surface tension"),
   ("g", np.array([0, -9.81, 0.]), "Gravitational acceleration vector"),
   ("outlet_level", 2.5, "outlet level")

]

# geometry options
context_options +=[
    ("scale", 50, "scale to reduce the structure dimensions by"),
    ("structure_draft", 30., "unscaled draft of the structure"),
    ("free_x", (1., 1., 1.), "Translational DOFs"),
    ("free_r", (1., 1., 1.), "Rotational DOFs"),
    ("Surge_Small_Displace", 0., "Small Displacement in x, full scale"),
    ("Sway_Small_Displace", 0., "Small Displacement in y, full scale"),
    ("Heave_Small_Displace", 0., "Small Displacement in z, full scale"),

]

# moorings
context_options +=[
    ("moorings", True, "moorings"),
    ("mooring_type", "simple_spring", "Catenary or Tendon"),
    ("spring_K", 1718.5, "Spring Stiffness"),
    ("spring_R", 273.6, "Spring Damping"),

]

#wave options
context_options +=[
    ("waveHeight", 0.2, "Wave Height"), 
    ("period", 1.1318, "Wave Period"),
    ("mwl", 4., "Water Level")

]
    
# run time options
context_options += [
    ("T", 1, "Simulation time in s"),
    ("dt_init", 0.0001 ,"Value of initial time step"),
    ("dt_fixed", None, "Value of maximum time step"),
    ("dt_output", 0.1, "number of saves per second"),
    ("cfl", 0.5 ,"Target CFL value"),
    ("chrono_dt", 1e-4, "time step in chrono")

]

#mesh refinement
context_options += [
    ("he", .1, "Set characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file 	    exists"),
    ("use_gmsh", False, "use_gmsh"),
    ("refinement", True, "ref"),
    ("refinement_freesurface", 0.05, "ref"),
    ("refinement_grading", 1.2, "ref"),
    ("movingDomain", False, "True/False")
]

# instantiate context options
opts=Context.Options(context_options)

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PiecewiseLinearComplexDomain()

# Wave
g = np.array([0., -9.81, 0.])
nu_0 = 1.004e-6
period=opts.period

wave = wt.MonochromaticWaves(period=period,
                             waveHeight=opts.waveHeight,
                             mwl=opts.mwl,
                             depth=opts.mwl,
                             g=g,
                             waveDir=(1,0,0),
                             waveType='Linear',
                             fast=True)

#ShapeSTL
#Default size of the domain is a 4x6x2 tank, with a 2x6x2 generation zone, and a 4x6x2 absorbtion zone. (Note that y is the vertical axis)

SG=st.ShapeSTL(domain,'Domain.stl')

boundaryTags= SG.boundaryTags

SG.regions=np.array([[1., 3., 1.], [4., 3., 1.], [8., 3., 1.]])
SG.regionFlags=np.array([1, 2, 3])

dragAlpha = 5*(2*np.pi/period)/nu_0
smoothing = 3. * opts.he
nd = domain.nd
xTop = max(SG.vertices[:,1])

#Structure
structure = st.ShapeSTL(domain, 'CapstoneSTL.stl')
scale=opts.scale
location = [4., opts.mwl-opts.structure_draft/opts.scale, 1.]
mx = 0.
my = (min(structure.vertices[:,2])+max(structure.vertices[:,2]))/2.
mz = 0.

structure.setRegions([[mx, my, mz]], [1])
structure.setHoles([[mx, my, mz]])
KG = 64.06/scale  # from Maine report
structure.setBarycenter(np.array([0.0, min(structure.vertices[:,2])+KG, 0.]))
structure.setRegions([location], [1])
structure.setHoles([location])
structure.holes_ind = np.array([0])

SG.setChildShape(structure, 0)

for key, bc in structure.BC.items():
    bc.setNoSlip()

# Chrono

g = np.array([0., -9.81, 0.])
system = fsi.ProtChSystem()
system.ChSystem.Set_G_acc(pc.ChVectorD(g[0], g[1], g[2]))
system.setTimeStep(opts.chrono_dt)
system.build_kdtree = True
#system.setTimestepperType("Euler".encode('utf-8'))   #Euler is default
body = fsi.ProtChBody(system)
body.attachShape(structure)
body.setConstraints(free_x=np.asarray(opts.free_x), free_r=np.asarray(opts.free_r))
mass = 1361000/scale**3                 #Floater+TowerTop+Tower
#mass = 24.76                             #At Water Density
Ixx = Iyy = ((52.61+52.69)/scale/2.)**2*mass
Izz = (9.4/scale)**2*mass               #Assumed Turbine makes no difference
body.ChBody.SetMass(mass)

inert = pc.ChVectorD(Ixx, Iyy, Izz)
body.ChBody.SetInertiaXX(inert)
body.setRecordValues(all_values=True)


##MOORINGS##


##Relaxation Zones
SG.setGenerationZones(flags=1,
                      epsFact_solid=2.,
                      center=np.array([1., 2.5, 0.5]),
                      orientation=np.array([1., 0., 0.]),
                      waves=wave,
                      dragAlpha=dragAlpha,
                      vert_axis=1,
                      smoothing=smoothing,
                      )

SG.setAbsorptionZones(flags=3,
                      epsFact_solid=4.,
                      center=np.array([8., 2.5, 0.5]),
                      orientation=np.array([-1., 0., 0.]),
                      dragAlpha=dragAlpha,
                      vert_axis=1
                      )

# Boundary Conditions

SG.BC['Top_Gen'].setAtmosphere(orientation=np.array([0,+1,0]))
SG.BC['Bed_Gen'].setFreeSlip()
SG.BC['Wall_Gen'].setFreeSlip()
SG.BC['Top_Tank'].setAtmosphere(orientation=np.array([0,+1,0]))
SG.BC['Bed_Tank'].setFreeSlip()
SG.BC['Wall_Tank'].setFreeSlip()
SG.BC['Top_Abs'].setAtmosphere(orientation=np.array([0,+1,0]))
SG.BC['Bed_Abs'].setFreeSlip()
SG.BC['Wall_Abs'].setFreeSlip()

SG.BC['Inlet_Tank'].setNonMaterial()
SG.BC['Outlet_Tank'].setNonMaterial()

SG.BC['Inlet_Gen'].setUnsteadyTwoPhaseVelocityInlet(wave=wave, 
                                                    vert_axis=1,
                                                    smoothing=3.*opts.he,
                                                    orientation=np.array([-1.,0.,0.]),
                                                    vof_air=1.,
                                                    vof_water=0.)

SG.BC['Outlet_Abs'].setHydrostaticPressureOutletWithDepth(seaLevel= opts.mwl,
                                                          rhoUp=opts.rho_1,
                                                          rhoDown=opts.rho_0,
                                                          g=opts.g,
                                                          refLevel=xTop,
                                                          smoothing=smoothing,
                                                          orientation=np.array([1.,0.,0.]),
                                                          vert_axis=1)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#

vert_axis = np.where(abs(opts.g)>0)[0][0]

class At_Rest:
    def uOfXT(self, x, t):
        return 0.0
class PHI_IC:
    def uOfXT(self, x, t):
        return x[vert_axis]-opts.outlet_level
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(smoothing,opts.outlet_level-x[vert_axis]) 
class P_IC:
    def uOfXT(self, x, t):        
        p_top = 0.0
        phi_top = xTop
        phi = x[vert_axis] - opts.outlet_level
        return p_top - opts.g[vert_axis] * (opts.rho_0 * (phi_top - phi) +
                                       (opts.rho_1 - opts.rho_0) *
                                       (smoothedHeaviside_integral(smoothing, phi_top)-
                                        smoothedHeaviside_integral(smoothing, phi)))

initialConditions = {'pressure': P_IC(),
                     'vel_u': At_Rest(),
                     'vel_v': At_Rest(),
                     'vel_w':At_Rest(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}
outputStepping = TpFlow.OutputStepping(final_time=opts.T,
                                       dt_init=opts.dt_init,
                                       dt_output=opts.dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=opts.he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None)

params = myTpFlowProblem.Parameters
myTpFlowProblem.useSuperLu=False#True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma_01

m = params.Models

m.rdls.p.coefficients.epsFact=0.75
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4
m.rdls.n.maxLineSearches=0
m.rdls.n.maxNonlinearIts=50
#pdb.set_trace()


#Assemble domain
domain.MeshOptions.he = opts.he
mesh_fileprefix = 'mesh_he_'+str(int(opts.he*1000))
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)
st.assembleDomain(domain)
#pdb.set_trace()
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']

