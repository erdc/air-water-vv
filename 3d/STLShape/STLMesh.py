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
    ("U", np.array([0.1,0,0]), "inlet velocity"),
    ("rampTime", 0.1, "inlet ramp"),
   ("outlet_level", 2.1, "outlet level")
    

]
    
# run time options
context_options += [
    ("T", 0.1 ,"Simulation time in s"),
    ("dt_init", 0.0001 ,"Value of initial time step"),
    ("dt_fixed", None, "Value of maximum time step"),
    ("dt_output", 0.1, "number of saves per second"),
    ("cfl", 0.5 ,"Target CFL value")
    ]

context_options += [
   
	# mesh refinement
    ("he", 1., "Set characteristic element size"),
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

SG=st.ShapeSTL(domain,'Blocks.stl')
boundaryTags = SG.boundaryTags


current=wt.SteadyCurrent(U=opts.U,
                         mwl=opts.outlet_level,
                         rampTime=opts.rampTime)

smoothing = 3. * opts.he
nd = domain.nd
xTop = max(SG.vertices[:,1])



SG.BC['Top0'].setAtmosphere(orientation=np.array([0,+1,0]))
SG.BC['Wall0'].setFreeSlip()
SG.BC['Bed0'].setFreeSlip()
SG.BC['Concrete0'].setFreeSlip()



SG.BC['Inlet0'].setUnsteadyTwoPhaseVelocityInlet(wave=current, 
                                                 vert_axis=1,
					         smoothing=3.*opts.he,
						 orientation=np.array([-1.,0.,0.]),
						 vof_air=1.,
                                                 vof_water=0.)

SG.BC['Outlet0'].setHydrostaticPressureOutletWithDepth(seaLevel= opts.outlet_level,
                                                       rhoUp=opts.rho_1,
                                                       rhoDown=opts.rho_0,
                                                       g=opts.g,
                                                       refLevel=xTop,
                                                       smoothing=smoothing,
						       orientation=np.array([1.,0.,0.]),
					               vert_axis=1)

#SG.BC['Inlet0'].setFreeSlip()




#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

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

m.rdls.p.CoefficientsOptions.epsFact=0.75
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4
m.rdls.n.maxLineSearches=0
m.rdls.n.maxNonlinearIts=50

#Assemble domain
domain.MeshOptions.he = opts.he
st.assembleDomain(domain)
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
