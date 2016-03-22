from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.315, "Height of free surface above bottom"), # Choose waterLevels=[0.315, 0.276]
    # initial tank
    #("tank_dim", [10.,10.], "Initial tank dimensions"),
    # waves
    ("wave_period", 2., "Period of the waves"), # Choose periods=[1., 1.5, 2.]
    ("wave_height", 0.03, "Height of the waves"), # Choose waveHeights=[0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12]
    # breakwater
    ("b", 0.25, "Top base of the breakwater"), # Choose b=[0.25, 0.80]
    ("hs", 0.25, "Height of the breakwater"),
    ("slope1", 1./2., "Slope1 of the breakwater"),
    ("slope2", 1./1.5, "Slope2 of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', 0.058, "Mean diameter of the medium"),
    # numerical options
    ("refinement_level", 50 ,"he=walength/refinement_level"),
    ("cfl", 0.9 ,"Target cfl"),
    ("freezeLevelSet", False, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing"),
    ])


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
period=opts.wave_period

waterLevel=opts.water_level

waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface

waveHeight=opts.wave_height

inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])


# ----- Phisical constants ----- #

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))


# ----- WAVE input ----- #

waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=None, # if wave is linear I can use None
                                  waveType="Linear")

#---------Domain Dimension
nd = 2


# ----- SHAPES ----- #

L_leftSpo  = waveinput.wavelength
L_rightSpo = waveinput.wavelength
hs=opts.hs
b=opts.b
slope1=opts.slope1
slope2=opts.slope2

x1=L_leftSpo
x2=x1+2.*waveinput.wavelength
x3=x2+(hs/slope1)
x4=x3+b
x5=x4+(hs/slope2)
x6=x5+waveinput.wavelength
x7=x6+L_rightSpo

#tank_dim=opts.tank_dim
tank_dim = [x7, 0.75]

boundaryOrientations = {'bottom': [0., -1.,0.],
                        'right': [1., 0.,0.],
                        'top': [0., 1.,0.],
                        'left': [-1., 0.,0.],
                        'sponge': None,
                       }
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,
                'sponge': 5,
               }


vertices=[[0.0, 0.0],#0
          [x1,  0.0],#1
          [x2,  0.0],#2
          [x3,  hs ],#3
          [x4,  hs ],#4
          [x5,  0.0],#5
          [x6,  0.0],#6
          [x7,  0.0],#7
          [x7,    tank_dim[1]],#8
          [x6,    tank_dim[1]],#9
          [x1,    tank_dim[1]],#10
          [0.0,   tank_dim[1]],#11
          ]

vertexFlags=np.array([1, 1, 1, 1, 1, 1, 1, 1,
                      3, 3, 3, 3,
                     ])


segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],

          [7,8],

          [8,9],
          [9,10],
          [10,11],

          [11,0],

          [1,10],
          [6,9],
          [2,5],
         ]

segmentFlags=np.array([1, 1, 5, 5, 5, 1, 1,
                       2,
                       3, 3, 3,
                       4,
                       5, 5, 1,
                      ])


regions = [ [ 0.95*x1 , 0.10*tank_dim[1] ],
            [ 0.95*x2 , 0.10*tank_dim[1] ],
            [ x3+0.05 , 0.95*hs ],
            [ 0.95*x7 , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3, 4])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


##########################################
# POROUS MEDIA
##########################################

porosity=opts.porosity
voidFrac=1.0-porosity
d50=opts.d50
d15=d50/1.2

term1=3.12*(10**-3.)
term2=(gAbs/(nu_0**2.))**(2./3.)
term3=(d15**2.)
Alpha1=1684+term1*term2*term3 #Shih
#Alpha1=150 #Ergun
#Alpha1=360 #Engelund

term1=-5.10*(10**-3.)
term2=(gAbs/(nu_0**2.))**(1./3.)
term3=(d15)
Beta1=1.72+1.57*exp(term1*term2*term3) #Shih
#Beta1=1.75 #Ergun
#Beta1=3.6 #Engelund

#Alpha=Alpha1*nu_0*(voidFrac**3)/((porosity**2)*(d15**2))  #Engelund
Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2))  #Ergun
Beta=Beta1*voidFrac/((porosity**3)*d15)

#Proteus scale in viscosity, so i need to divide alpha and beta by nu_0
dragAlpha=porosity*Alpha/nu_0
dragBeta=Beta/nu_0


# ----- BOUNDARY CONDITIONS ----- #

for bc in tank.BC_list:
    bc.setFreeSlip()

tank.BC.top.setOpenAir()

tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)

tank.BC.bottom.setFreeSlip()

tank.BC.right.setFreeSlip()

tank.BC.sponge.setNonMaterial()



# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #


tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setPorousZones(flags=3, epsFact_solid=float((x5-x3)/2.),
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(x7-L_rightSpo/2.), 0., 0.),
                        )

# ----- Output Gauges ----- #
T = 40*period
gauge_dx=0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)

LG=[]
#PG=[]
PG0=[]

for i in probes:
    LG.append(((i,0.,0.), (i,tank_dim[1],0.)),)
    #PG.append((i,0.5,0.),)
    PG0.append((i, 0., 0.),)

point_output=ga.PointGauges(gauges=((('p'),PG0),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='point_gauges.csv')


#line_output=ga.LineGauges(gauges=((('u', 'v'), LG),
 #                                 (('p'), LG),
  #                               ),
   #                       activeTime = (0., T),
    #                      sampleRate=0.,
     #                     fileName='line_gauges.csv')

columnLines=tuple(map(tuple,LG))
fields=(('vof',))
integral_output=ga.LineIntegralGauges(gauges=((fields, columnLines),),
                                      activeTime = (0., T),
                                      sampleRate=0.,
                                      fileName='line_integral_gauges.csv')

domain.auxiliaryVariables += [point_output, integral_output
                             ]

##########################################
# Numerical Options and other parameters #
##########################################

he = waveinput.wavelength/opts.refinement_level
domain.MeshOptions.he = he 


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

st.assembleDomain(domain)

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = 0.1
dt_init = min(0.1*dt_fixed,0.001)
T = T
nDTout= int(round(T/dt_fixed))
runCFL = opts.cfl

#----------------------------------------------------
#  Discretization -- input options
#----------------------------------------------------

checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=opts.freezeLevelSet
useOnlyVF = False # if TRUE  proteus uses only these modules --> twp_navier_stokes_p + twp_navier_stokes_n
                  #                                              vof_p + vof_n
movingDomain=False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988

genMesh=True

# By DEFAULT on the other files.py -->  fullNewtonFlag = True 
#                                       multilevelNonlinearSolver & levelNonlinearSolver == NonlinearSolvers.Newton

useOldPETSc=False # if TRUE  --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.PETSc
                  # if FALSE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.KSP_petsc4py

useSuperlu = False #if TRUE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.LU

spaceOrder = 1
useHex     = False # used for discretization, if 1.0 --> CubeGaussQuadrature
                   #                          ELSE   --> SimplexGaussQuadrature

useRBLES   = 0.0 # multiplied with subGridError
useMetrics = 1.0 # if 1.0 --> use of user's parameters as (ns_shockCapturingFactor, ns_lag_shockCapturing, ecc ...)
useVF = opts.useVF # used in the smoothing functions as (1.0-useVF)*smoothedHeaviside(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf))


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
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.25 # magnifies numerical viscosity in NS (smoothening velocity fields)
    ns_lag_shockCapturing = True # lagging numerical viscosity speedsup Newton but destabilzes the solution
    ns_lag_subgridError = True # less nonlinear but less stable
    ls_shockCapturingFactor  = 0.25 # numerical diffusion of level set (smoothening phi)
    ls_lag_shockCapturing = True # less nonlinear but less stable
    ls_sc_uref  = 1.0 # reference gradient in numerical solution (higher=more diffusion)
    ls_sc_beta  = 1.0 # 1 is fully nonlinear, 2 is linear
    vof_shockCapturingFactor = 0.25 # numerical diffusion of level set (smoothening volume of fraction)
    vof_lag_shockCapturing = True # less nonlinear but less stable
    vof_sc_uref = 1.0 
    vof_sc_beta = 1.0 
    rd_shockCapturingFactor  = 0.25
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0 # control width of water/air transition zone
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0 # affects smoothing diffusion in mass conservation
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True # False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
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
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = max(1.0e-12,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)

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
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))

tank.BC.top.p_dirichlet = twpflowPressure_init
