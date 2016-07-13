from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.315, "Height of free surface above bottom"), # Choose waterLevels=[0.315]
    # Geometry of the tank - left lower boundary at (0.,0.,0.)
    ("Ls", 2., "Distance of the front toe of the structure end from generation zone in wavelengths"),
    ("Lend", 3.5, "Distance of the back toe of the structure end from absorption zone in meters"),
    ("Lgen", 1., "Length of generation zone in wavelegths"),
    ("Labs", 2., "Length of absorption zone in wavelegths"),
    ("h", 0.75, "Height of domain in meters"),
    # Geometry of a trapezoidal breakwater
    ("b", 0.25, "Width of breakwater at the top"), # Choose b=[0.25]
    ("hs", 0.25, "Height of the breakwater"),
    ("slope1", 1./2., "Slope1 of the breakwater"),
    ("slope2", 1./1.5, "Slope2 of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', 0.058, "Mean diameter of the medium"),
    # waves
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 1., "Period of the waves"), # Choose periods=[1., 1.5, 2.]
    ("wave_height", 0.03, "Height of the waves"), # Choose waveHeights=[0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12]
    ('wavelength', 1.393, 'Wavelength only if Fenton is activated'),
    ('Ycoeff', [0.06742994, 0.00358822, 0.00023583, 0.00001839, 0.00000158, 0.00000014, 0.00000001, 0.00000000], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [0.07136193, 0.00097421,-0.00000844,-0.00000031, 0.00000001, 0.00000000, 0.00000000, 0.00000000], 'Bcoeff only if Fenton is activated'),
    # numerical options
    ("he", .01 ,"he mesh size"),
    ("cfl", 0.9 ,"Target cfl"),
    ("freezeLevelSet", False, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing - For porous media set to 0.0"),
    ("conservative_Flux", False, 'For porous interface problem should be set eqaul to False'),
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
if opts.waveType=='Linear':
    waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=None, # if wave is linear I can use None
                                  waveType=opts.waveType)

if opts.waveType=='Fenton':
    waveinput = wt.MonochromaticWaves(period=period,
                                      waveHeight=waveHeight,
                                      mwl=mwl,
                                      depth=waterLevel,
                                      g=g,
                                      waveDir=waveDir,
                                      wavelength=opts.wavelength, # if wave is linear I can use None
                                      waveType=opts.waveType,
                                      Ycoeff=opts.Ycoeff,
                                      Bcoeff=opts.Bcoeff,
                                      )


#---------Domain Dimension
nd = 2
wl = waveinput.wavelength

# ----- SHAPES ----- #

L_leftSpo  = opts.Lgen*wl
L_rightSpo = opts.Labs*wl
hs=opts.hs
b=opts.b
slope1=opts.slope1
slope2=opts.slope2

x1=L_leftSpo
x2=x1+opts.Ls*wl
x3=x2+(hs/slope1)
x4=x3+b
x5=x4+(hs/slope2)
x6=x5+opts.Lend*wl
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


##############################################################################################################################
# POROUS MEDIA
##############################################################################################################################

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
dragAlpha=(porosity**2)*Alpha/nu_0
dragBeta=0.0 #(porosity**3)*Beta/nu_0


# ----- BOUNDARY CONDITIONS ----- #

for bc in tank.BC_list:
    bc.setFreeSlip()

tank.BC.top.setOpenAir()

tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)

tank.BC.bottom.setFreeSlip()

tank.BC.right.setNoSlip() #FreeSlip()

tank.BC.sponge.setNonMaterial()



# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #


tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setPorousZones(flags=3, dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(x7-L_rightSpo/2.), 0., 0.),
                        )

##############################################################################################################################
# ----- Output Gauges ----- #
##############################################################################################################################

T = 50*period

# ----- Layout 1
probes1=[]
probes1.append(((vertices[2][0]-0.47, 0.0, 0.0), (vertices[2][0]-0.47, tank_dim[1], 0.0)),)
probes1.append(((vertices[2][0]-0.35, 0.0, 0.0), (vertices[2][0]-0.35, tank_dim[1], 0.0)),)
for i in range(20):
    probes1.append(((vertices[2][0]+1.125+0.85+i*0.1, 0.0, 0.0), (vertices[2][0]+1.125+0.85+i*0.1, tank_dim[1], 0.0)),)


# ----- Layout 2
probes2=[]
probes2.append(((vertices[2][0]+1.125+1.50, 0.0, 0.0),                (vertices[2][0]+1.125+1.50, tank_dim[1], 0.0)),)   
probes2.append(((vertices[2][0]+1.125+1.50+0.56, 0.0, 0.0),           (vertices[2][0]+1.125+1.50+0.56, tank_dim[1], 0.0)),)  
probes2.append(((vertices[2][0]+1.125+1.50+0.56+0.30, 0.0, 0.0),      (vertices[2][0]+1.125+1.50+0.56+0.30, tank_dim[1], 0.0)),)  
probes2.append(((vertices[2][0]+1.125+1.50+0.56+0.30+0.16, 0.0, 0.0), (vertices[2][0]+1.125+1.50+0.56+0.30+0.16, tank_dim[1], 0.0)),)  


# ----- GenZone
genProbes=[]
dx=L_leftSpo/10.
for j in range(10+1):
    genProbes.append(((j*dx, 0.0, 0.0), (j*dx, tank_dim[1], 0.0)),)

columnLines1=tuple(map(tuple,probes1))
columnLines2=tuple(map(tuple,probes2))
columnLinesG=tuple(map(tuple,genProbes))
fields=(('vof',))
integral_output1=ga.LineIntegralGauges(gauges=((fields, columnLines1),),
                                      activeTime = (0., T),
                                      sampleRate=0.,
                                      fileName='ProbesConfiguration1.csv')


integral_output2=ga.LineIntegralGauges(gauges=((fields, columnLines2),),
                                      activeTime = (0., T),
                                      sampleRate=0.,
                                      fileName='ProbesConfiguration2.csv')


integral_output3=ga.LineIntegralGauges(gauges=((fields, columnLinesG),),
                                      activeTime = (0., T),
                                      sampleRate=0.,
                                      fileName='ProbesGeneration.csv')



domain.auxiliaryVariables += [integral_output1,integral_output2,integral_output3]

##########################################
# Numerical Options and other parameters #
##########################################

he = opts.he
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
    ns_shockCapturingFactor  = 0.5 # magnifies numerical viscosity in NS (smoothening velocity fields)
    ns_lag_shockCapturing = True # lagging numerical viscosity speedsup Newton but destabilzes the solution
    ns_lag_subgridError = True # less nonlinear but less stable
    ls_shockCapturingFactor  = 0.5 # numerical diffusion of level set (smoothening phi)
    ls_lag_shockCapturing = True # less nonlinear but less stable
    ls_sc_uref  = 1.0 # reference gradient in numerical solution (higher=more diffusion)
    ls_sc_beta  = 1.5 # 1 is fully nonlinear, 2 is linear
    vof_shockCapturingFactor = 0.5 # numerical diffusion of level set (smoothening volume of fraction)
    vof_lag_shockCapturing = True # less nonlinear but less stable
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0 # control width of water/air transition zone
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0 # affects smoothing diffusion in mass conservation
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True # False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
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

# Initial condition
waterLine_x = 2*tank_dim[0]
waterLine_z = waterLevel


def waveHeight(x,t):
    waterDepth = waveinput.eta(x, t) + waveinput.mwl
    return waterDepth


def wavePhi(x,t):
    [nd-1]- waveHeight(x,t)
    

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[nd-1]-waterLine_z

    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)

