from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 1., "Height of free surface above bottom"), 
    # waves
    ('waveType', 'Linear', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 1.94, "Period of the waves"), 
    ("wave_height", 0.025, "Height of the waves"),
    # porous medium 
    ('porosity', 0.4, "Porosity of the medium"),
    ('d_grain', 0.0102, "grain size"),
    ('aDarcy', 150, "Ergun Recommended values"),
    ('bForc', 0.0, "Forchheimer term is 0 for this case"),
    # numerical options
    ("refinement_level", 100. ,"he=walength/refinement_level"),
    ("cfl", 0.5 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
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

if opts.waveType=='Linear':
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
he = waveinput.wavelength/opts.refinement_level # MESH SIZE	


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = waveinput.wavelength
L_rightSpo = waveinput.wavelength


#-Tank
x1=L_leftSpo
x2=x1+1.*waveinput.wavelength
x3=x2+3.*waveinput.wavelength
tank_dim = [x3, 2.0]

boundaryOrientations = {'bottom': [0., -1.,0.],
                            'right': [1., 0.,0.],
                            'top': [0., 1.,0.],
                            'left': [-1., 0.,0.],
                            'sponge': None,
                            'porousLayer': None,
                            'moving_porousLayer': None,
                           }
boundaryTags = {'bottom': 1,
                    'right': 2,
                    'top': 3,
                    'left': 4,
                    'sponge': 5,
                    'porousLayer': 6,
                    'moving_porousLayer': 7,
                       }

##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2,  0.0],#2
              [x3,  0.0],#3
              [x3,  tank_dim[1] ],#4
              [x2,  tank_dim[1]],#5
              [x1,  tank_dim[1]],#6
              [0.0,  tank_dim[1]],#7
              ]

vertexFlags=np.array([1, 1, 1, 1,
                          3, 3, 3, 3,
                         ])

segments=[[0,1],
              [1,2],
              [2,3],
              [3,4],
              [4,5],
              [5,6],
              [6,7],
              [7,0],

              [1,6],
              [2,5],
             ]

segmentFlags=np.array([1, 1, 1,
                       2, 3, 3, 3, 4,
                       5, 6,
                          ])

regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.10*tank_dim[1] ],
            [ 0.90*x3 , 0.90*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


##################################################################################################################################################################################################################
# POROUS MEDIA
##################################################################################################################################################################################################################

"""
porosity=opts.porosity
voidFrac=1.0-porosity

d50=opts.d50
if d50==None:
    d15=opts.d15
else:
    d15=d50/1.2
aDarcy = 150.
bFoc = 0.

#----- SHIH
if opts.Resistance=='Shih':
    term1=3.12*(10**-3.)
    term2=(gAbs/(nu_0**2.))**(2./3.)
    term3=(d15**2.)
    Alpha1=1684+term1*term2*term3 #Shih
    Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2)) 

    term1=-5.10*(10**-3.)
    term2=(gAbs/(nu_0**2.))**(1./3.)
    term3=(d15)
    Beta1=1.72+1.57*exp(term1*term2*term3) #Shih
    Beta=Beta1*voidFrac/((porosity**3)*d15)

#----- ERGUN
if opts.Resistance=='Ergun':
    Alpha1=150 #Ergun
    Beta1=1.75 #Ergun
    Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2))
    Beta=Beta1*voidFrac/((porosity**3)*d15)

#----- ENGELUND
if opts.Resistance=='Engelund':
    Alpha1=360 #Ergun
    Beta1=3.6 #Ergun
    Alpha=Alpha1*nu_0*(voidFrac**3)/((porosity**2)*(d15**2))
    Beta=Beta1*voidFrac/((porosity**3)*d15)

#Proteus scale in viscosity, so i need to divide alpha and beta by nu_0
dragAlpha=(porosity**2)*1.96*(10**6) #(porosity**2)*Alpha/nu_0
dragBeta=0.0#(porosity**3)*Beta/nu_0
"""

#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

tank.BC.top.setOpenAir()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)
tank.BC.bottom.setFreeSlip()
tank.BC.right.setFreeSlip()
tank.BC.sponge.setNonMaterial()

tank.BC.porousLayer.setNonMaterial()

########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################


tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
"""                        )
tank.setPorousZones(flags=3, dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
"""
############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = 15.*period

gauge_dx=0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)
PG=[]
zProbes=waterLevel*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

point_output=ga.PointGauges(gauges=((('p'),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='point_gauges.csv')

domain.auxiliaryVariables += [ point_output,
                             ]

######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

he = he
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
movingDomain= not True
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

