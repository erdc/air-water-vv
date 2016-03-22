from proteus import Domain, Context
#from proteus.mprans
import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.463, "Height of free surface above bottom"), # Choose waterLevels=[0.425, 0.463, 0.538]
    # waves
    ("wave_period", 0.5, "Period of the waves"), # Choose periods=[0.5, 0.8, 1.1, 1.5]
    ("wave_height", 0.025, "Height of the waves"), # Choose waveHeights=[0.025, 0.075, 0.125, 0.234]
    # breakwater
    ("hs", 0.075, "Height of the breakwater"),
    ("slope", 1./2., "Slope of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', None, "Mean diameter of the medium"),
    ('d15', 0.038, "15% grading curve diameter of the medium"),
    # caisson
    ("caisson", True, "Switch on/off caisson"),
    ('dimx', 0.44, 'X-dimension of the caisson'), # Choose dimx=[0.25,0.44,0.64]
    ('dimy', 0.4, 'Y-dimension of the caisson'), 
    ("rotation", False, "Initial position for free oscillation"),
    # numerical options
    ("refinement_level", 100. ,"he=walength/refinement_level"),
    ("cfl", 0.9 ,"Target cfl"),
    ("freezeLevelSet", False, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing"),
    ('movingDomain', False, "Moving domain and mesh option"),
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

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = waveinput.wavelength
L_rightSpo = waveinput.wavelength

hs=opts.hs
slope=opts.slope

#-Caisson
dimx=opts.dimx
dimy=opts.dimy
b=dimx

#-Tank
x1=L_leftSpo
x2=x1+2.*waveinput.wavelength
x3=x2+(hs/slope)

xc1=x3+0.13
xc2=xc1+b
yc1=yc2=hs

x4=xc2+0.13
x5=x4+(hs/slope)
x6=x5+2.*waveinput.wavelength
x7=x6+L_rightSpo
tank_dim = [x7, 1.0]

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

##############################################################################################################################################################################################################
# Caisson 
############################################################################################################################################################################################################

if opts.caisson:
    dimx=dimx
    dimy=dimy
    dim=(dimx,dimy)
    coords=[xc1+b/2., hs+dimy/2.] # For bodyDimensions and barycenter
    VCG=dim[1]/2.                 # For barycenter
    width=1.                      # The 3rd dimension
    mass=1.
    It=(dimx**2.+dimy**2.)/12.

    rotation=opts.rotation        # Initial position for free oscillation

    free_x=(0.0, 0.0, 0.0) # Translational DOFs
    free_r=(0.0, 0.0, 0.0) # Rotational DOFs
    if opts.movingDomain==True:
        free_x=(1.0, 1.0, 0.0) # Translational DOFs
        free_r=(0.0, 0.0, 1.0) # Rotational DOFs
    
    caisson2D = st.Rectangle(domain, dim=dim, coords=coords)

# Here I need to over-ride the points I define previously on the caisson's ones
# In order to avoid segmentation violation and the mixing of the regions!
# The reason is that st.Rectangle creates points with hs+0.0000000000000001!!! 
    caisson2D.vertices[0][0]=xc1
    caisson2D.vertices[0][1]=yc1
    caisson2D.vertices[1][0]=xc2
    caisson2D.vertices[1][1]=yc2

    caisson2D.setRigidBody()
    caisson2D.setMass(mass)
    caisson2D.setConstraints(free_x=free_x, free_r=free_r)
    if rotation:
        caisson2D.rotate(rotation)
    caisson2D.It= It
    caisson2D.setRecordValues(pos=True, rot=True, F=True, M=True)

##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

if opts.caisson==False:

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

    segmentFlags=np.array([1, 1,
                           5, 5, 5,
                           1, 1,
                           2, 3, 3, 3, 4,
                           5, 5, 1,
                          ])
else:
    
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
              [xc1, yc1],#12
              [xc2, yc2],#13
              ]

    vertexFlags=np.array([1, 1, 1, 1, 1, 1, 1, 1,
                          3, 3, 3, 3,
                          1, 1,
                         ])

    segments=[[0,1],
              [1,2],
              [2,3],

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
              [3,12],
              [13,4],              
             ]

    segmentFlags=np.array([1, 1,
                           5, 5,
                           1, 1,
                           2, 3, 3, 3, 4,
                           5, 5, 1,
                           5, 5,
                         ])

regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.10*tank_dim[1] ],
            [ xc1 , 0.50*hs ],
            [ 0.95*x7 , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3, 4])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


##################################################################################################################################################################################################################
# POROUS MEDIA
##################################################################################################################################################################################################################


porosity=opts.porosity
voidFrac=1.0-porosity

d50=opts.d50
if d50==None:
    d15=opts.d15
else:
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


#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################


if opts.caisson:
    for bc in caisson2D.BC_list:
        bc.setNoSlip()
    caisson2D.BC.bottom.setFreeSlip()


for bc in tank.BC_list:
    bc.setFreeSlip()
tank.BC.top.setOpenAir()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)
tank.BC.bottom.setFreeSlip()
tank.BC.right.setFreeSlip()
tank.BC.sponge.setNonMaterial()


########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################



tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setPorousZones(flags=3, epsFact_solid=float((x5-x2)/2.),
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(x7-L_rightSpo/2.), 0., 0.),
                        )


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################

T = 40*period

PG=LG=LG1=LG2=LG3=LG4=[]


#-------------Pressure gauges for WaveHeight--------------------------------#

PG=[(x3-0.60,0.5,0.), (x5,0.5,0.), (x5+0.50,0.5,0.), (x5+0.70,0.5,0.)]
PG=tuple(map(tuple,PG))
pressureGauges=ga.PointGauges(gauges=((('p'),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='pressureGauges.csv')

if opts.caisson:
    xc1=caisson2D.vertices[0][0]
    yc1=caisson2D.vertices[0][1]-1*(10**-10) #to avoid floating point error
    xc2=caisson2D.vertices[1][0]
    yc2=caisson2D.vertices[1][1]-1*(10**-10) #to avoid floating point error
    xc3=caisson2D.vertices[2][0]
    yc3=caisson2D.vertices[2][1]+1*(10**-10) #to avoid floating point error
    xc4=caisson2D.vertices[3][0]
    yc4=caisson2D.vertices[3][1]+1*(10**-10) #to avoid floating point error

#-------------Wave overtopping-----------------------------------------------#

    LG=(((xc4,yc4,0.), (xc4,tank_dim[1],0.)),)

    overtoppingGauges=ga.LineGauges(gauges=((('u','v'), LG),
                                           ),
                                    activeTime = (0., T),
                                    sampleRate=0.,
                                    fileName='overtoppingGauges.csv')

    vofGauges=ga.LineGauges(gauges=((('vof'), LG),
                                   ),
                            activeTime = (0., T),
                            sampleRate=0.,
                            fileName='vofGauges.csv')


#-------------Wave loading---------------------------------------------------#


    LG1=(((xc1,yc1,0.), (xc2,yc2,0.)),)
    LG2=(((xc2,yc2,0.), (xc3,yc3,0.)),)
    LG3=(((xc3,yc3,0.), (xc4,yc4,0.)),)
    LG4=(((xc4,yc4,0.), (xc1,yc1,0.)),)

    loadingsGauges=ga.LineGauges(gauges=((('p'), LG1),
                                         (('p'), LG2),
                                         (('p'), LG3),
                                         (('p'), LG4),
                                     ),
                              activeTime = (0., T),
                              sampleRate=0.,
                              fileName='loadingsGauges.csv')


    domain.auxiliaryVariables += [pressureGauges,overtoppingGauges,loadingsGauges,vofGauges,
                                 ]
else:
    domain.auxiliaryVariables += [pressureGauges,
                                 ]


######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

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
