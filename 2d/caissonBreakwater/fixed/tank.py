from __future__ import print_function
from __future__ import division
from builtins import map
from past.utils import old_div
from proteus import Domain, Context
#from proteus.mprans
import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.425, "Height of free surface above bottom"), # Choose waterLevels=[0.425, 0.463]
    # waves
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 0.5, "Period of the waves"), # Choose periods=[0.5, 0.8, 1.1, 1.5, 2.8, 3.9, 4.0]
    ("wave_height", 0.025, "Height of the waves"), # # Choose for d=0.425-->[0.025, 0.075, 0.125, 0.234]. Choose for d=0.463-->[0.025, 0.075, 0.125, 0.254].
    ('wavelength', 0.4, 'Wavelength only if Fenton is activated'), # Choose for d=0.425-->[0.4, 1.0, 1.8, 2.9]. Choose for d=0.463-->[0.4, 1.0, 1.8, 2.9, 3.0, 5.7, 5.9, 8.8, 9.4].
    ('Ycoeff', [0.19167938     ,  0.01943414     ,  0.00299676     ,  0.00055096     ,  0.00011165  ,  0.00002413     ,  0.00000571     ,  0.00000251], 'Ycoeff only if Fenton is activated'), 
    ('Bcoeff', [0.19063009     ,  0.00072851     ,  0.00002905     ,  0.00000131     ,  0.00000006  ,  0.00000000     ,  0.00000000     ,  0.00000000], 'Bcoeff only if Fenton is activated'), 
    # Geometry of the tank - left lower boundary at (0.,0.,0.)
    ("Ls",   2.0, "Distance of the front toe of the structure end from generation zone in wavelengths"),
    ("Lend", 2.0, "Distance of the back toe of the structure end from absorption zone in wavelengths"),
    ("Lgen", 1., "Length of generation zone in wavelegths"),
    ("Labs", 2., "Length of absorption zone in wavelegths"),
    ("h", 1.0, "Height of domain in meters"),
    # breakwater
    ("hs", 0.075, "Height of the breakwater"),
    ("slope", old_div(1.,2.), "Slope of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', None, "Mean diameter of the medium"),
    ('d15', 0.038, "15% grading curve diameter of the medium"),
    # caisson
    ("caisson", True, "Switch on/off caisson"),
    ('dimx', 0.44, 'X-dimension of the caisson'), # Choose dimx=[0.44]
    ('dimy', 0.4, 'Y-dimension of the caisson'), 
    ("rotation", not True, "Initial position for free oscillation"),
    ("friction", not True, "Switch on/off friction module"),
    ("m_static", 0.0, "Static friction factor between caisson and rubble mound"),
    ("m_dynamic", 0.0, "Dynamic friction factor between caisson and rubble mound"),
    # numerical options
    ("refinement_level", 200. ,"he=walength/refinement_level"),
    ("cfl", 0.9 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing"),
    ('movingDomain', not True, "Moving domain and mesh option"),
    ('conservativeFlux', not True,'Fix post-processing velocity bug for porous interface'),
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

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

L_leftSpo  = opts.Lgen*wl
L_rightSpo = opts.Labs*wl

hs=opts.hs
slope=opts.slope

#-Caisson
dimx=opts.dimx
dimy=opts.dimy
b=dimx

#-Tank
x1=L_leftSpo
x2=x1+opts.Ls*wl
x3=x2+(old_div(hs,slope))

xc1=x3+0.13
xc2=xc1+b
yc1=yc2=hs

x4=xc2+0.13
x5=x4+(old_div(hs,slope))
x6=x5+opts.Lend*wl
x7=x6+L_rightSpo
tank_dim = [x7, opts.h]

boundaryOrientations = {'bottom': [0., -1.,0.],
                        'right': [1., 0.,0.],
                        'top': [0., 1.,0.],
                        'left': [-1., 0.,0.],
                        'sponge': None,
                        'moving_sponge': None,
                       }
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,
                'sponge': 5,
                'moving_sponge': 6,
               }

##############################################################################################################################################################################################################
# Caisson 
############################################################################################################################################################################################################

if opts.caisson:
    dimx=dimx
    dimy=dimy
    dim=(dimx,dimy)
    coords=[xc1+old_div(b,2.), hs+old_div(dimy,2.)] # For bodyDimensions and barycenter
    VCG=old_div(dim[1],2.)                 # For barycenter
    width=1.0                      # The 3rd dimension
    density=100000 #kg/m3
    volume=dimx*dimy*width
    mass=density*volume
    It=old_div((dimx**2.+dimy**2.),12.)

    caisson2D = st.Rectangle(domain, dim=dim, coords=coords)
    caisson2D.vertices[0][0]=xc1
    caisson2D.vertices[0][1]=yc1
    caisson2D.vertices[1][0]=xc2
    caisson2D.vertices[1][1]=yc2   

    free_x=(1.0, 1.0, 1.0) # Translational DOFs
    free_r=(0.0, 0.0, 1.0) # Rotational DOFs
    m_static=opts.m_static # Static friction
    m_dynamic=opts.m_dynamic # Dynamic friction

    if opts.movingDomain==True:
        free_x=(1.0, 1.0, 0.0) # Translational DOFs
        free_r=(0.0, 0.0, 1.0) # Rotational DOFs
   
    caisson2D.setMass(mass)
    caisson2D.setConstraints(free_x=free_x, free_r=free_r)
    caisson2D.setFriction(friction=opts.friction, m_static=m_static, m_dynamic=m_dynamic)

    if opts.rotation==True: # Initial position for free oscillation
        caisson2D.rotate(rotation)

    caisson2D.It= It/caisson2D.mass/width
    caisson2D.setRecordValues(all_values=True)
    caisson2D.setRigidBody()

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

    vertexFlags=np.array([1, 1, 1, 
                          5, 5, 
                          1, 1, 1,
                          3, 3, 3, 3,
                          6, 6,
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
                           5, 5,
                           1,
                           6, 6,
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
    d15=old_div(d50,1.2)

term1=3.12*(10**-3.)
term2=(old_div(gAbs,(nu_0**2.)))**(old_div(2.,3.))
term3=(d15**2.)
Alpha1=1684+term1*term2*term3 #Shih
#Alpha1=150 #Ergun
#Alpha1=360 #Engelund

term1=-5.10*(10**-3.)
term2=(old_div(gAbs,(nu_0**2.)))**(old_div(1.,3.))
term3=(d15)
Beta1=1.72+1.57*exp(term1*term2*term3) #Shih
#Beta1=1.75 #Ergun
#Beta1=3.6 #Engelund

#Alpha=Alpha1*nu_0*(voidFrac**3)/((porosity**2)*(d15**2))  #Engelund
Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2))  #Ergun
Beta=Beta1*voidFrac/((porosity**3)*d15)

#Proteus scale in viscosity, so i need to divide alpha and beta by nu_0
dragAlpha=(porosity**2)*Alpha/nu_0
dragBeta=0.0#(porosity**3)*Beta/nu_0


#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

if opts.caisson:
    for bc in caisson2D.BC_list:
        bc.setFreeSlip()

tank.BC.top.setOpenAir()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)
#tank.BC.bottom.setFreeSlip()
tank.BC.bottom.setNoSlip()
#tank.BC.right.setFreeSlip()
tank.BC.right.setNoSlip()
tank.BC.sponge.setNonMaterial()
tank.BC.moving_sponge.setNonMaterial()

if opts.movingDomain==True:
    for tb in [tank.BC.right, tank.BC.left, tank.BC.top, tank.BC.bottom, tank.BC.sponge]:
        tb.hx_dirichlet= lambda x, t: 0.0
        tb.hy_dirichlet= lambda x, t: 0.0
        tb.hz_dirichlet= lambda x, t: 0.0
        tb.u_stress=None
        tb.v_stress=None
        tb.w_stress=None
    ms=tank.BC.moving_sponge
    ms.hx_dirichlet= None
    ms.hy_dirichlet= None
    ms.hz_dirichlet= lambda x, t: 0.0
    ms.u_stress=None
    ms.v_stress=None
    ms.w_stress=None



########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################



tank.setGenerationZones(flags=1, epsFact_solid=float(old_div(L_leftSpo,2.)),
                        orientation=[1., 0.], center=(float(old_div(L_leftSpo,2.)), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setPorousZones(flags=3, epsFact_solid=float(old_div((x5-x2),2.)),
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=float(old_div(L_rightSpo,2.)),
                        orientation=[-1., 0.], center=(float(x7-old_div(L_rightSpo,2.)), 0., 0.),
                        )


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = 20.*period

PG=[]
LG=[]
LG1=[]
LG2=[]
LG3=[]
LG4=[]

#-------------Pressure and vof gauges for WaveHeight--------------------------------#


z_probes = waterLevel*0.5
PG=[(0.0,z_probes,0.), (xc1-0.73,z_probes,0.), (xc1+0.53,z_probes,0.), (xc1+1.03,z_probes,0.), (xc1+1.23,z_probes,0.)]
pressureGauges=ga.PointGauges(gauges=((('p',),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='pressureProbes.csv')


VG=[((0.0,0.,0.),(0.0,tank_dim[1],0.)),
    ((xc1-0.73,0.,0.),(xc1-0.73,tank_dim[1],0.)),
    ((xc1+0.53,0.,0.),(xc1+0.53,tank_dim[1],0.)),
    ((xc1+1.03,0.,0.),(xc1+1.03,tank_dim[1],0.)), 
    ((xc1+1.23,0.,0.),(xc1+1.23,tank_dim[1],0.))]

VG=tuple(map(tuple,VG))
fields=(('vof',))
vof_probes=ga.LineIntegralGauges(gauges=((fields, VG),
                                               ),
                                      activeTime = (0., T),
                                      sampleRate=0.,
                                      fileName='waveProbes.csv')


if opts.caisson:
    xc1=caisson2D.vertices[0][0]-1*(10**-5) #to avoid floating point error
    yc1=caisson2D.vertices[0][1]-1*(10**-5) #to avoid floating point error
    xc2=caisson2D.vertices[1][0]+1*(10**-5) #to avoid floating point error
    yc2=caisson2D.vertices[1][1]-1*(10**-5) #to avoid floating point error
    xc3=caisson2D.vertices[2][0]+1*(10**-5) #to avoid floating point error
    yc3=caisson2D.vertices[2][1]+1*(10**-5) #to avoid floating point error
    xc4=caisson2D.vertices[3][0]-1*(10**-5) #to avoid floating point error
    yc4=caisson2D.vertices[3][1]+1*(10**-5) #to avoid floating point error
    
#-------------Wave overtopping-----------------------------------------------#
    
    dx=0.01
    probes=np.linspace(yc1, tank_dim[1], old_div((tank_dim[1]-yc1),dx)+1)    
    for i in probes:
        LG.append((xc1,i,0.),)
    overtoppingGauges=ga.PointGauges(gauges=((('u',), LG),),
                              activeTime = (0., T),
                              sampleRate=0.,
                              fileName='overtoppingVelGauges.csv')
    
    vofGauges=ga.PointGauges(gauges=(((('vof'),), LG),),
                              activeTime = (0., T),
                              sampleRate=0.,
                              fileName='overtoppingVofGauges.csv')


#-------------Wave loading---------------------------------------------------#

    probes1=np.linspace(xc1, xc2, old_div((xc2-xc1),dx)+1) 
    for i in probes1:
        LG1.append((i,yc1,0.),)    

    probes2=np.linspace(yc2, yc3, old_div((yc3-yc2),dx)+1) 
    for i in probes2:
        LG2.append((xc2,i,0.),)   

    probes3=np.linspace(xc4, xc3, old_div((xc3-xc4),dx)+1) 
    for i in probes3:
        LG3.append((i,yc3,0.),)    

    probes4=np.linspace(yc1, yc4, old_div((yc4-yc1),dx)+1) 
    for i in probes4:
        LG4.append((xc4,i,0.),)  

    loadingsGauges=ga.PointGauges(gauges=((('p',), LG1),
                                         (('p',), LG2),
                                         (('p',), LG3),
                                         (('p',), LG4),),
                              activeTime = (0., T),
                              sampleRate=0.,
                              fileName='loadingGauges.csv')


    domain.auxiliaryVariables += [pressureGauges,vof_probes,overtoppingGauges,vofGauges,
                                  loadingsGauges,
                                 ]
else:
    domain.auxiliaryVariables += [pressureGauges, vof_probes,
                                 ]



######################################################################################################################################################################################################################
# Numerical Options and other parameters #
######################################################################################################################################################################################################################

he = old_div(waveinput.wavelength,opts.refinement_level)
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
weak_bc_penalty_constant = old_div(10.0,nu_0) #100
dt_fixed = 0.1
dt_init = min(0.1*dt_fixed,0.001)
T = T
nDTout= int(round(old_div(T,dt_fixed)))
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
movingDomain=opts.movingDomain
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
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES)
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
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
