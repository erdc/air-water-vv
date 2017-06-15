from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.4, "Height of free surface above bottom"),
    # Geometry
    ('Lgen', 0.5, 'Genaration zone in terms of wave lengths'),
    ('Labs', 0.5, 'Absorption zone in terms of wave lengths'),
    ('Ls', 0.5, 'Length of domain from genZone to the front toe of rubble mound in terms of wave lengths'),
    ('Lend', 0.5, 'Length of domain from absZone to the back toe of rubble mound in terms of wave lengths'),
    ('th', 0.6, 'Total height of numerical tank'),
    ('wl', 0.75, 'Used instead of wavelength from wavetools. Spatial parameter to build-up the doamin only if wave is False'),    
    # waves
    ('wave', not True, 'Turn on wave generation'),
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 0., "Period of the waves"),
    ("wave_height", 0, "Height of the waves"), 
    ('wavelength', 0., 'Wavelength only if Fenton is activated'), 
    ('Ycoeff', [ 0. ], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [ 0. ], 'Bcoeff only if Fenton is activated'),
    ('Nf', 0 ,'Number of frequency components for fenton waves'),
    ('meanVelocity', [ 0.1, 0., 0.],'Velocity used for currents'),
    ('phi0', 0.0 ,'Initial phase for waves'),
    ('Uwind', [0.0, 0.0, 0.0], 'Set air velocity'),
    # bridge2D
    ("bridge2D", True, "Switch on/off bridge2D"),
    ('maxLength', 0.3381502, 'Maximum length of the bridge'),
    ('maxHeight', 0.05447284, 'Maximum height of the bridge'),  
    ('modelScale', 1.0, 'Use it for scaling the model'),
    ('draught', 0.04, 'Draught of the bridge (enter a positive number)'),    
    ('width', 1.0, 'Z-dimension of the bridge2D'),
    ('mass', 5.0, 'Mass of the caisson2D [kg]'),    # density of polypropylene is 946 kg/m3
    ("rotation", not True, "Initial position for free oscillation"),
    ('bridgeBC', 'WallFunction', 'bridge2D boundaries: NoSlip or FreeSlip or WallFunction'),
    ('InputMotion', not True, 'If True, set a motion as input rather than calculate it'),
    ('scheme', 'Forward_Euler', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Forward_Euler)'),
    ("springs", not True, "Switch on/off soil module"),
    ("Kx", 0.0, "Horizontal stiffness in N/m"),
    ("Ky", 0.0, "Vertical stiffness in N/m"),
    ("Krot", 0.0, "Rotational stiffness in N*m"),
    ("Cx", 0.0, "Damping factor in N/m*s "),
    ("Cy", 0.0, "Damping factor in N/m*s "),
    ("Crot", 0.0, "Rotational damping factor in N*m*s "),
    # Turbulence
    ("useRANS", 3, "Switch ON turbulence models: 0-None, 1-K-Epsilon, 2-K-Omega1998, 3-K-Omega1988"),
    ("ns_closure", 4, "1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega"),
    ("c_mu", 0.09, "mu coefficient for the turbulence model"), 
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("sigma_e", 1.0, "sigma_e coefficient for the turbulence model"),
    ("Y", 0.01, "Y used for y+ calculation"),
    # numerical options
    ("GenZone", True, 'Turn on generation zone at left side'),
    ("AbsZone", True, 'Turn on absorption zone at right side'),
    ('duration', 5., 'Simulation duration'),
    ("refinement_level", 0.0,"he=walength/refinement_level"),
    ("he", 0.01,"Mesh size"),
    ("cfl", 0.90 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 0.0, "For density and viscosity smoothing"), # used in the smoothing functions as (1.0-useVF)*smoothedHeaviside(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf))
    ('movingDomain', not True, "Moving domain and mesh option"),
    ('conservativeFlux', True,'Fix post-processing velocity bug for porous interface'),
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
nu_0=1.004e-6
rho_1=1.205
nu_1=1.500e-5

sigma_01=0.0
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))


# ----- WAVE input ----- #
if opts.wave == True:
    waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=opts.wavelength,       # used by fenton waves
                                  waveType=opts.waveType, 
                                  Ycoeff=np.array(opts.Ycoeff),     # used by fenton waves
                                  Bcoeff=np.array(opts.Bcoeff),     # used by fenton waves
                                  Nf=opts.Nf,                       # used by fenton waves
                                  meanVelocity = np.array(opts.meanVelocity),
                                  phi0 = opts.phi0,
                                      )

#---------Domain Dimension

nd = 2.
wl=opts.wl
th = opts.th

#---------MESH SIZE
if opts.he == 0.0:
    he = wl/opts.refinement_level
else:
    he = opts.he 

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

if opts.bridge2D:
    L_leftSpo  = opts.Lgen*wl
    L_rightSpo = opts.Labs*wl

#-bridge
    b=opts.maxLength*opts.modelScale
    h=opts.maxHeight*opts.modelScale
    
#-Tank
    x1=L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+b
    x4=x3+opts.Lend*wl
    x5=x4+L_rightSpo

    xc1 = x2+(b*0.5)
    yc1 = opts.water_level

    tank_dim = [x5, opts.th]

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1, 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'bridge': None,
                           }
boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                    'bridge':6,
                       }

##############################################################################################################################################################################################################
# bridge2D
############################################################################################################################################################################################################
if opts.bridge2D:
    width=opts.width           # The 3rd dimension
    mass=opts.mass             # kg
    I = ((b**2)+(h**2))*mass/12.

    offsetX = x2
    offsetY = opts.water_level - opts.draught
    modelScale = opts.modelScale
    maxLength = opts.maxLength
    maxHeight = opts.maxHeight

# --- Shape properties setup
    xb1	    ,	yb1 	= 	(	0.3381502	 + offsetX)*modelScale	,	(	0.0544728	 + offsetY)*modelScale
    xb2	    ,	yb2	    = 	(	0.3381502	 + offsetX)*modelScale	,	(	0.0374244	 + offsetY)*modelScale
    xb3	    ,	yb3	    = 	(	0.3311144	 + offsetX)*modelScale	,	(	0.0230073	 + offsetY)*modelScale
    xb4	    ,	yb4	    = 	(	0.3201899	 + offsetX)*modelScale	,	(	0.0112598	 + offsetY)*modelScale
    xb5	    ,	yb5	    = 	(	0.3063240	 + offsetX)*modelScale	,	(	0.0031953	 + offsetY)*modelScale
    xb6	    ,	yb6	    = 	(	0.2947060	 + offsetX)*modelScale	,	(	0.0000000	 + offsetY)*modelScale
    xb7	    ,	yb7	    = 	(	0.1690751	 + offsetX)*modelScale	,	(	0.0000000	 + offsetY)*modelScale
    xb8	    ,	yb8	    = 	(	0.0434442	 + offsetX)*modelScale	,	(	0.0000000	 + offsetY)*modelScale
    xb9	    ,	yb9	    = 	(	0.0318262	 + offsetX)*modelScale	,	(	0.0031953	 + offsetY)*modelScale
    xb10	,	yb10	= 	(	0.0179603	 + offsetX)*modelScale	,	(	0.0112598	 + offsetY)*modelScale
    xb11	,	yb11	= 	(	0.0070358	 + offsetX)*modelScale	,	(	0.0230073	 + offsetY)*modelScale
    xb12	,	yb12	= 	(	0.0000000	 + offsetX)*modelScale	,	(	0.0374244	 + offsetY)*modelScale
    xb13	,	yb13	= 	(	0.0000000	 + offsetX)*modelScale	,	(	0.0544728	 + offsetY)*modelScale

    vertices	  =  	[	[	xb1	,	yb1	] ,	#	0
               		    	[	xb2	,	yb2	] ,	#	1
        	        		[	xb3	,	yb3	] ,	#	2
                       		[	xb4	,	yb4	] ,	#	3
                       	    [	xb5	,	yb5	] ,	#	4
             		    	[	xb6	,	yb6	] ,	#	5
                    		[	xb7	,	yb7	] ,	#	6
        	        		[	xb8	,	yb8	] ,	#	7
                    		[	xb9	,	yb9	] ,	#	8
                           	[	xb10	,	yb10	] ,	#	9
            		    	[	xb11	,	yb11	] ,	#	10
                		    [	xb12	,	yb12	] ,	#	11
                           	[	xb13	,	yb13	] ,	#	12
        	        	]							


    vertexFlags=np.array([4, 4, 6, 6, 6, 6,
                          1, 
                          2, 2, 6, 6, 6, 6, 
                          ])

    segments = [ 	[	0	,	1	] , 
	                [	1	,	2	] , 
                	[	2	,	3	] , 
                	[	3	,	4	] , 
                    [	4	,	5	] , 
                    [	5	,	6	] , 
                    [	6	,	7	] , 
                    [	7	,	8	] , 
                    [	8	,	9	] , 
                    [	9	,	10	] , 
                    [	10	,	11	] , 
                    [	11	,	12	] , 
                    [	12	,	0	] , 
                    ]				

    segmentFlags=np.array([4, 6, 6, 6, 6,
                           1, 1, 
                           6, 6, 6, 6, 2,
                           3,
                         ])
    holes = [ [ xc1 , yc1 ],
            ]

    regions = [ [ xc1 , yc1 ],
            ]

    regionFlags=np.array([1])

    bridge = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                          segments=segments, segmentFlags=segmentFlags,
                          holes=holes,
                          regions=regions, regionFlags=regionFlags,
                          boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)

# --- Body properties setup
    bridge2D = bd.RigidBody(shape=bridge)
    if opts.movingDomain==True:
        free_x=(0.0, 0.0, 0.0) # Translational DOFs
        free_r=(0.0, 0.0, 1.0) # Rotational DOFs
    bridge2D.setMass(mass)
    bridge2D.It= I/bridge2D.mass/width
    bridge2D.setConstraints(free_x=free_x, free_r=free_r)
    bridge2D.setNumericalScheme(scheme=opts.scheme)
    bridge2D.setRecordValues(filename='bridge2D', all_values=True)

#----- Spring setup
    if opts.springs:
        Kx = opts.Kx
        Ky = opts.Ky
        Krot = opts.Krot
        Cx = opts.Cx
        Cy = opts.Cy
        Crot = opts.Crot
        bridge2D.setSprings(opts.springs, Kx, Ky, Krot, Cx, Cy, Crot)


##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

if opts.bridge2D==True:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2,  0.0],#2
              [x3,  0.0],#3
              [x4,  0.0],#4
              [x5,  0.0],#5
              [x5,  tank_dim[1]],#6
              [x4,  tank_dim[1]],#7
              [x1,  tank_dim[1]],#8
              [0.0, tank_dim[1]],#9
              ]

    vertexFlags=np.array([1, 1, 1, 1, 1, 1,
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
              [9,0],

              [1,8],
              [4,7],
             ]

    segmentFlags=np.array([ 1, 1, 1, 1,   
                            1,                        
                            2, 3, 3, 3, 4,
                            5, 5,
                         ])


regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.10*tank_dim[1] ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)

# Hole
#if opts.bridge2D:
#    # let gmsh know that the caisson is IN the tank
#    # index says which region the child is going to be attached to.
#    tank.setChildShape(bridge, 0)


############################################################################################################################################################################
# ----- Turbulence ----- #
###########################################################################################################################################################################

d = maxLength #tank_dim[1] / 2.
L = 0. #tank_dim[0]
U0 = opts.meanVelocity[0] # freestream velocity
c_mu = opts.c_mu
Re0 = U0*d/nu_0
ReL = U0*L/nu_0

Uwind = opts.Uwind # opts.meanVelocity
dwind = th - opts.water_level #d
Re0wind = Uwind[0]*dwind/nu_1

# Schlichting
cf = 0.045*(Re0**(-1./4.))
cfwind = 0.
if Re0wind >0.: cfwind = 0.045*(Re0wind**(-1./4.))

Y_ = opts.Y

Ut = U0*sqrt(cf/2.)
kappaP = (Ut**2)/sqrt(c_mu)
dissipationP = (Ut**3)/(0.41*Y_)

# air
Utwind = Uwind[0]*sqrt(cfwind/2.)
kappaPwind = (Utwind**2)/sqrt(c_mu)
if kappaPwind == 0.: kappaPwind = 1e-10
dissipationPwind = (Utwind**3)/(0.41*Y_)
if dissipationPwind == 0.: dissipationPwind = 1e-10

# ke or kw
useRANS = opts.useRANS  # 0 -- None
                        # 1 -- K-Epsilon
                        # 2 -- K-Omega, 1998
                        # 3 -- K-Omega, 1988
if useRANS >= 2:
    # in kw model w = e/k
    dissipationP = dissipationP/kappaP
    dissipationP = dissipationPwind/kappaPwind

# inlet values 
kInflow = kappaP #* 0.0001 # None
dissipationInflow = dissipationP #* 0.0001 # None

kInflowAir = kappaPwind #* 0.0001 # None
dissipationInflowAir = dissipationPwind #* 0.0001 # None

#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

if opts.bridge2D:
    for bc in bridge.BC_list:
        if opts.bridgeBC == 'FreeSlip':
            bc.setFreeSlip()
        elif opts.bridgeBC == 'NoSlip':
            bc.setNoSlip()
        elif opts.bridgeBC == 'WallFunction':
            bc.setWallFunction_kw(Y=Y_, U0=U0, d=d, skinf='turbulent', B=5.2)          
            bc.u_diffusive.resetBC()
    # top bridge boundary
    if opts.bridgeBC == 'FreeSlip': bridge.BC['y+'].setFreeSlip()
    if opts.bridgeBC == 'NoSlip':   bridge.BC['y+'].setNoSlip()    
            
if opts.wave==True:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1)
    tank.BC['x+'].setFreeSlip()
else:
    tank.BC['x-'].setTwoPhaseVelocityInlet(U=opts.meanVelocity,  waterLevel=opts.water_level,
                                           smoothing=he*3.0, Uwind=Uwind,
                                           kInflow=kInflow, dissipationInflow=dissipationInflow,
                                           kInflowAir=kInflowAir, dissipationInflowAir=dissipationInflowAir)
    tank.BC['x+'].setHydrostaticPressureOutletWithDepth( seaLevel=opts.water_level, rhoUp=rho_1,
                                                         rhoDown=rho_0, g=g,
                                                         refLevel=opts.th, smoothing=he*3.0,
                                                         U=opts.meanVelocity, Uwind=Uwind )

tank.BC['y+'].setAtmosphere()
#tank.BC['y+'].setTurbulentZeroGradient()
tank.BC['y-'].setWallFunction_kw(Y=Y_, U0=U0, d=d, skinf='turbulent', B=5.2)
tank.BC['y-'].u_diffusive.resetBC()
tank.BC['sponge'].setNonMaterial()

if opts.movingDomain==True:
    for tb in tank.BC_list:
        tb.setFixedNodes()
    for bc in bridge.BC_list:
        bc.hx_dirichlet.uOfXT, bc.hy_dirichlet.uOfXT, bc.hz_dirichlet.uOfXT = None, None, None  
        bc.u_stress.uOfXT, bc.v_stress.uOfXT, bc.w_stress.uOfXT = None, None, None


########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################

if opts.GenZone:
    if opts.wave != True: waveinput = None
    tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput,
                        )

if opts.AbsZone:
    tank.setAbsorptionZones(flags=3, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(tank_dim[0]-L_rightSpo/2.), 0., 0.),
                        )

############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = opts.duration


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
dt_fixed = 0.2
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
movingDomain=opts.movingDomain

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
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = ecH = epsFact_density
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

