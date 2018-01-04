from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.325, "Height of free surface above bottom"),
    # Geometry
    ('Lgen', 1.0, 'Genaration zone in terms of wave lengths'),
    ('Labs', 1.0, 'Absorption zone in terms of wave lengths'),
    ('Ls', 1.0, 'Length of domain from genZone to the front toe of rubble mound in terms of wave lengths'),
    ('Lend', 1.0, 'Length of domain from absZone to the back toe of rubble mound in terms of wave lengths'),
    # waves
    ('wave', True, 'Enable wave generation'),
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 1.30, "Period of the waves"),
    ("wave_height", 0.167, "Height of the waves"),
    ('wavelength', 2.121, 'Wavelength only if Fenton is activated'),
    ('Ycoeff', [0.21107604, 0.07318902, 0.02782228, 0.01234846, 0.00618291, 0.00346483, 0.00227917, 0.00194241], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [0.23112932, 0.03504843, 0.00431442, 0.00036993, 0.00004245, 0.00001877, 0.00000776, 0.00000196], 'Bcoeff only if Fenton is activated'),
    ('Nf', 8 ,'Number of frequency components for fenton waves'),
    ('meanVelocity', [ 0., 0., 0.],'Velocity used for currents'),
    ('phi0', 0.0 ,'Initial phase for waves'),
    ('Uwind', [0.0, 0.0, 0.0], 'Set air velocity'),
    ('fast', True ,'Switches ON fast cosh approximation'),
    # rubble mound
    ('porousMedia', True, 'Enable porus media region'),
    ("hs", 0.175, "Height of the breakwater"),
    ("slope1", 1./3., "Slope1 of the breakwater"),
    ("slope2", 1./2., "Slope2 of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', 0.030, "Mean diameter of the medium"),
    ('d15', None, "15% grading curve diameter of the medium"),
    ('Resistance', 'Shih', 'Ergun or Engelund or Shih'),
    # soil foundation
    ("springs", True, "Switch on/off soil module"),
    ("Kx", 541553.2, "Horizontal stiffness in Pa"),
    ("Ky", 582633.7, "Vertical stiffness in Pa"),
    ("Krot", 16246.6, "Rotational stiffness in N"),
    ("Cx", 1694.2, "Damping factor in Pa s "),
    ("Cy", 1757.32, "Damping factor in Pa s "),
    ("Crot", 69.61, "Rotational damping factor in N s "),
    # caisson
    ("caisson2D", True, "Switch on/off caisson2D"),
    ('dimx', 0.300, 'X-dimension of the caisson2D'),
    ('dimy', 0.385, 'Y-dimension of the caisson2D'),
    ('width', 1.0, 'Z-dimension of the caisson2D'),
    ('mass', 64.8/0.4, 'Mass of the caisson2D [kg]'),
    ('caissonBC', 'FreeSlip', 'caisson2D boundaries: NoSlip or FreeSlip'),
    ("rotation", False, "Initial position for free oscillation"),
    ("friction", True, "Switch on/off friction module for sliding"),
    ("overturning", True, "Switch on/off overturning module"),
    ("m_static", 0.500, "Static friction factor between caisson2D and rubble mound"),
    ("m_dynamic", 0.500, "Dynamic friction factor between caisson2D and rubble mound"),
    ('scheme', 'Runge_Kutta', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Central_Difference)'),
    # numerical options
    ("GenZone", True, 'Turn on generation zone at left side'),
    ("AbsZone", True, 'Turn on absorption zone at right side'),
    ("refinement_level", 0.0,"he=walength/refinement_level"),
    ("he", 0.01,"he=walength/refinement_level"),
    ("cfl", 0.450 ,"Target cfl"),
    ("duration", 20., "Durarion of the simulation"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', True, "Moving domain and mesh option"),
    ('conservativeFlux', True,'Fix post-processing velocity bug for porous interface'),
    ])


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
period=opts.wave_period

omega=2*np.pi/opts.wave_period

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
                                  fast = opts.fast,
                                      )

#---------Domain Dimension

nd = 2
wl = waveinput.wavelength

#---------MESH SIZE
if opts.he == 0.0:
    he = wl/opts.refinement_level
else:
    he = opts.he 


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# ----- SHAPES ----- #
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

if opts.caisson2D:
    L_leftSpo  = opts.Lgen*wl
    L_rightSpo = opts.Labs*wl

    hs=opts.hs
    slope1=opts.slope1
    slope2=opts.slope2

#-caisson2D
    dimx=opts.dimx
    dimy=opts.dimy
    b=dimx

#-Tank
    x1=L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+(hs/slope1)

    xc1=x3+0.20
    xc2=xc1+b
    yc1=yc2=hs

    x4=xc2+0.20
    x5=x4+(hs/slope2)
    x6=x5+opts.Lend*wl
    x7=x6+L_rightSpo
    tank_dim = [x7, 1.0]
 
    boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                            'x+': np.array([0., -1.,0.]),
                            'y+': np.array([0., -1.,0.]),
                            'x-': np.array([-1., 0.,0.]),
                            'sponge': None,
                            'porousLayer': None,
                            'moving_porousLayer': None,
                           }     
    boundaryTags = {'y-' : 1,
                    'x+' : 2,
                    'y+' : 3,
                    'x-' : 4,
                    'sponge' : 5,
                    'porousLayer' : 6,
                    'moving_porousLayer' : 7,
                       }

else:
    L_leftSpo  = opts.Lgen*wl
    L_rightSpo = opts.Labs*wl

    #-Tank
    x1=L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+L_rightSpo

    tank_dim = [x3, 1.0]
    
    boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                            'x+': np.array([0., -1.,0.]),
                            'y+': np.array([0., -1.,0.]),
                            'x-': np.array([-1., 0.,0.]),
                            'sponge': None,
                           }
    boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                       }


##############################################################################################################################################################################################################
# caisson2D
############################################################################################################################################################################################################

if opts.caisson2D:
    dimx=dimx
    dimy=dimy
    dim=(dimx,dimy)
    coords=[xc1+b/2., hs+dimy/2.] # For bodyDimensions and barycenter
    VCG=dim[1]/2.                 # For barycenter
    width=opts.width                     # The 3rd dimension
    mass=opts.mass #kg
    volume=float(dimx*dimy*width)
    density=float(mass/volume) #kg/m3
    I=mass*(dimx**2.+dimy**2.)/12.
    # It=(dimx**2.+dimy**2.)/12.

# --- Shape properties setup
    caisson = st.Rectangle(domain, dim=dim, coords=coords)
    caisson.vertices[0][0]=xc1
    caisson.vertices[0][1]=yc1
    caisson.vertices[1][0]=xc2
    caisson.vertices[1][1]=yc2


# --- Body properties setup
    caisson2D = bd.CaissonBody(shape=caisson, substeps=20)
    free_x=(0.0, 0.0, 0.0) # Translational DOFs
    free_r=(0.0, 0.0, 0.0) # Rotational DOFs
    m_static=opts.m_static # Static friction
    m_dynamic=opts.m_dynamic # Dynamic friction
    if opts.movingDomain==True:
        free_x=(1.0, 1.0, 0.0) # Translational DOFs
        if opts.overturning==True:
            free_r=(0.0, 0.0, 1.0) # Rotational DOFs
    caisson2D.setMass(mass)
    caisson2D.setConstraints(free_x=free_x, free_r=free_r)
    caisson2D.setFriction(friction=opts.friction, m_static=m_static, m_dynamic=m_dynamic,
                          tolerance=he/(float(10**6)), grainSize=opts.d50)
    overturning=opts.overturning
    caisson2D.setOverturning(overturning)
    if opts.rotation==True: # Initial position for free oscillation
        caisson2D.rotate(rotation)
    caisson2D.It= I/caisson2D.mass/width
    caisson2D.setNumericalScheme(scheme=opts.scheme)
    caisson2D.setRecordValues(filename='caisson2D', all_values=True)


##############################################################################################################################################################################################################
# Tank
#########################################################################################################################################################################################################

if opts.caisson2D==False:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2, 0.0], #2
              [x3,  0.0 ],#3
              [x3,  tank_dim[1] ],#4
              [x2,  tank_dim[1] ],#5
              [x1,  tank_dim[1] ],#6
              [0.0,  tank_dim[1] ],#7 
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
                           5, 5,
                          ])

    regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.90*tank_dim[1] ],
            [ 0.95*x3 , 0.95*tank_dim[1] ] ]

    regionFlags=np.array([1, 2, 3])

else:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2, 0.0], #2
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
                          6, 6,
                          1, 1, 1,
                          3, 3, 3, 3,
                          7, 7,
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

              [2,5],
              [1,10],
              [6,9],
              [3,12],
              [13,4],
             ]

    segmentFlags=np.array([1, 1,
                           6, 6,
                           1, 1,
                           2, 3, 3, 3, 4,
                           1,
                           5, 5,
                           7, 7,
                          ])

    regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.90*tank_dim[1] ],
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
dragAlpha=(porosity**2)*Alpha/nu_0
dragBeta=(porosity**3)*Beta/nu_0

#----- Spring setup

springs=opts.springs
Kx = opts.Kx
Ky = opts.Ky
Krot = opts.Krot
Cx = opts.Cx
Cy = opts.Cy
Crot = opts.Crot

if opts.caisson2D:
    caisson2D.setSprings(springs, Kx, Ky, Krot, Cx, Cy, Crot)


#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

if opts.caisson2D:
    # Caisson boundaries
    for bc in caisson.BC_list:
        if opts.caissonBC == 'FreeSlip':
            bc.setFreeSlip()
        if opts.caissonBC == 'NoSlip':
            bc.setNoSlip()

# Tank Boundaries
tank.BC['y+'].setAtmosphere()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, smoothing=3.0*he)
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

if opts.caisson2D:
    # Porous media buondaries
    tank.BC['porousLayer'].reset()
    tank.BC['moving_porousLayer'].reset()

    # Moving Mesh Options
    if opts.movingDomain==True:
        for tb in [tank.BC['x+'], tank.BC['x-'], tank.BC['y+'], tank.BC['y-'], tank.BC['sponge'], tank.BC['porousLayer']]:
            tb.hx_dirichlet.uOfXT= lambda x, t: 0.0
            tb.hy_dirichlet.uOfXT= lambda x, t: 0.0
            tb.hz_dirichlet.uOfXT= lambda x, t: 0.0
            tb.u_stress.uOfXT=None
            tb.v_stress.uOfXT=None
            tb.w_stress.uOfXT=None
        ms=tank.BC['moving_porousLayer']
        ms.hx_dirichlet.uOfXT= None
        ms.hy_dirichlet.uOfXT= None
        ms.hz_dirichlet.uOfXT= lambda x, t: 0.0
        ms.u_stress.uOfXT=None
        ms.v_stress.uOfXT=None
        ms.w_stress.uOfXT=None



########################################################################################################################################################################################################################################################################################################################################################
# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #
########################################################################################################################################################################################################################################################################################################################################################


# Waves and Generation zone
if opts.GenZone and opts.wave:
    tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, smoothing=3.0*he, dragAlpha=10.*omega/nu_0)

# Only Generation zone                        
elif opts.GenZone:
    tank.setAbsorptionZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        dragAlpha=10.*omega/nu_0)

# Porous zone
if opts.porousMedia:
    tank.setPorousZones(flags=3,
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,)

# Absorption zone
if opts.AbsZone:
    if opts.caisson2D:
        tank.setAbsorptionZones(flags=4, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(tank_dim[0]-L_rightSpo/2.), 0., 0.),
                        dragAlpha=10.*omega/nu_0)
    else:
        tank.setAbsorptionZones(flags=3, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(tank_dim[0]-L_rightSpo/2.), 0., 0.),
                        dragAlpha=10.*omega/nu_0)       



############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################
T = opts.duration

gauge_dx=0.25
tank_dim_x=int(tank_dim[0])
nprobes=int(tank_dim_x/gauge_dx)+1
probes=np.linspace(0., tank_dim_x, nprobes)
PG=[]
if opts.caisson2D:
    zProbes=hs*0.5
else:
    zProbes=opts.water_level*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

if opts.caisson2D:
    gauge_dy=0.01
    tol=np.array([1*(10**-5),1*(10**-5),0.])
    i_point_f=np.array([caisson.vertices[0][0],caisson.vertices[0][1],0.])
    i_point_f += -tol #to avoid floating point error
    i_point_b=np.array([caisson.vertices[1][0],caisson.vertices[1][1],0.])
    i_point_b += tol #to avoid floating point error
    yProbes = np.linspace(i_point_f[1],i_point_f[1]+dimy, int(dimy/gauge_dy)+1)
    LG1=[]
    LG2=[]
    for j in yProbes:
        LG1.append((i_point_f[0],j,0.),)
        LG2.append((i_point_b[0],j,0.),)

#point_output=ga.PointGauges(gauges=((('p'),PG),
#                                 ),
#                          activeTime = (0., T),
#                          sampleRate=0.,
#                         fileName='point_gauges.csv')

#loadingsGauges=ga.PointGauges(gauges=((('p'),LG1),
#                                      (('p'),LG2),
#                                 ),
#                          activeTime = (0., T),
#                          sampleRate=0.,
#                          fileName='loadingsGauges.csv')

levelset_output=ga.PointGauges(gauges=((('phi',),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='levelset_gauges.csv')


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
dt_fixed = 1
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

