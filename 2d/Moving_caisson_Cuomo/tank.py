from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd


opts=Context.Options([
    # predefined test cases
    ("water_level", 1.044, "Height of free surface above bottom"),
    # Geometry
    ('Lgen', 1.0, 'Genaration zone in terms of wave lengths'),
    ('Labs', 2.0, 'Absorption zone in terms of wave lengths'),
    ('Ls', 2.0, 'Length of domain from genZone to the front toe of rubble mound in terms of wave lengths'),
    ('Lend', 1.0, 'Length of domain from absZone to the back toe of rubble mound in terms of wave lengths'),
    ('th', 3.0, 'Total height of numerical tank'),
    # waves
    ('waveType', 'Fenton', 'Wavetype for regular waves, Linear or Fenton'),
    ("wave_period", 3.04, "Period of the waves"),
    ("wave_height", 0.428, "Height of the waves"), # [0.428, 0.524, 0.569, 0.619]
    ('wavelength', 9.328, 'Wavelength only if Fenton is activated'), # [9.328, 9.485, 9.568, 9.666]
    ('Ycoeff', [ 0.12585567    ,0.04571332    ,0.01554661    ,0.00565350    ,0.00222585,0.00093237    ,0.00040825    ,0.00018459    ,0.00008553   ,0.00004041,0.00001943   ,0.00000955   ,0.00000490   ,0.00000286   ,0.00000228], 'Ycoeff only if Fenton is activated'),
    ('Bcoeff', [0.15535477    ,0.02753642    ,0.00441347    ,0.00049197    ,0.00000772,-0.00001065    ,-0.00000195    ,0.00000017    ,0.00000018   ,0.00000005,0.00000001   ,0.00000000   ,0.00000000   ,0.00000000   ,0.00000000], 'Bcoeff only if Fenton is activated'),
    # rubble mound
    ("hs1", 0.644, "Height of the breakwater"),
    ("hs2", 0.244, "Height of the breakwater"),
    ("slope1", 1./1.5, "Slope1 of the breakwater"),
    ("slope2", 1./2., "Slope2 of the breakwater"),
    ('porosity', 0.4, "Porosity of the medium"),
    ('d50', 0.020, "Mean diameter of the medium"),
    ('d15', None, "15% grading curve diameter of the medium"),
    ('Resistance', 'Shih', 'Ergun or Engelund or Shih'),
    # soil foundation
    ("springs", True, "Switch on/off soil module"),
    ("Kx", 1.5*(10**8)/0.78, "Horizontal stiffness in Pa"),
    ("Ky", 1.4*(10**8)/0.78, "Vertical stiffness in Pa"),
    ("Krot", 4.91*(10**7)/0.78, "Rotational stiffness in N"),
    ("Cx", 4.8*(10**5)/0.78, "Damping factor in Pa s "),
    ("Cy", 10.*4.8*(10**5)/0.78, "Damping factor in Pa s "),
    ("Crot", 10.*4.8*(10**5)/0.78, "Rotational damping factor in N s "),
    # caisson
    ("caisson2D", True, "Switch on/off caisson2D"),
    ('dimx', 1.0, 'X-dimension of the caisson2D'),
    ('dimy', 1.12, 'Y-dimension of the caisson2D'),
    ('width', 1.0, 'Z-dimension of the caisson2D'),
    ('mass', 2117.15, 'Mass of the caisson2D [kg]'),
    ('caissonBC', 'FreeSlip', 'caisson2D boundaries: NoSlip or FreeSlip'),
    ("rotation", False, "Initial position for free oscillation"),
    ("friction", True, "Switch on/off friction module for sliding"),
    ("overturning", True, "Switch on/off overturning module"),
    ("m_static", 0.600, "Static friction factor between caisson2D and rubble mound"),
    ("m_dynamic", 0.400, "Dynamic friction factor between caisson2D and rubble mound"),
    ('scheme', 'Runge_Kutta', 'Numerical scheme applied to solve motion calculation (Runge_Kutta or Central_Difference)'),
    # numerical options
    ("refinement_level", 200.,"he=walength/refinement_level"),
    ("he", 0.00,"Mesh size"),
    ("cfl", 0.30 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', True, "Moving domain and mesh option"),
    ('conservativeFlux', False,'Fix post-processing velocity bug for porous interface'),
    ])

# Fenton 
# Y coeff
#   0.12585567    ,0.04571332    ,0.01554661    ,0.00565350    ,0.00222585,0.00093237    ,0.00040825    ,0.00018459    ,0.00008553   ,0.00004041,0.00001943   ,0.00000955   ,0.00000490   ,0.00000286   ,0.00000228
#   0.14172859    ,0.06063868    ,0.02467932    ,0.01080073    ,0.00513114 ,0.00259759    ,0.00137627    ,0.00075390    ,0.00042382   ,0.00024361,0.00014329   ,0.00008704   ,0.00005609   ,0.00004063   ,0.00003594
#   0.14686514    ,0.06695851    ,0.02928978    ,0.01382733    ,0.00709924 ,0.00388864    ,0.00223155    ,0.00132551    ,0.00080928   ,0.00050648,0.00032568   ,0.00021742   ,0.00015435   ,0.00012147   ,0.00011127
#   0.15047047    ,0.07299799    ,0.03434850    ,0.01753599    ,0.00976622 ,0.00581525    ,0.00363496    ,0.00235717    ,0.00157582   ,0.00108422,0.00077024   ,0.00057031   ,0.00044763   ,0.00038120   ,0.00036018

# B coeff
#   0.15535477    ,0.02753642    ,0.00441347    ,0.00049197    ,0.00000772,-0.00001065    ,-0.00000195    ,0.00000017    ,0.00000018   ,0.00000005,0.00000001   ,0.00000000   ,0.00000000   ,0.00000000   ,0.00000000
#   0.17300514    ,0.03615329    ,0.00718074    ,0.00112536    ,0.00010505 ,-0.00000271    ,-0.00000101    ,0.00000133    ,0.00000090   ,0.00000032,0.00000008   ,0.00000001   ,0.00000000   ,0.00000000   ,0.00000000
#   0.17828208    ,0.03971863    ,0.00862872    ,0.00156338    ,0.00021009 ,0.00001958    ,0.00000488    ,0.00000380    ,0.00000207   ,0.00000081,0.00000025   ,0.00000007   ,0.00000002   ,0.00000001   ,0.00000000
#   0.18152464    ,0.04305215    ,0.01026377    ,0.00216538    ,0.00039704 ,0.00007496    ,0.00002349    ,0.00001149    ,0.00000557   ,0.00000237,0.00000091   ,0.00000034   ,0.00000013   ,0.00000006   ,0.00000002

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

if opts.waveType=='Fenton':
    waveinput = wt.MonochromaticWaves(period=period,
                                      waveHeight=waveHeight,
                                      mwl=mwl,
                                      depth=waterLevel,
                                      g=g,
                                      waveDir=waveDir,
                                      wavelength=opts.wavelength, # if wave is linear I can use None
                                      waveType="Fenton",
                                      Ycoeff=opts.Ycoeff,
                                      Bcoeff=opts.Bcoeff,
                                      )

#---------Domain Dimension

nd = 2.
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

    hs1=opts.hs1
    hs2=opts.hs2
    slope1=opts.slope1
    slope2=opts.slope2

#-Caisson
    dimx=opts.dimx
    dimy=opts.dimy
    b=dimx

#-Tank
    x1=L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+(hs1/slope1)
    x4=x3+1.14
    x5=x4
    x6=x5+b
    x7=x6+2.*hs2
    x8=x7+(hs2/slope2)
    x9=x8+opts.Lend*wl
    x10=x9+L_rightSpo

    y3=hs1
    y5=hs2

    xc1=x5
    xc2=xc1+b
    yc1=yc2=hs2

    tank_dim = [x10, opts.th]

else:
    L_leftSpo  = opts.Lgen*wl
    L_rightSpo = opts.Labs*wl

    hs1=opts.hs1
    hs2=opts.hs2
    slope1=opts.slope1
    slope2=opts.slope2

#-Tank
    x1=L_leftSpo
    x2=x1+opts.Ls*wl
    x3=x2+(hs1/slope1)
    x4=x3+1.14
    x5=x4
    x6=x5+b
    x7=x6+2.*hs2
    x8=x7+(hs2/slope2)
    x9=x8+opts.Lend*wl
    x10=x9+L_rightSpo

    y3=hs1
    y5=hs2

    xc1=x5
    xc2=xc1+b
    yc1=yc2=hs2

    tank_dim = [x10, opts.th]


boundaryOrientations = {'y-': [0., -1.,0.],
                            'x+': [1., 0.,0.],
                            'y+': [0., 1.,0.],
                            'x-': [-1., 0.,0.],
                            'sponge': None,
                            'porousLayer': None,
                            'moving_porousLayer': None,
                           }
boundaryTags = {'y-': 1,
                    'x+': 2,
                    'y+': 3,
                    'x-': 4,
                    'sponge': 5,
                    'porousLayer': 6,
                    'moving_porousLayer': 7,
                       }


##############################################################################################################################################################################################################
# caisson2D
############################################################################################################################################################################################################

if opts.caisson2D:
    dimx=dimx
    dimy=dimy
    dim=(dimx,dimy)
    coords=[xc1+b/2., hs2+dimy/2.] # For bodyDimensions and barycenter
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

if opts.caisson2D==True:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2,  0.0],#2
              [x3,  hs1],#3
              [x4,  hs1],#4
              [x5,  hs2],#5
              [x6,  hs2],#6
              [x7,  hs2],#7
              [x8,  0.0],#8
              [x9,  0.0],#9
              [x10, 0.0],#10
              [x10,   tank_dim[1]],#11
              [x9,    tank_dim[1]],#12
              [x1,    tank_dim[1]],#13
              [0.0,   tank_dim[1]],#14
              ]

    vertexFlags=np.array([1, 1, 1,
                          6,
                          7, 7, 7, 
                          6,
                          1, 1, 1,
                          3, 3, 3, 3,
                          ])
    segments=[[0,1],
              [1,2],
              [2,3],
              [3,4],

              [6,7],
              [7,8],
              [8,9],
              [9,10],
              [10,11],
              [11,12],
              [12,13],
              [13,14],
              [14,0],

              [1,13],
              [9,12],
              [2,8],
             ]

    segmentFlags=np.array([ 1, 1,
                            6, 7, 7, 6,
                            1, 1,
                            2, 3, 3, 3, 4,
                            5, 5,
                            1,
                         ])

else:

    vertices=[[0.0, 0.0],#0
              [x1,  0.0],#1
              [x2,  0.0],#2
              [x3,  hs1],#3
              [x4,  hs1],#4
              [x5,  hs2],#5
              [x6,  hs2],#6
              [x7,  hs2],#7
              [x8,  0.0],#8
              [x9,  0.0],#9
              [x10, 0.0],#10
              [x10,   tank_dim[1]],#11
              [x9,    tank_dim[1]],#12
              [x1,    tank_dim[1]],#13
              [0.0,   tank_dim[1]],#14
              ]

    vertexFlags=np.array([1, 1, 1,
                          6,
                          7, 7, 7, 
                          6,
                          1, 1, 1,
                          3, 3, 3, 3,
                          ])
    segments=[[0,1],
              [1,2],
              [2,3],
              [3,4],

              [6,7],
              [7,8],
              [8,9],
              [9,10],
              [10,11],
              [11,12],
              [12,13],
              [13,14],
              [14,0],

              [1,13],
              [9,12],
              [2,8],

              [4,5],
              [5,6],
             ]

    segmentFlags=np.array([ 1, 1,
                            6, 6, 6, 6, 
                            1, 1,
                            2, 3, 3, 3, 4,
                            5, 5,
                            1,
                            6, 6,                            
                         ])


regions = [ [ 0.90*x1 , 0.10*tank_dim[1] ],
            [ 0.90*x2 , 0.90*tank_dim[1] ],
            [ xc1 , 0.50*hs2 ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]

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
dragBeta=0.0#(porosity**3)*Beta/nu_0

#----- Spring setup

springs=opts.springs
Kx = opts.Kx
Ky = opts.Ky
Krot = opts.Krot
Cx = opts.Cx
Cy = opts.Cy
Crot = opts.Crot

caisson2D.setSprings(springs, Kx, Ky, Krot, Cx, Cy, Crot)


#############################################################################################################################################################################################################################################################################################################################################################################################
# ----- BOUNDARY CONDITIONS ----- #
#############################################################################################################################################################################################################################################################################################################################################################################################

if opts.caisson2D:
    for bc in caisson.BC_list:
        if opts.caissonBC == 'FreeSlip':
            bc.setFreeSlip()
        if opts.caissonBC == 'NoSlip':
            bc.setNoSlip()

tank.BC['y+'].setAtmosphere()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1)
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['porousLayer'].reset()
tank.BC['moving_porousLayer'].reset()


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


tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput,
                        )
tank.setPorousZones(flags=3,
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=float(L_rightSpo/2.),
                        orientation=[-1., 0.], center=(float(tank_dim[0]-L_rightSpo/2.), 0., 0.),
                        )


############################################################################################################################################################################
# ----- Output Gauges ----- #
############################################################################################################################################################################

T = 30.*period

gauge_dx=0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)
PG=[]
zProbes=hs2*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

gauge_dx=0.01
gauge_dy=0.01
tol=np.array([1*(10**-5),1*(10**-5),0.])

# front
i_point_f = np.array([caisson.vertices[0][0],caisson.vertices[0][1],0.])
i_point_f += -tol #to avoid floating point error
# bottom
i_point_bo = np.array([caisson.vertices[0][0],caisson.vertices[0][1],0.])
i_point_bo += -tol #to avoid floating point error
# back
i_point_ba = np.array([caisson.vertices[1][0],caisson.vertices[1][1],0.])
i_point_ba += tol #to avoid floating point error
# top
i_point_t = np.array([caisson.vertices[3][0],caisson.vertices[3][1],0.])
i_point_t += tol #to avoid floating point error

xProbes = np.linspace(i_point_t[0],i_point_t[0]+dimx, (dimx/gauge_dx)+1.)
yProbes = np.linspace(i_point_f[1],i_point_f[1]+dimy, (dimy/gauge_dy)+1.)

LG1=[]
LG2=[]
LG3=[]
LG4=[]

for ii in xProbes:
    LG2.append((ii,i_point_bo[1],0.),)
    LG4.append((ii,i_point_t[1],0.),)
for jj in yProbes:
    LG1.append((i_point_f[0],jj,0.),)
    LG3.append((i_point_ba[0],jj,0.),)
   

point_output=ga.PointGauges(gauges=((('p'),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='point_gauges.csv')

loadingsGauges=ga.PointGauges(gauges=((('p'),LG1),
                                      (('p'),LG2),
                                      (('p'),LG3),
                                      (('p'),LG4),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='loadingsGauges.csv')


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

