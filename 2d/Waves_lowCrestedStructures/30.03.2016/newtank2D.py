from math import *
import numpy as np

from proteus import Domain
from proteus.mprans import SpatialTools as st

from proteus import MeshTools, AuxiliaryVariables

from proteus import Gauges as ga
from proteus import WaveTools as wt

from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
periods=[1., 1.5, 2.]
period=periods[0]

waterLevels=[0.315, 0.276]
waterLevel = waterLevels[0]

waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface

waveHeights=[0.03, 0.04, 0.05,
             0.06, 0.07, 0.08,
             0.09, 0.10, 0.11, 0.12]
waveHeight=waveHeights[0]

omega = float(2.0*pi/period)

wavelengths=[1.39, 1.34] #from dispersion law
wavelength = wavelengths[0]

inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])
k = float(2.0*pi/wavelength)

#---------Domain Dimension
nd = 2


# ----- SHAPES ----- #

tank_dim = [16.2, 0.75]
L_leftSpo  = 1.4 # 1wavelength
L_rightSpo = 2.8 # 2wavelength
hs=0.25
b=0.25
x2=6.5
x3=x2+hs*2
x4=x3+b
x5=x4+hs*1.5
x6=tank_dim[0]-L_rightSpo


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


vertices=[[0., 0.],#0
          [L_leftSpo, 0.],#1
          [x2, 0.],#2
          [x3, hs],#3
          [x4, hs],#4
          [x5, 0.],#5
          [tank_dim[0]-L_rightSpo, 0.],#6
          [tank_dim[0], 0.],#7
          [tank_dim[0], tank_dim[1]],#8
          [tank_dim[0]-L_rightSpo, tank_dim[1]],#9
          [L_leftSpo, tank_dim[1]],#10
          [0., tank_dim[1]],#11
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
          [2,5],
          [6,9],
         ]

segmentFlags=np.array([1, 1, 5, 5, 5, 1, 1,
                       2, 3, 3, 3, 4,
                       5, 1, 5,
                      ])


regions = [ [ 0.10*tank_dim[0] , 0.1*tank_dim[1] ],
	    [ 0.50*tank_dim[0] , 0.95*tank_dim[1] ],
            [ 0.99*x4 , 0.99*hs ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ],
          ]

regionFlags=np.array([1, 2, 3, 4,
                     ])


tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


##########################################
# Numerical Options and other parameters #
##########################################

nd = 2

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5

# Surface tension
sigma_01 = 0.0

# Gravity
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))

##########################################
# POROUS MEDIA
##########################################

porosity=0.4
voidFrac=1.0-porosity
d50=0.058
d15=d50/1.2

Alpha1=1684+3.12*(10**-3)*((gAbs/(nu_0**2))**(2/3))*(d15**2) #Shih
#Alpha1=150 #Ergun
Beta1=1.72+1.57*exp(-5.10*(10**-3)*((gAbs/(nu_0**2))**(1/3))*d15) #Shih
#Beta1=1.75 #Ergun

#Alpha=Alpha1*nu_0*(voidFrac**3)/((porosity**2)*(d15**2))  #Engelund
Alpha=Alpha1*nu_0*(voidFrac**2)/((porosity**3)*(d15**2))  #Ergun
Beta=Beta1*voidFrac/((porosity**3)*d15)

#Proteus scale in viscosity, so i need to divide alpha and beta by nu_0
dragAlpha=Alpha/nu_0
dragBeta=0.#Beta/nu_0


#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------

weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = period/20
dt_init = min(0.1*dt_fixed,0.001)
T = 50*period
nDTout= int(round(T/dt_fixed))
runCFL = 0.3 #0.9




# ----- WAVE input ----- #

"""
waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=wavelength,
                                  waveType='Fenton',
                                  Ycoeff = Ycoeff,
                                  Bcoeff = Bcoeff,
                                  )
"""
waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=wavelength,
                                  waveType='Linear',
                                 )



# ----- BOUNDARY CONDITIONS ----- #

tank.BC.top.setOpenAir()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)
#tank.BC.left.setFreeSlip()
tank.BC.bottom.setFreeSlip()
tank.BC.bottom.vof_advective=0.
tank.BC.right.setNoSlip()
tank.BC.right.vof_advective=0.
#tank.BC.sponge.setParallelFlag0()
tank.BC.sponge.setNonMaterial()


# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #

tank.setGenerationZones(flags=1, epsFact_solid=L_leftSpo/2, orientation=[1, 0, 0], center=(L_leftSpo/2, 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setPorousZones(flags=3, epsFact_solid=(x5-x2)/2,
                    dragAlpha=dragAlpha, dragBeta=dragBeta,
                    porosity=porosity,
                   )
tank.setAbsorptionZones(flags=4, epsFact_solid=L_rightSpo/2, orientation=[-1, 0, 0,],
                        center=(tank_dim[0]-L_rightSpo/2, 0., 0.))

#
# ----- Output Gauges ----- #

gauge_dx=1.0 #0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)#
LG=[]
PG=[]
PG0=[]
for i in probes:
    #LG.append(((i,0.,0.), (i,tank_dim[1],0.)),)
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


#fields=(('vof',))
#gauges_locations=(((0., 0., 0.), (0., tank_dim[1], 0.)),
#                  ((23., 0., 0.), (23., tank_dim[1], 0.)))
#integral_output=ga.LineIntegralGauges(gauges=((fields, gauges_locations), #,LG),
 #                                             ),
  #                                    activeTime = (0., T),
   #                                   sampleRate=0.,
    #                                  fileName='line_integral_gauges.csv')

domain.auxiliaryVariables += [point_output, #integral_output, #line_output
                              ]


# ----- Mesh ----- #

he=wavelength/100
domain.MeshOptions.elementSize(he)
st.assembleDomain(domain)

domain.writePoly("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

quad_order = 3


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------


checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True #False



#  Discretization -- input options


useOnlyVF = False
movingDomain=False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988

genMesh=True
useOldPETSc=False
useSuperlu = False
#timeDiscretization='be' #'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0

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
ns_forceStrongDirichlet = False #True
backgroundDiffusionFactor=0.01


""" Default configuration by Tristan
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True #False
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
"""


if useMetrics:
    ns_shockCapturingFactor  = 0.25
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.25
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.25
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.25
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
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
rd_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he)
mcorr_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)




#turbulence
ns_closure=2 #2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4

# Initial condition
waterLine_x = 2*tank_dim[0]
waterLine_z = inflowHeightMean



def waveHeight(x,t):
   # waterDepth=waveinput.waterDepth(x[0], x[1], x[2], t)
    waterDepth=waveinput.eta(x, t)+waveinput.mwl
    return waterDepth


def wavePhi(x,t):
    if nd==2:
        return x[1]- waveHeight(x,t)
    else:
        return x[2]- waveHeight(x,t)



def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def signedDistance(x):
    phi_x = x[0]-waterLine_x # (tank_dim[0]-L_rightSpo) # x0_por
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

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))


