from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
period=1.94
waterLevel = 1.
waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface
waveHeight=0.025
omega = float(2.0*pi/period)
wavelength = 4.997 #from dispersion law
inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])
k = float(2.0*pi/wavelength)

#---------Domain Dimension
nd = 2


# ----- SHAPES ----- #

tank_dim = [20.0, 1.5]

L_leftSpo  = 5.# wavelength
L_rightSpo  = 5.# wavelength


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


vertices=[[0.0,0.0],#0
          [L_leftSpo,0.0],#1
          [tank_dim[0]-L_rightSpo,0.0],#2
          [tank_dim[0],0.0],#3
          [tank_dim[0],tank_dim[1]],#4
          [tank_dim[0]-L_rightSpo,tank_dim[1]],#5
          [L_leftSpo,tank_dim[1]],#6
          [0.0,tank_dim[1]],#7
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

segmentFlags=np.array([1, 1, 1, 2,
                       3, 3, 3, 4,
                       5, 5,
                      ])


regions = [ [ 0.1*tank_dim[0] , 0.1*tank_dim[1] ],
            [ 0.5*tank_dim[0], 0.95*tank_dim[1] ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2, 3])



tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


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
                                  wavelength=wavelength,
                                  waveType="Linear")


# ----- BOUNDARY CONDITIONS ----- #

for bc in tank.BC_list:
    bc.setFreeSlip()

tank.BC.top.setOpenAir()

tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)

tank.BC.bottom.setNoSlip()

tank.BC.right.setFreeSlip()

tank.BC.sponge.setNonMaterial()



# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #

dragAlpha=1.96*(10**6)
dragBeta=0.
porosity=0.4


tank.setGenerationZones(flags=1, epsFact_solid=float(L_leftSpo/2.),
                        orientation=[1., 0., 0.], center=(float(L_leftSpo/2.), 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )
tank.setAbsorptionZones(flags=3, epsFact_solid=float(L_rightSpo/2.),
                       orientation=[-1., 0., 0.], center=(float(L_rightSpo/2.), 0., 0.),
                       )

#
# ----- Output Gauges ----- #
T = 40*period
gauge_dx=1.0 #0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)
PG0=[]
for i in probes:
    PG0.append((i, 0., 0.),)


point_output=ga.PointGauges(gauges=((('p'),PG0),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='point_gauges.csv')

domain.auxiliaryVariables += [point_output]



##########################################
# Numerical Options and other parameters #
##########################################

he = wavelength/400 #H/he should be not so small beacause of possible dissipations in genZone!!!
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
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=False
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = period/19.4
dt_init = min(0.1*dt_fixed,0.001)
T = T
nDTout= int(round(T/dt_fixed))
runCFL = 0.9

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 0.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
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
    redist_Newton = True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True
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
    redist_Newton = False#True
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
