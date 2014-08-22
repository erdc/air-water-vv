from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
import numpy as np
import matplotlib.pyplot as plt

transect = np.loadtxt('transectBox.txt',skiprows=0,delimiter="\t")
#for  debugging
#plt.plot(transect[:,0],transect[:,1],'g-o')
#plt.title("transect")
#plt.show()

#wave generator
t_h_u = np.loadtxt('t_h_u.txt')
from scipy.interpolate import InterpolatedUnivariateSpline
surface_displacement_spline = InterpolatedUnivariateSpline(t_h_u[:,0],t_h_u[:,1])
velocity_spline = InterpolatedUnivariateSpline(t_h_u[:,0],t_h_u[:,2])
#plt.plot(t_h_u[:100,0],surface_displacement_spline(t_h_u[:100,0]),'b',t_h_u[:100,0],velocity_spline(t_h_u[:100,0]),'r')
#plt.title("wave generator input")
#plt.show()

windVelocity = (0.0,0.0)
inflowHeightMean = 0.60959#m
inflowVelocityMean = (0.0,0.0)
amplitude = 0.5*(max(t_h_u[:,1])-min(t_h_u[:,1]))
print amplitude
#  Discretization -- input options  
genMesh=True
useOldPETSc=False
useSuperlu=False
timeDiscretization='vbdf'#'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
applyCorrection=True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
useMetrics=1.0
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()    
    
#  Discretization   
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3) 	    
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
    
# Domain and mesh
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
maxZ = max(transect[:,1])
minZ = min(transect[:,1])
maxX = max(transect[:,0])
minX = min(transect[:,0])
Lx = maxX - minX
Lz = max(inflowHeightMean+7.0*amplitude,maxZ)
L = (Lx,Lz,1.0)
he = Lz/float(21)#amplitude*0.1
boundaries=['left','right','bottom','top','front','back','obstacle']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
vertices=[[transect[-1,0],inflowHeightMean-5.0*amplitude],#0
          [transect[-1,0],inflowHeightMean+5.0*amplitude],#1
          [transect[-1,0],minZ+Lz],#2
          [transect[-2][0],minZ+Lz],#3
          [transect[0,0],minZ+Lz],#4
          [transect[0,0],inflowHeightMean+5.0*amplitude],#5
          [transect[0,0],inflowHeightMean-5.0*amplitude]]#6
vertexFlags=[boundaryTags['right'],#0
             boundaryTags['right'],#1
             boundaryTags['right'],#2
             boundaryTags['top'],#3
             boundaryTags['top'],#4
             boundaryTags['left'],#5
             boundaryTags['left']]#6
for p in transect:
    vertices.append([p[0],p[1]])
    vertexFlags.append(boundaryTags['bottom'])
nSegments = len(vertices)
segments=[]
segmentFlags=[]
for i in range(nSegments):
    segments.append([i,(i+1) % nSegments])
    if i in [0,1,nSegments-1]:
        segmentFlags.append(boundaryTags['right'])
    elif i in [2,3]:
        segmentFlags.append(boundaryTags['top'])
    elif i in [4,5,6]:
        segmentFlags.append(boundaryTags['left'])
    else:
        segmentFlags.append(boundaryTags['bottom'])
domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags)
#go ahead and add a boundary tags member 
domain.boundaryTags = boundaryTags
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
# Time stepping
T=t_h_u[-1,0]
dt_fixed =t_h_u[1,0]
dt_init = min(0.01*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False
ns_shockCapturingFactor  = 0.9
ns_lag_shockCapturing = True
ns_lag_subgridError = True
ls_shockCapturingFactor  = 0.9
ls_lag_shockCapturing = True
ls_sc_uref  = 1.0
ls_sc_beta  = 1.5
vof_shockCapturingFactor = 0.9
vof_lag_shockCapturing = True
vof_sc_uref = 1.0
vof_sc_beta = 1.5
rd_shockCapturingFactor  = 0.9
rd_lag_shockCapturing = False
epsFact_density    = 1.5
epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
epsFact_redistance = 0.33
epsFact_consrv_diffusion = 10.0
redist_Newton = True
kappa_shockCapturingFactor = 0.9
kappa_lag_shockCapturing = True#False
kappa_sc_uref = 1.0
kappa_sc_beta = 1.0
dissipation_shockCapturingFactor = 0.1
dissipation_lag_shockCapturing = True#False
dissipation_sc_uref = 1.0
dissipation_sc_beta = 1.0

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
rd_nl_atol_res = max(1.0e-12,0.1*he)
mcorr_nl_atol_res = max(1.0e-12,0.001*he**2)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,-9.8]

# Initial condition
waterLine_x = 2*L[0]
waterLine_z = inflowHeightMean

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[1]-waterLine_z 
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

def waveHeight(x,t):
    return inflowHeightMean + surface_displacement_spline(t)

def waveVelocity_u(x,t):
    return velocity_spline(t)

def waveVelocity_v(x,t):
    return 0.0

#solution variables

def wavePhi(x,t):
    return x[1] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = 0.0
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)
