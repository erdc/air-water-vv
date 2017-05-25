from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
from proteus import Comm
from proteus import Context
from proteus.WaveTools import TimeSeries
from proteus.Domain import InterpolatedBathymetryDomain, PiecewiseLinearComplexDomain
from proteus.MeshTools import InterpolatedBathymetryMesh

comm = Comm.init()
opts=Context.Options([
    ("wave_type", 'single-peaked', "type of waves generated: 'linear', 'Nonlinear', 'single-peaked', 'double-peaked', 'time-series'"),
    ("depth", 7.25, "water depth [m]"),
    ("wave_height", 4.0, "wave height [m]"),
    ("peak_period", 1.0/0.09, "Peak period [s]"),
    ("peak_period2", 6.0, "Second peak period (only used in double-peaked case)[s]"),
    ("peak_wavelength",10.0,"Peak wavelength in [m]"),
    ("parallel", True, "Run in parallel"),
    ("gauges", False, "Enable gauges")])

# Wave generator
windVelocity = [0., 0., 0.]
depth = 7.0
inflowVelocityMean = [0., 0., 0.]
period = opts.peak_period
omega = 2.0*math.pi/period
waveheight = opts.wave_height
amplitude = waveheight/2.0
wavelength = opts.peak_wavelength
k = -2.0*math.pi/wavelength


# Discretization -- input options
genMesh = True #False
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False #not opts.parallel
timeDiscretization = 'be'  #'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega


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


# Discretization
nd = 3
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
he = 0.2  #wavelength/50.0     

nLevels = 1
weak_bc_penalty_constant = 100.0


#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

structured=False

if genMesh:
    if useHex:
        nnx=4*Refinement+1
        nny=2*Refinement+1
        hex=True
        domain = Domain.RectangularDomain(L)
    else:
        bathy=np.genfromtxt("FRF_FRF_20150915_1116_NAVD88_LARC_GPS_UTC.csv", delimiter=",",skiprows=1)
        xmax=bathy[:,7].max()
        xmin=bathy[:,7].min()
        ymax=bathy[:,8].max()
        ymin=bathy[:,8].min()
        zmax=bathy[:,9].max()
        zmin=bathy[:,9].min()
        
        #reset to domain of interest
        xmin = 60.0 #425.0
        xmax = 450.0
        ymin = 900.0
        ymax = 900.0 + 5.0*he
        #
        boundaries=['empty','left','right','bottom','top','front','back']
        boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
        #
        #process 2D domain
        #
        bathy_points = np.vstack((bathy[:,7],bathy[:,8],bathy[:,9])).transpose()
        domain2D = InterpolatedBathymetryDomain(vertices=[[xmin,ymin],[xmin,ymax],[xmax,ymax],[xmax,ymin]],
                                              vertexFlags=[boundaryTags['left'],boundaryTags['left'],boundaryTags['right'],boundaryTags['right']],
                                              segments=[[0,1],[1,2],[2,3],[3,0]],
                                              segmentFlags=[boundaryTags['left'],boundaryTags['back'],boundaryTags['right'],boundaryTags['front']],
                                              regions=[(0.5*(xmax+xmin),0.5*(ymax+ymin))],
                                              regionFlags=[1],
                                              name="frfDomain2D",
                                              units='m',
                                              bathy = bathy_points)
        domain2D.writePoly(domain2D.name)

        #now adaptively refine 2D mesh to interpolate bathy to desired accuracy
        mesh2D = InterpolatedBathymetryMesh(domain2D,
                                            triangleOptions="gVApq30Dena%8.8f" % ((1000.0**2)/2.0,),
                                            atol=1.0e-1,
                                            rtol=1.0e-1,
                                            maxLevels=25,
                                            maxNodes=50000,
                                            bathyType="points",
                                            bathyAssignmentScheme="interpolation",
                                            errorNormType="Linfty")
        #
        # Done processing 2D
        #
        fineMesh = mesh2D.meshList[-1]
        newNodes = {}
        verticalEdges = {}
        nN_start = fineMesh.nodeArray.shape[0]
        nN = nN_start
        zTop = max(fineMesh.nodeArray[:,2].max(),fineMesh.nodeArray[:,2].min()+depth+2*waveheight)
        for nN_bottom, n,f in zip(range(nN_start), fineMesh.nodeArray, fineMesh.nodeMaterialTypes):
            if f > 0:
                newNodes[nN] = (n[0],n[1],zTop)
                verticalEdges[nN_bottom] = nN
                nN += 1
        newFacets = []
        newFacetFlags = []
        for t in fineMesh.elementNodesArray:
            newFacets.append([[t[0],t[1],t[2]]])
            newFacetFlags.append(6)
        topConnectivity = {}
        for edge, edgeF in zip(fineMesh.elementBoundaryNodesArray, fineMesh.elementBoundaryMaterialTypes):
            if edgeF > 0:
                n00 = edge[0]
                n10 = edge[1]
                n11 = verticalEdges[edge[1]]
                n01 = verticalEdges[edge[0]]
                newFacets.append([[n00,
                                  n10,
                                  n11,
                                  n01]])
                newFacetFlags.append(edgeF)
                if topConnectivity.has_key(n11):
                    topConnectivity[n11].append(n01)
                else:
                    topConnectivity[n11] = [n01]
                if topConnectivity.has_key(n01):
                    topConnectivity[n01].append(n11)
                else:
                    topConnectivity[n01] = [n11]
        topFacet = [topConnectivity[nN_start][0],nN_start,topConnectivity[nN_start][1]]

        while len(topFacet) < len(newNodes):
            nN = topFacet[-1]
            if topFacet[-2] == topConnectivity[nN][0]:
                topFacet.append(topConnectivity[nN][1])
            else:
                topFacet.append(topConnectivity[nN][0])
        newFacets.append([topFacet])
        newFacetFlags.append(5)
        newVertices = []
        newVertexFlags = []
        for n,nF in zip(fineMesh.nodeArray, fineMesh.nodeMaterialTypes):
            newVertices.append([n[0],n[1],n[2]])
            if nF > 0:
                newVertexFlags.append(nF)
            else:
                newVertexFlags.append(6)
        for nN in range(len(newNodes.keys())):
            newVertices.append(newNodes[nN+nN_start])
            newVertexFlags.append(5)
        xmin_new = fineMesh.nodeArray[:,0].min()
        xmax_new = fineMesh.nodeArray[:,0].max()
        ymin_new = fineMesh.nodeArray[:,1].min()
        ymax_new = fineMesh.nodeArray[:,1].max()

        regions = [[0.5*(xmin+xmax),0.5*(ymin+ymax),zTop-0.001]]
        regionFlags = [1.0]
        domain = PiecewiseLinearComplexDomain(vertices=newVertices,
                                              facets=newFacets,
                                              facetHoles=None,
                                              holes=None,
                                              regions=regions,
                                              vertexFlags=newVertexFlags,
                                              facetFlags=newFacetFlags,
                                              regionFlags=regionFlags,
                                              regionConstraints=None,
                                              name="frfDomain3D",
                                              units="m")
        domain.writePoly("frfDomain3D")
        domain.writePLY("frfDomain3D")
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)
        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
        porosityTypes = numpy.array([1.0])
        dragAlphaTypes = numpy.array([0.0])
        dragBetaTypes = numpy.array([0.0])
        epsFact_solidTypes = np.array([0.0])
        comm.barrier()
else:
    boundaries=['empty','left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    from proteus.Domain import  PiecewiseLinearComplexDomain
    domain = PiecewiseLinearComplexDomain(fileprefix="frfDomain3D")
    domain.boundaryTags = boundaryTags


zmin = np.array(domain.vertices)[:,2].min()
inflowHeightMean = zmin + depth


# Time stepping
T = 70.
dt_fixed = 1.
dt_init = 0.001
runCFL=0.33
nDTout = int(round(T/dt_fixed))


# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.01
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
    epsFact_consrv_diffusion = 10.
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
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0


ns_nl_atol_res = max(1.0e-10,0.001*he**2)
vof_nl_atol_res = max(1.0e-10,0.001*he**2)
ls_nl_atol_res = max(1.0e-10,0.001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)


# Turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4


# Physical Parameters
rho_0 = 998.2
nu_0  = 1.004e-6
rho_1 = 1.205
nu_1  = 1.500e-5
sigma_01 = 0.
g = np.array([0., 0., -9.8])


# Initial condition
waterLine_z =  inflowHeightMean


def signedDistance(x):
    phi_z = x[2]-waterLine_z
    return phi_z


def theta(x,t):
    return k*x[0] - omega*t + math.pi/2.0


def z(x):
    return x[2] - inflowHeightMean


def ramp(t):
  t0=.1 #ramptime
  if t<t0:
    return t/t0
  else:
    return 1


domain_vertices = numpy.array(domain.vertices)
h = inflowHeightMean - domain_vertices[:,2].min()# - transect[0][1] if lower left hand corner is not at z=0
sigma = omega - k*inflowVelocityMean[0]


timeSeriesFile = "Duck_series.txt"
skiprows = 0
N = 32
mwl = inflowHeightMean
waveDir = np.array([-1,0,0])
rec_direct = True
window_params = None
timeSeriesPosition = [0., 0., 0.]

tseries = TimeSeries(timeSeriesFile,
                     skiprows,
                     timeSeriesPosition,
                     depth,     #Need to set the depth
                     N,         #Dummy
                     mwl,       #mean water level
                     waveDir,
                     g,         #acceleration of gravity
                     rec_direct,
                     window_params,
                     )


def waveHeight(x,t):
    return inflowHeightMean + tseries.etaDirect(x,t)*ramp(t)#+ tseries


def waveVelocity_u(x,t):
    return  tseries.uDirect(x,t)[0]*ramp(t)#+ tseries


def waveVelocity_v(x,t):
    return  tseries.uDirect(x,t)[1]*ramp(t)#+ tseries


def waveVelocity_w(x,t):
    return  tseries.uDirect(x,t)[2]*ramp(t)#+ tseries    return waves.u(x[0],x[1],x[2],t,"z")


def wavePhi(x,t):
    return x[2] - waveHeight(x,t)


def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u


def twpflowVelocity_v(x,t):
    return 0.0


def twpflowVelocity_w(x,t):
    waterspeed = waveVelocity_w(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[2]+(1.0-H)*waterspeed


