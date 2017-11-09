import numpy as np
import matplotlib.pyplot as plt
from proteus.mprans import SpatialTools as st
from proteus import Domain
from proteus.mbd import ChRigidBody as crb
from proteus import Context
from proteus import WaveTools as wt
import copy

opts=Context.Options([
    ('water_level', 0.9, 'Height of free surface above bottom'),
    # tank
    ('tank_wavelength_scale', True, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_x', 2., 'Length of tank'),
    ('tank_y', 2., 'Width of tank'),
    ('tank_z', 1.8, 'Height of tank'),
    ('tank_sponge_lambda_scale', True, 'True: sponge length relative to wavelength'),
    ('tank_sponge_xn', 1., 'length of sponge layer x-'),
    ('tank_sponge_xp', 2., 'length of sponge layer x+'),
    ('tank_sponge_yn', 0., 'length of sponge layer y-'),
    ('tank_sponge_yp', 0., 'length of sponge layer y+'),
    ('tank_sponge_gen', 'x-y-y+', 'sponge layers to be gen zones (if waves)'),
    ('tank_sponge_abs', 'x+', 'sponge layers to be abs zones'),
    # bridge
    ('bridge_scale', 10., 'scaling factor of bridge'),
    ('bridge_length', 0.7*10.5625*25*0.0254, 'length of bridge segments'),
    ('bridge_mass', 6350., 'mass of bridge segment'),
    ('bridge_Ixx', [0., 0., 0.], 'moment of inertia of bridge'),
    ('bridge_VCG', 2.1446/2., 'VCG (vertical centre of gravity)'),
    ('bridge_free_x', [1., 1., 1.], 'translational constraints'),
    ('bridge_free_r', [1., 1., 1.], 'rotational constaints'),
    ('bridge_x', 1., 'x position of bridge'),
    ('bridge_y', 1., 'y position of bridge'),
    ('bridge_z', 0.9, 'z position of bridge'),
    # chrono
    ('chrono_dt', 1e-4, 'time step for Chrono'),
    # waves
    ('waves', True, 'Generate waves (True/False)'),
    ('wave_period', 1.4, 'Period of the waves'),
    ('wave_height', 0.08, 'Height of the waves'),
    ('wave_dir', (1., 0., 0.), 'Direction of the waves (from left boundary)'),
    # mesh refinement
    ('he', 0.1, 'Set characteristic element size'),
    # numerical options
    ('genMesh', True, 'True: generate new mesh every time. False: do not generate mesh if file exists'),
    ('use_gmsh', False, 'use_gmsh'),
    ('movingDomain', True, 'True/False'),
    ('T', 50.0, 'Simulation time'),
    ('dt_init', 0.001, 'Initial time step'),
    ('dt_fixed', None, 'Fixed (maximum) time step'),
    ('chrono_dt', 1e-5, 'time step in chrono'),
    ('timeIntegration', 'backwardEuler', 'Time integration scheme (backwardEuler/VBDF)'),
    ('cfl', 0.4, 'Target cfl'),
    ('nsave', 5, 'Number of time steps to save per second'),
    ('useRANS', 0, 'RANS model'),
    ])

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = np.array([0., 0., -9.81])

water_level = opts.water_level

# WAVES

# TO DO: scaling from prototype scale
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = np.array(opts.wave_dir)
    period = opts.wave_period
    BCoeffs = np.zeros(3)
    YCoeffs = np.zeros(3)
    wave = wt.MonochromaticWaves(period=period,
                                 waveHeight=height,
                                 mwl=mwl,
                                 depth=depth,
                                 g=g,
                                 waveDir=direction,
                                 #wavelength=wavelength,
                                 waveType='Linear',
                                 Ycoeff=YCoeffs,
                                 Bcoeff=BCoeffs,
                                 Nf=len(BCoeffs),
                                 fast=False)
    wavelength = wave.wavelength


# DOMAIN
domain = Domain.PiecewiseLinearComplexDomain()

# TANK
tank_dim = [opts.tank_x, opts.tank_y, opts.tank_z]
tank = st.Tank3D(domain, dim=tank_dim)
sponges = {'x-': opts.tank_sponge_xn,
           'x+': opts.tank_sponge_xp,
           'y-': opts.tank_sponge_yn,
           'y+': opts.tank_sponge_yp}
if opts.tank_sponge_lambda_scale and opts.waves:
    for key in sponges:
        sponges[key] *= wave.wavelength
tank.setSponge(x_n=sponges['x-'], x_p=sponges['x+'], y_n=sponges['y-'], y_p=sponges['y+'])
# boundary conditions
tank.BC['z+'].setAtmosphere()
tank.BC['z-'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['y-'].setNoSlip()
tank.BC['y+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()
dragAlpha = 0.5/nu_0
gen = opts.tank_sponge_gen
abso = opts.tank_sponge_abs
if opts.waves is True:
    dragAlpha = 5*(2*np.pi/period)/nu_0
    smoothing = opts.he*3
    tank.setGenerationZones(x_n=('x-' in gen and sponges['x-'] != 0),
                            x_p=('x+' in gen and sponges['x+'] != 0),
                            y_n=('y-' in gen and sponges['y-'] != 0),
                            y_p=('y+' in gen and sponges['y+'] != 0),
                            waves=wave,
                            smoothing=smoothing,
                            dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_n=('x-' in abso and sponges['x-'] != 0),
                        x_p=('x+' in abso and sponges['x+'] != 0),
                        y_n=('y-' in abso and sponges['y-'] != 0),
                        y_p=('y+' in abso and sponges['y+'] != 0),
                        dragAlpha=dragAlpha)


# BRIDGE
b = np.genfromtxt('bridge.csv', delimiter=',', names=True)
# real dim
vertices = []
segments = []
facets = []
for i in range(len(b['x'])):
    vertices += [[0., b['x'][i], b['y'][i]]]
    segments += [[i, i+1]]
segments[-1][1] = 0
facets += [[[i for i in range(len(vertices))]]]
vertices = np.array(vertices)
segments = np.array(segments)
bt = boundaryTags = {'bridge': 1}
vertices *= 25 * 0.0254  # prescale from JIMAjares

def extrude_geom(extrusion, vertices, segments, facets):
    nv = len(vertices)
    ns = len(segments)
    vertices = np.vstack(( vertices, vertices+extrusion ))
    segments_side = [[i, i+nv] for i in range(nv)]
    segments_extruded = segments+[ns, ns]
    segments = np.vstack((segments, segments_side, segments_extruded))
    facets_side = [[[i, i+1, i+nv+1, i+nv]] for i in range(nv)]
    facets_side[-1][0][1] = 0
    facets_side[-1][0][2] = nv
    facets_extruded = copy.deepcopy(facets)
    for facet in facets_extruded:
        for subfacet in facet:
            for i in range(len(subfacet)):
                subfacet[i] += nv
    facets = facets+facets_side+facets_extruded
    volumes = [[[i for i in range( len(facets) )]]]
    return vertices, segments, facets, volumes
bridge_length = opts.bridge_length
extrusion = [bridge_length, 0., 0.]
vertices, segments, facets, volumes = extrude_geom(extrusion, vertices, segments, facets)
vertexFlags = np.array([bt['bridge'] for i in range(len(vertices))])
segmentFlags = np.array([bt['bridge'] for i in range(len(segments))])
facetFlags = np.array([bt['bridge'] for i in range(len(facets))])
offset_x = 0.5 * (max(vertices[:, 0]) + min(vertices[:, 0]))
offset_y = 0.5 * (max(vertices[:, 1]) + min(vertices[:, 1]))
offset_z = 0.5 * (max(vertices[:, 2]) + min(vertices[:, 2]))
regions = np.array([[offset_x, offset_y, offset_z]])
regionFlags = np.array([1])
holes = regions
barycenter = np.array([0., 0., 0.])
bridge = st.CustomShape(domain,
                        barycenter=barycenter,
                        vertices=vertices,
                        vertexFlags=vertexFlags,
                        segments=segments,
                        segmentFlags=segmentFlags,
                        facets=facets,
                        facetFlags=facetFlags,
                        regions=regions,
                        regionFlags=regionFlags,
                        boundaryTags=boundaryTags,
                        volumes=volumes,
                        holes=holes)
bridge.holes_ind = np.array([0])
tank.setChildShape(bridge)  # volume information gmsh


# translate to be centered at 0.
bridge.translate([-offset_x, -offset_y, -offset_z])
# set barycenter according to VCG
range_z = (max(bridge.vertices[:, 2]-min(bridge.vertices[:, 2])))
VCG = opts.bridge_VCG
bridge.setBarycenter([0., 0., VCG-range_z/2.])
# scale it
scale = opts.bridge_scale
bridge.vertices *= 1./scale
#translate to desired height
brige_height = water_level
bridge.translate([opts.bridge_x, opts.bridge_y, opts.bridge_z])

# boundary conditions
bridge.BC['bridge'].setNoSlip()

# ASSEMBLE GEOMETRIES
# necessary once all egometries have been created
he = opts.he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setElementSize(he)
domain.MeshOptions.setTriangleOptions()
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)
domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
domain.use_gmsh = opts.use_gmsh
st.assembleDomain(domain)

# CHRONO
system = crb.ProtChSystem(g)
system.setTimeStep(opts.chrono_dt)

body = crb.ProtChBody(system, shape=bridge)
# by accessing the ChBody attribute of ProtChBody, standard chrono functions can be used
body.ChBody.SetMass(opts.bridge_mass/opts.bridge_scale)
Ixx = opts.bridge_Ixx
range_x = (max(bridge.vertices[:, 0]-min(bridge.vertices[:, 0])))
range_y = (max(bridge.vertices[:, 1]-min(bridge.vertices[:, 1])))
range_z = (max(bridge.vertices[:, 2]-min(bridge.vertices[:, 2])))
if Ixx[0] == 0:
    Ixx[0] = (range_y**2+range_z**2)/12.
if Ixx[1] == 0:
    Ixx[1] = (range_x**2+range_z**2)/12.
if Ixx[2] == 0:
    Ixx[2] = (range_y**2+range_x**2)/12.
from proteus.mbd import pyChronoCore as pych
Ixx_vec = pych.ChVector(Ixx[0], Ixx[1], Ixx[2])
body.ChBody.SetInertiaXX(Ixx_vec)
body.setConstraints(free_x=np.array(opts.bridge_free_x),
                    free_r=np.array(opts.bridge_free_r))
body.setRecordValues(all_values=True)



from MeshRefinement import geometry_to_gmsh
mesh = geometry_to_gmsh(domain)
mesh.Options.Mesh.CharacteristicLengthMax = he
mesh.writeGeo(mesh_fileprefix+'.geo')










##########################################
# Numerical Options and other parameters #
##########################################

from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10./nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = opts.useRANS # 0 -- None
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
nd = domain.nd
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
sc = 0.5
sc_beta = 1.5
epsFact_consrv_diffusion = 10.0
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
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

ns_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-6,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-6,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-6,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-6,opts.mesh_tol*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - water_level
    phi = x[nd-1] - water_level
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))

#isosurface = None
from proteus.Isosurface import Isosurface
#isosurface = Isosurface(isosurfaces=(('phi', (0.,)),), domain=domain, format='h5', sampleRate=opts.sampleRate)
isosurface = None


def load_simulation_globals():
    """Put some variables we need in engine namespace.

    These can then be retrieved by clients for inspection, visualization, etc.
    """
    nodes = isosurface.nodes_array
    triangles = isosurface.elements_array
    x = nodes[:,0]
    y = nodes[:,1]
    z = nodes[:,2]
    vertices = nodes
    nn = len(x)
