import matplotlib.pyplot as plt
import numpy as np
from proteus.iproteus import *
import proteus.default_p as physics
import proteus.default_n as numerics
from proteus.TransportCoefficients import PoissonEquationCoefficients
from proteus  import Profiling
from proteus.Profiling import logEvent
from proteus import WaveTools as wt
Profiling.logLevel=7
Profiling.verbose=True


opts=Context.Options([
    # predefined test cases
    ("water_level", 0.5, "Height of free surface above bottom"),
    # tank
    ('tank_wavelength_scale', False, 'if True: tank_x=value*wavelength, tank_y=value*wavelength'),
    ('tank_x', 1., 'Length of tank'),
    ('tank_y', 1., 'Width of tank'),
    ('tank_sponge_lambda_scale', True, 'True: sponge length relative to wavelength'),
    ('tank_sponge_xn', 0., 'length of sponge layer x-'),
    ('tank_sponge_xp', 0., 'length of sponge layer x+'),
    ("tank_sponge", (0., 0.), "Length of absorption zones (front/back, left/right)"),
    ('tank_sponge_gen', 'x-', 'sponge layers to be gen zones (if waves)'),
    ('tank_sponge_abs', 'x+', 'sponge layers to be abs zones'),
    # waves
    ("waves", False, "Generate waves (True/False)"),
    ("wave_period", 0.4, "Period of the waves"),
    ("wave_height", 0.01, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Ditankion of the waves (from left boundary)"),
    # numerical options
    ("he", 0.02, "Set characteristic element size"),
    ("addedMass", False, "Set characteristic element size"),
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ])

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = np.array([0., 0., -9.81])

waterLevel = water_level = opts.water_level

wavelength=1.
# general options

if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = np.array(opts.wave_dir)
    period = opts.wave_period
    wave = wt.MonochromaticWaves(period=period,
                                 waveHeight=height,
                                 mwl=mwl,
                                 depth=depth,
                                 g=g,
                                 waveDir=direction,
                                 wavelength=wavelength,
                                 waveType='Fenton',
                                 Nf=8,
                                 fast=False)
    wavelength = wave.wavelength

def get_center_area(e_nodes):
    detJ = (e_nodes[1][0] - e_nodes[0][0]) * (e_nodes[2][1] -
                                              e_nodes[0][1]) - (e_nodes[1][1] - e_nodes[0][1]) * (e_nodes[2][0] -
                                                                                                  e_nodes[0][0])
    # since the orientation is clockwise
    center = ( e_nodes[0]+e_nodes[1]+e_nodes[2] )/3.
    area = 0.5 * np.abs(detJ)
    return area, center

import copy

he_max=1.
he_min=0.05
r = 0.1
nSmooth = 10



nd = 2

# use structured mesh
# domain =  Domain.tankangularDomain(L=[1.0,1.0], x=[0.,0.])
# nn = 50
from proteus.mprans import SpatialTools as st
domain = Domain.PlanarStraightLineGraphDomain()


# TANK
tank_dim = np.array([opts.tank_x, opts.tank_y])
if opts.tank_wavelength_scale and opts.waves:
    tank_dim[0:2] *= wavelength
tank = st.Tank2D(domain, tank_dim)
sponges = {'x-': opts.tank_sponge_xn,
           'x+': opts.tank_sponge_xp}
if opts.tank_sponge_lambda_scale and opts.waves:
    for key in sponges:
        sponges[key] *= wave.wavelength
tank.setSponge(x_n=sponges['x-'], x_p=sponges['x+'])
# to change individual BC functions, example:
# tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 3.*t
for key, bc in tank.BC.items():
    # fix the nodes on the wall of tank
    # in case of moving domain
    bc.setFixedNodes()
tank.BC['x-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['y-'].setFreeSlip()
tank.BC['y+'].setAtmosphere()
tank.BC['sponge'].setNonMaterial()
dragAlpha = 0.5/nu_0
gen = opts.tank_sponge_gen
abso = opts.tank_sponge_abs
if opts.waves is True:
    dragAlpha = 5*(2*np.pi/period)/nu_0
    smoothing = opts.he*3
    tank.setGenerationZones(x_n=('x-' in gen and sponges['x-'] != 0),
                            x_p=('x+' in gen and sponges['x+'] != 0),
                            waves=wave,
                            smoothing=smoothing,
                            dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_n=('x-' in abso and sponges['x-'] != 0),
                        x_p=('x+' in abso and sponges['x+'] != 0),
                        dragAlpha=dragAlpha)

boundaryNormals = {tank.boundaryTags['x-']: np.array([-1.,0.]),
                   tank.boundaryTags['x+']: np.array([1.,0.]),
                   tank.boundaryTags['y-']: np.array([0.,-1.]),
                   tank.boundaryTags['y+']: np.array([0.,1.])}
boundaryNormals = None

# center1 = [1.5-2*r, 1.5-0.1]
# center1 = [tank.dim[0]/2., tank.dim[1]/2.+0.2]
# center2 = [tank.dim[0]/2., tank.dim[1]/2.-0.2]
# center3 = [tank.dim[0]/2.+0.2, tank.dim[1]/2.]
# center4 = [tank.dim[0]/2.-0.2, tank.dim[1]/2.]

center1 = [tank_dim[0]/2., tank_dim[1]/2.]

st.assembleDomain(domain)
he = opts.he
domain.MeshOptions.he = he
domain.MeshOptions.setTriangleOptions()
domain.MeshOptions.setOutputFiles(name="mesh")
domain.use_gmsh = False
genMesh = True

use_gmsh = opts.use_gmsh
if use_gmsh:
    from MeshRefinement import geometry_to_gmsh
    import py2gmsh
    domain.use_gmsh = True
    mesh = geometry_to_gmsh(domain)
    field_list = []
    me01 = py2gmsh.Fields.MathEval(mesh=mesh)
    me01.F = "{he}+{he}*10*abs(sqrt((y-{water_level})^2))".format(he=he/5., water_level=water_level)
    field_list += [me01]
    # me02 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me02.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center2[0], center_y=center2[1], radius=r)
    # field_list += [me02]
    # me03 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me03.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center3[0], center_y=center3[1], radius=r)
    # field_list += [me03]
    # me04 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me04.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center4[0], center_y=center4[1], radius=r)
    # field_list += [me04]
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)
    mesh.Options.Mesh.CharacteristicLengthMax = he*(he_max/he_min)
    mesh.Options.Mesh.CharacteristicLengthMax = he*5
    mesh.writeGeo("mesh.geo")
    domain.geofile = "mesh"


def my_func(x):
    return min(he_max, max(np.abs(np.sqrt((x[0]-0.5)**2+(x[1]-0.5)**2)-r), he_min))
scale = 0.5
def my_func(x, t):
    t_start = 1000.
    if t>t_start:
        dist1 = np.sqrt((x[0]-(center1[0]))**2+(x[1]-(center1[1]-(t-t_start)/100))**2)-r
        # dist2 = np.sqrt((x[0]-(center2[0]))**2+(x[1]-(center2[1]+(t-t_start)/100))**2)-r
        # dist3 = np.sqrt((x[0]-(center3[0])+(t-t_start)/100)**2+(x[1]-(center3[1]))**2)-r
        # dist4 = np.sqrt((x[0]-(center4[0])-(t-t_start)/100)**2+(x[1]-(center4[1]))**2)-r
        circle1 = abs(dist1)/scale
        # circle2 = min(he_max, max(abs(dist2)/scale, he_min))
        # circle3 = min(he_max, max(abs(dist3)/scale, he_min))
        # circle4 = min(he_max, max(abs(dist4)/scale, he_min))
    else:
        dist1 = np.sqrt((x[0]-center1[0])**2+(x[1]-center1[1])**2)-r
        # dist2 = np.sqrt((x[0]-center2[0])**2+(x[1]-center2[1])**2)-r
        # dist3 = np.sqrt((x[0]-center3[0])**2+(x[1]-center3[1])**2)-r
        # dist4 = np.sqrt((x[0]-center4[0])**2+(x[1]-center4[1])**2)-r
        circle1 = abs(dist1)/scale
        # circle2 = min(he_max, max(abs(dist2)/scale, he_min))
        # circle3 = min(he_max, max(abs(dist3)/scale, he_min))
        # circle4 = min(he_max, max(abs(dist4)/scale, he_min))
    #border = min(he_max, min(max(abs(x[0]-0.)/scale, he_min),
    #                         max(abs(x[0]-domain.L[0])/scale, he_min),
    #                         max(abs(x[1]-0.)/scale, he_min),
    #                         max(abs(x[1]-domain.L[1])/scale, he_min)))
    # return min(circle1, circle2, circle3, circle4) #min(border, circle)
    return he_max#circle1#min(circle1, circle2) #min(border, circle)


movingDomain = True




from petsc4py import PETSc
OptDB = PETSc.Options()
#OptDB.setValue("ksp_type", "cg")
#OptDB.setValue("pc_type", "jacobi")
#OptDB.setValue("pc_factor_mat_solver_package", "superlu_dist")
OptDB.setValue("ksp_constant_null_space", 1)
#OptDB.setValue("ksp_converged_true_residual", 1)
OptDB.setValue("info", 1)
#OptDB.setValue("pc_factor_shift_type","NONZERO")
#OptDB.setValue("pc_factor_shift_amount", 1e-10)
print(OptDB.getAll())



from proteus.default_n import *
movingDomain = True
T = opts.T


spaceOrder = 1
useHex = False
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


# model = ns.modelList[0].levelModelList[-1]

# x0 = model.mesh_array0[:,0]
# y0 = model.mesh_array0[:,1]
# x1 = model.mesh.nodeArray[:,0]
# y1 = model.mesh.nodeArray[:,1]
# z = model.u[0].dof
# qx0 = model.q['x'][:,:,0].flatten()
# qy0 = model.q['x'][:,:,1].flatten()
# qzx = model.q[('grad(u)',0)][:,:,0].flatten()
# qzy = model.q[('grad(u)',0)][:,:,1].flatten()
# Nlevels = 1000

# import matplotlib.cm as cm

# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax1.set_aspect('equal')
# ax2.tricontourf(x1,
#                 y1,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax2.set_aspect('equal')
# print('min/max',min(model.u[0].dof), max(model.u[0].dof))
# c0 = CS.collections[0]
# f.colorbar(CS, ax=ax1)
# f.savefig('sol.png', bbox_inches='tight', dpi=100)

# # plt.figure()
# # plt.plot(model.q['x'][:,:,0], model.q['x'][:,:,1], 'b.')

# import matplotlib.tri as tri
# # Create the Triangulation; no triangles so Delaunay triangulation created.
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
# ax1.triplot(x0, y0, model.mesh.elementNodesArray, lw=0.8)
# ax1.set_aspect('equal')
# ax2.triplot(x1, y1, model.mesh.elementNodesArray, lw=0.8)
# ax2.set_aspect('equal')
# f.savefig('mesh.png', bbox_inches='tight', dpi=100)

# from mpl_toolkits.mplot3d import Axes3D
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,subplot_kw=dict(projection='3d'))
# fig = plt.figure()
# #ax = fig.add_subplot(111, projection='3d')
# ax1.plot_trisurf(x0, y0, z)
# ax2.plot_trisurf(x1, y1, z)
# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax1.set_axis_off()
# ax1.set_aspect('equal')
# f.savefig('a.png', bbox_inches='tight', dpi=100)

# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 model.grads[:,0],
#                 cmap=cm.jet)
# f.savefig('gradx.png', bbox_inches='tight', dpi=100)

# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                     y0,
#                     model.mesh.elementNodesArray,
#                     model.grads[:,1],
#                     Nlevels,
#                     cmap=cm.jet)
# f.savefig('grady.png', bbox_inches='tight', dpi=100)


# def my_func(x, y):
#     return np.minimum(he_max, np.maximum(np.abs(np.sqrt((x-0.5)**2+(y-0.5)**2)-0.25)/0.25, he_min))
# x = np.linspace(0., 1., 100)
# y = np.linspace(0., 1., 100)
# X, Y = np.meshgrid(x, y)
# Z = my_func(X, Y)
# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10))
# CS = ax1.contourf(X,
#                     Y,
#                     Z,
#                     Nlevels,
#                     cmap=cm.jet)
# f.colorbar(CS, ax=ax1)
# f.savefig('func.png', bbox_inches='tight', dpi=100)

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
    dt_out= (T-opts.dt_init)/nDTout
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
am_nl_atol_res = 1e-6 #max(1.0e-6,opts.mesh_tol*he**2)



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

