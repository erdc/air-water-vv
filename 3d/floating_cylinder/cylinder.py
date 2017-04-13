from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.mbd import ChRigidBody as crb
import numpy as np

from proteus import Comm
comm=Comm.init()

opts=Context.Options([
    ("water_level", 1., "Height of free surface above bottom"),
    # tank
    ("tank_x", 2., "Length of tank"),
    ("tank_y", 2., "Width of tank"),
    ("tank_z", 2., "Height of tank"),
    # cylinder
    ("cylinder", True, "cylinder"),
    ("cylinder_radius", 0.515/2., "radius of cylinder"),
    ("cylinder_height", 0.344, "radius of cylinder"),
    ("cylinder_coords", (1., 1., 1.), "coordinates of cylinder"),
    ("cylinder_mass", 35.85, "mass of cylinder"),
    ("free_x", (1., 1., 1.), "Translational DOFs"),
    ("free_r", (1., 1., 1.), "Rotational DOFs"),
    ("VCG", 0.0758, "VCG"),
    ("moorings", True, "moorings"),
    # mesh refinement
    ("he", 0.2, "Set characteristic element size"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "use_gmsh"),
    ("refinement_freesurface", 0.1, "ref"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("chrono_dt", 1e-5, "time step in chrono"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("sc", 0.25, "shockCapturing factor"),
    ("parallel", True ,"Run in parallel")])


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., 0., -9.81]

# ----- CONTEXT ------ #

# general options
water_level = opts.water_level

# tank options
tank_dim = [opts.tank_x, opts.tank_y, opts.tank_z]
cylinder_coords = opts.cylinder_coords
cylinder_radius = opts.cylinder_radius
cylinder_height = opts.cylinder_height
free_x = np.array(opts.free_x)
free_r = np.array(opts.free_r)
chrono_dt = opts.chrono_dt


# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()


# ----- SHAPES ----- #

tank = st.Cuboid(domain, dim=tank_dim, coords=np.array(tank_dim)/2.)
#tank = st.Tank3D(domain, dim=(5., 5., 1.8))
barycenter = cylinder_coords-np.array([0.,0.,cylinder_height/2.+opts.VCG])
cylinder = st.Cylinder(domain, radius=cylinder_radius, height=cylinder_height, coords=cylinder_coords, barycenter=barycenter, nPoints=20)
#cylinder.translate(cylinder_coords)
cylinder.setHoles([cylinder_coords])
cylinder.holes_ind = np.array([0])
cylinder.setBarycenter(barycenter)

tank.setChildShape(cylinder, 0)



# ----- CHRONO ----- #

g = np.array([0., 0., -9.81])
system = crb.System(g)
system.setTimeStep(chrono_dt)
mesh = crb.Mesh(system)

if opts.cylinder is True:
    body = crb.RigidBody(system, cylinder, free_x=free_x, free_r=free_r)
    body.setMass(opts.cylinder_mass)
    body.setRecordValues(all_values=True)

if opts.moorings is True:
    mooring_X = 6.660
    cylinder_radius = 0.5
    cylinder_coords = [1.,1.,1.]
    water_level=1.

    anchor = np.array([6.660, 0., 0.9])
    fairlead1_offset = np.array([cylinder_radius*np.cos(np.radians(120)), cylinder_radius*np.sin(np.radians(120)), 0.])
    fairlead2_offset = np.array([cylinder_radius,0.,0.])
    fairlead3_offset = np.array([cylinder_radius*np.cos(np.radians(-120)), cylinder_radius*np.sin(np.radians(-120)), 0.])
    fairlead_center = np.array([cylinder_coords[0], cylinder_coords[1], water_level])

    fairlead1 = fairlead_center + fairlead1_offset
    fairlead2 = fairlead_center + fairlead2_offset
    fairlead3 = fairlead_center + fairlead3_offset

    anchor1 = fairlead1-[0.,0.,water_level] + np.array([mooring_X*np.cos(np.radians(120)), mooring_X*np.sin(np.radians(120)), 0.])
    anchor2 = fairlead2-[0.,0.,water_level] + np.array([mooring_X, 0., 0.])
    anchor3 = fairlead3-[0.,0.,water_level] + np.array([mooring_X*np.cos(np.radians(-120)), mooring_X*np.sin(np.radians(-120)), 0.])
        
    from catenary import MooringLine

    L = 6.950
    w = 1.243
    EA = 1e20
    d = 4.786e-3
    nb_elems = 20
    dens = 0.1446/(np.pi*d**2/4)
    E = 1e9
    l1 = MooringLine(L=L, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1, tol=1e-8)
    l2 = MooringLine(L=L, w=w, EA=EA, anchor=anchor2, fairlead=fairlead2, tol=1e-8)
    l3 = MooringLine(L=L, w=w, EA=EA, anchor=anchor3, fairlead=fairlead3, tol=1e-8)

    l1.setVariables()
    l2.setVariables()
    l3.setVariables()

    L = np.array([L])
    d = np.array([d])
    nb_elems = np.array([nb_elems])
    E = np.array([E])
    dens = np.array([dens])
    m1 = crb.Moorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E)
    m2 = crb.Moorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E)
    m3 = crb.Moorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=dens, E=E)

    m1.setName('mooring1')
    m2.setName('mooring2')
    m3.setName('mooring3')

    m1.setNodesPositionFunction(l1.s)
    m2.setNodesPositionFunction(l2.s)
    m3.setNodesPositionFunction(l3.s)
    m1.setNodesPosition()
    m2.setNodesPosition()
    m3.setNodesPosition()

    m1.buildNodes()
    m2.buildNodes()
    m3.buildNodes()

    m1.fixFrontNode(True)
    m2.fixFrontNode(True)
    m3.fixFrontNode(True)

    #m1.fixBackNode(True)
    #m2.fixBackNode(True)
    #m3.fixBackNode(True)

    if opts.cylinder is True:
        m1.attachBackNodeToBody(body)
        m2.attachBackNodeToBody(body)
        m3.attachBackNodeToBody(body)

    # FLOOR
    pos1 = m1.getNodesPosition()
    pos2 = m2.getNodesPosition()
    pos3 = m3.getNodesPosition()

    material = crb.ChMaterialSurfaceDEM()
    material.SetYoungModulus(2e4)
    material.SetFriction(0.3)
    material.SetRestitution(0.2)
    material.SetAdhesion(0)

    m1.setContactMaterial(material)
    m2.setContactMaterial(material)
    m3.setContactMaterial(material)

    box_pos = np.array([0.,0.,-0.11])
    box_dim = np.array([20.,1.,0.2])
    vec = crb.ChVector(0., 0., -0.11)
    box = crb.ChBodyEasyBox(system, box_dim[0], box_dim[1], box_dim[2], 1000, True)
    box.SetPos(vec)
    box.SetMaterialSurface(material)
    box.SetBodyFixed(True)


# ----- BOUNDARY CONDITIONS ----- #

# to change individual BC functions, example:
# tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 3.*t
tank.BC['z+'].setAtmosphere()
tank.BC['z-'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['y-'].setNoSlip()
tank.BC['y+'].setNoSlip()
for key, bc in tank.BC.items():
    # fix the nodes on the wall of tank
    # in case of moving domain
    bc.setFixedNodes()

for key, bc in cylinder.BC.items():
    bc.setNoSlip()
# moving mesh BC were created automatically when
# making a chrono rigid body for the cylinder

domain.MeshOptions.he = opts.he
domain.MeshOptions.setTriangleOptions()
st.assembleDomain(domain)  # must be called after defining shapes

he = opts.he
domain.writePoly('mesh')








tank_sponge = [0,0]



from py2gmsh2 import Fields, Entity
from MeshRefinement import geometry_to_gmsh
mesh = geometry_to_gmsh(domain)
grading = 1.03
he = opts.he
he_max = 10.
he_max_water = 10.
ecH = 3.
if opts.refinement_freesurface > 0:
    box = opts.refinement_freesurface
else:
    box = ecH*he
field_list = []

# refinement free surface
box1 = Fields.Box(mesh=mesh)
box1.VIn = he
box1.VOut = he_max
box1.XMin = -tank_sponge[0]
box1.XMax = tank_dim[0]+tank_sponge[1]
box1.YMin = water_level-box
box1.YMax = water_level+box
field_list += [box1]

#p0 = Entity.Point([0., 0., water_level+box], mesh=mesh)
#p1 = Entity.Point([0., tank_dim[1], water_level+box], mesh=mesh)
#p2 = Entity.Point([tank_dim[0], tank_dim[1], water_level+box], mesh=mesh)
#p3 = Entity.Point([tank_dim[0], 0., water_level+box], mesh=mesh)
#l1 = Entity.Line([p0, p1], mesh=mesh)
#l2 = Entity.Line([p1, p2], mesh=mesh)
#l3 = Entity.Line([p2, p3], mesh=mesh)
#l4 = Entity.Line([p3, p0], mesh=mesh)
#ll1 = Entity.LineLoop([l1, l2, l3, l4], mesh=mesh)
#s1 = Entity.PlaneSurface([ll1], mesh=mesh)
#
#p10 = Entity.Point([0., 0., water_level-box], mesh=mesh)
#p11 = Entity.Point([0., tank_dim[1], water_level-box], mesh=mesh)
#p12 = Entity.Point([tank_dim[0], tank_dim[1], water_level-box], mesh=mesh)
#p13 = Entity.Point([tank_dim[0], 0., water_level-box], mesh=mesh)
#l11 = Entity.Line([p10, p11], mesh=mesh)
#l12 = Entity.Line([p11, p12], mesh=mesh)
#l13 = Entity.Line([p12, p13], mesh=mesh)
#l14 = Entity.Line([p13, p10], mesh=mesh)
#ll2 = Entity.LineLoop([l11, l12, l13, l14], mesh=mesh)
#s2 = Entity.PlaneSurface([ll2], mesh=mesh)
#
#bl2 = Fields.BoundaryLayer(mesh=mesh)
#bl2.hwall_n = he
#bl2.ratio = grading
#bl2.EdgesList = [l1, l2, l3, l4]
#field_list += [bl2]
#
#bl2 = Fields.BoundaryLayer(mesh=mesh)
#bl2.hwall_n = he
#bl2.ratio = grading
#bl2.EdgesList = [l11, l12, l13, l14]
#field_list += [bl2]

# max element size in water phase
box2 = Fields.Box(mesh=mesh)
box2.VIn = he_max_water
box2.VOut = he_max
box2.XMin = -tank_sponge[0]
box2.XMax = tank_dim[0]+tank_sponge[1]
box2.YMin = 0
box2.YMax = water_level
field_list += [box2]

if opts.cylinder:
    # boundary layer on caisson
    bl1 = Fields.BoundaryLayer()
    bl1.hwall_n = he
    bl1.ratio = grading
    bl1.EdgesList = mesh.getLinesFromIndex([i+1 for i in range(len(cylinder.segments))])
    mesh.addField(bl1)
    field_list += [bl1]

#if opts.refinement_caisson:
#    # create circle (non-physical) around caisson
#    refinement_caisson = opts.refinement_caisson
#    p0 = .Entity.Point([caisson_coords[0], caisson_coords[1], 0.], mesh=mesh)
#    p1 = .Entity.Point([caisson_coords[0]-refinement_caisson, caisson_coords[1], 0.], mesh=mesh)
#    p2 = .Entity.Point([caisson_coords[0]+refinement_caisson, caisson_coords[1], 0.], mesh=mesh)
#    p3 = .Entity.Point([caisson_coords[0]-refinement_caisson+0.00001, caisson_coords[1], 0.], mesh=mesh)
#    c1 = .Entity.Circle(p1, p0, p2, nb=100, mesh=mesh)
#    c2 = .Entity.Circle(p2, p0, p3, nb=101, mesh=mesh)
#
#    # refined circle around caisson
#    b1 = .Fields.Ball(mesh=mesh)
#    b1.VIn = he
#    b1.VOut = he_max
#    b1.Radius = refinement_caisson
#    b1.XCenter = caisson_coords[0]
#    b1.YCenter = caisson_coords[1]
#    b1.ZCenter = 0.
#    field_list += [b1]
#
#    # boundary layer on circle around caisson
#    bl3 = .Fields.BoundaryLayer(mesh=mesh)
#    bl3.hwall_n = he
#    bl3.ratio = grading
#    bl3.EdgesList = [c1, c2]
#    field_list += [bl3]

# background field
fmin = Fields.Min(mesh=mesh)
fmin.FieldsList = field_list
mesh.setBackgroundField(fmin)

# max element size
mesh.Options.Mesh.CharacteristicLengthMax = he_max

geofile = 'mesh'
mesh.writeGeo(geofile+'.geo')

domain.use_gmsh = opts.use_gmsh
domain.geofile = geofile+'.geo'
domain.MeshOptions.use_gmsh = opts.use_gmsh



















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

ns_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-12,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-12,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-12,opts.mesh_tol*he**2)

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
