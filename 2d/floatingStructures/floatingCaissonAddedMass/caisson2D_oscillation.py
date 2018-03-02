from proteus import Domain, Context, Comm
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from proteus.Profiling import logEvent
#import ChRigidBody as crb
from proteus.mbd import ChRigidBody as crb
from math import *
import numpy as np

opts=Context.Options([
    ("cf",1.0,"Use corrected form of added mass"),
    ("nc_form",True,"Use non-conservative NSE form"),
    ("use_chrono", True, "use chrono (True) or custom solver"),
    # predefined test cases
    ("water_level", 1.5, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (3., 3.,), "Dimensions of the tank"),
    ("tank_sponge", (0., 0.), "Length of absorption zones (front/back, left/right)"),
    # caisson
    ("addedMass", True, "added mass"),
    ("caisson", True, "caisson"),
    ("caisson_dim", (0.5, 0.5), "Dimensions of the caisson"),
    ("caisson_coords", (1.5, 1.5+0.5), "Dimensions of the caisson"),
    ("free_x", (0.0, 1.0, 0.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 0.0), "Rotational DOFs"),
    ("VCG", 0.05, "vertical position of the barycenter of the caisson"),
    ("caisson_mass", 125., "Mass of the caisson"),
    ("caisson_inertia", 4.05, "Inertia of the caisson"),
    ("rotation_angle", 0., "Initial rotation angle (in degrees)"),
    ("chrono_dt", 0.00001, "time step of chrono"),
    # mesh refinement
    ("refinement", True, "Gradual refinement"),
    ("he", 0.03, "Set characteristic element size"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("genMesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave", 5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("sc", 0.25, "shockCapturing factor"),
    ("weak_factor", 10., "weak bc penalty factor"),
    ("strong_dir", False, "strong dirichlet (True/False)"),
    ])



# ----- CONTEXT ------ #

# general options
waterLevel = water_level = opts.water_level
rotation_angle = np.radians(opts.rotation_angle)


# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options
if opts.caisson is True:
    free_x = opts.free_x
    free_r = opts.free_r
    rotation = np.radians(opts.rotation_angle)
    inertia = opts.caisson_inertia

    caisson_dim = opts.caisson_dim
    caisson_name='caisson_chrono'
    if opts.addedMass:
        caisson_name+='_am'
    if opts.nc_form:
        caisson_name+='_nc'
    caisson = st.Rectangle(domain, dim=opts.caisson_dim, coords=[0.,0.], barycenter=np.zeros(3))
    caisson.name=caisson_name
    ang = rotation_angle
    caisson.setHoles([[0., 0.]])
    caisson.holes_ind = np.array([0])
    print(caisson.regions+np.ones(2))
    print(caisson.regionFlags)

    trans = np.array([opts.caisson_coords[0], opts.caisson_coords[1]])
    print(trans+np.ones(2))
    caisson.translate(trans)
    # system = crb.System(np.array([0., -9.81, 0.]))
    # rotation = np.array([1, 0., 0., 0.])
    rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)
    for bc in caisson.BC_list:
        bc.setNoSlip()
    if opts.use_chrono:
        system = crb.ProtChSystem(np.array([0., -9.81, 0.]))
        system.setTimeStep(opts.chrono_dt)
        system.step_start = 10
        body = crb.ProtChBody(shape=caisson,
                              system=system)
        chbod = body.ChBody
        from proteus.mbd import pyChronoCore as pych
        x, y, z = caisson.barycenter
        pos = pych.ChVector(x, y, z)
        e0, e1, e2, e3 = rotation_init
        rot = pych.ChQuaternion(e0, e1, e2, e3)
        inertia = pych.ChVector(1., 1., inertia)
        chbod.SetPos(pos)
        chbod.SetRot(rot)
        chbod.SetMass(opts.caisson_mass)
        chbod.SetInertiaXX(inertia)
        body.setConstraints(free_x=np.array(opts.free_x), free_r=np.array(opts.free_r))
        # body.setInitialRot(rotation_init)
        # body.rotation_init=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
        body.setRecordValues(all_values=True)
    else:
        from proteus import AuxiliaryVariables
        class BodyAddedMass(AuxiliaryVariables.AV_base):
            def __init__(self):
                AuxiliaryVariables.AV_base.__init__(self)

            def attachModel(self, model, ar):
                self.model=model
                return self

            def calculate(self):
                pass
        from proteus.mprans import BodyDynamics as bd
        body = bd.RigidBody(shape=caisson)
        system = body  # so it is added in twp_navier_stokes_p.py
        bodyAddedMass = BodyAddedMass()
        body.bodyAddedMass = body.ProtChAddedMass = bodyAddedMass
        body.setMass(opts.caisson_mass)
        body.setConstraints(free_x=free_x,
                            free_r=free_r)
        body.It = np.array([[1., 0., 0.],
                            [0., 1., 0.],
                            [0., 0., inertia]])
        filename='record_caisson'
        if opts.addedMass:
            filename+='_am'
        if opts.nc_form:
            filename+='_nc'
        if opts.cf:
            filename+='_cf'
        body.setRecordValues(filename=filename, all_values=True)
        body.coords_system = caisson.coords_system  # hack
        body.last_Aij = np.zeros((6,6),'d')
        body.last_a = np.zeros((6,),'d')
        body.a = np.zeros((6,),'d')
        body.last_Omega = np.zeros((3,3),'d')
        body.last_velocity = np.zeros((3,),'d')
        body.last_Q = np.eye(3)
        body.last_h = np.zeros((3,),'d')
        body.last_mom = np.zeros((6,),'d')
        body.last_u = np.zeros((18,),'d')
        body.last_t = 0.0
        body.t = 0.0
        body.h = np.zeros((3,),'d')
        body.velocity = np.zeros((3,),'d')
        body.mom = np.zeros((6,),'d')
        body.Omega = np.zeros((3,3),'d')
        body.Q = np.eye(3,3)
        body.u = np.zeros((18,),'d')
        body.u[9:18]= body.Q.flatten()
        body.last_u[:] = body.u
        body.last_Q[:] = body.Q
        body.FT = np.zeros((6,),'d')
        body.last_FT = np.zeros((6,),'d')
        body.free_dof = np.zeros((6,),'d')
        body.free_dof[:3] = free_x
        body.free_dof[3:] = free_r
        # cek -- I  did this with a different hack, see dummy AddedMassBody above
        # body.Aij = np.zeros((6,6))

        # # we create a ProtChAddedMass auxiliary variable to add to added_mass_n.py
        # # the argument here should be a ProtChSystem, so we give it an empty system (hack)
        # # this is just to get access to the AddedMass model from the body functions
        # system = crb.ProtChSystem(g)
        # body.ProtChAddedMass = crb.ProtChAddedMass(system)

        def step(dt):
            from math import ceil
            logEvent("Barycenter "+str(body.barycenter))
            n = max(1.0,ceil(dt/opts.chrono_dt))
            DT=dt/n
            def F(u,theta):
                """The residual and Jacobian for discrete 6DOF motion with added mass"""
                v = u[:3]
                omega = u[3:6]
                h = u[6:9]
                Q = u[9:18].reshape(3,3)
                Omega = np.array([[      0.0, -omega[2],  omega[1]],
                                  [ omega[2],       0.0, -omega[0]],
                                  [-omega[1],  omega[0],      0.0]])
                I = np.matmul(np.matmul(Q, body.It), Q.transpose())
                body.Aij = np.zeros((6,6),'d')
                if opts.addedMass:
                    for i in range(1,5):#number of rigid body facets
                        body.Aij += body.bodyAddedMass.model.levelModelList[-1].Aij[i]
                avg_Aij=False
                if avg_Aij:
                    M = body.Aij*theta + body.last_Aij*(1-theta)
                else:
                    M = body.Aij.copy()
                for i in range(6):
                    for j in range(6):
                        M[i,j]*=body.free_dof[j]#only allow j accelerations to contribute to i force balance if j is free
                        M[j,i]*=body.free_dof[j]#only allow j added mass forces if j is free
                        body.Aij[i,j]*=body.free_dof[j]#only allow j accelerations to contribute to i force balance if j is free
                        body.Aij[j,i]*=body.free_dof[j]#only allow j added mass forces if j is free
                body.FT[:3] = body.F
                body.FT[3:] = body.M
                #cek debug
                #body.FT[:] = 0.0
                #body.FT[1] = body.mass * 0.125*math.pi**2 * math.cos(body.t*math.pi)
                for i in range(3):
                    M[i, i] += body.mass
                    for j in range(3):
                        M[3+i, 3+j] += I[i, j]
                r = np.zeros((18,),'d')
                BE=True
                CF=opts.cf
                if BE:
                    r[:6] = np.matmul(M, u[:6]) - np.matmul(body.Aij, body.last_u[:6]) - body.last_mom - DT*body.FT - CF*np.matmul(body.Aij,body.last_a)
                    r[6:9] = h - body.last_h - DT*v
                    rQ = Q - body.last_Q - DT*np.matmul(Omega,Q)
                else:
                    r[:6] = np.matmul(M, u[:6]) - np.matmul(body.Aij, body.last_u[:6]) - body.last_mom - DT*(body.FT*theta+body.last_FT*(1.0-theta)) - CF*np.matmul(body.Aij,body.last_a)
                    r[6:9] = h - body.last_h - DT*0.5*(v + body.last_velocity)
                    rQ = Q - body.last_Q - DT*0.25*np.matmul((Omega + body.last_Omega),(Q+body.last_Q))
                r[9:18] = rQ.flatten()
                J = np.zeros((18,18),'d')
                J[:6,:6] = M
                #neglecting 0:6 dependence on Q
                for i in range(3):
                    if BE:
                        J[6+i,i] = -DT
                    else:
                        J[6+i,i] = -DT*0.5
                    J[6+i,6+i] = 1.0
                for i in range(9):
                    J[9+i,9+i] = 1.0
                for i in range(3):
                    for j in range(3):
                        if BE:
                            J[9+i*3+j, 9+i+j*3] -= DT*Omega[i,j]
                        else:
                            J[9+i*3+j, 9+i+j*3] -= DT*0.25*(Omega+body.last_Omega)[i,j]
                #neglecting 9:18 dependence on omega
                body.Omega[:] = Omega
                body.velocity[:] = v
                body.Q[:] = Q
                body.h[:] = h
                body.u[:] = u
                body.mom[:3] = body.mass*u[:3]
                body.mom[3:6] = np.matmul(I,u[3:6])
                body.a[:] = u[:6] - body.last_u[:6]
                return r, J
            nd = body.nd
            Q_start=body.Q.copy()
            h_start=body.last_h.copy()
            for i in range(int(n)):
                theta = (i+1)*DT/dt
                body.t = body.last_t + DT
                logEvent("6DOF theta "+`theta`)
                u = np.zeros((18,),'d')
                u[:] = body.last_u
                r = np.zeros((18,),'d')
                r,J = F(u,theta)
                its=0
                maxits=100
                while ((its==0 or np.linalg.norm(r) > 1.0e-8) and its < maxits):
                    u -= np.linalg.solve(J,r)
                    r,J = F(u,theta)
                    its+=1
                logEvent("6DOF its "+`its`)
                logEvent("6DOF res "+`np.linalg.norm(r)`)
                body.last_Aij[:]=body.Aij
                body.last_FT[:] = body.FT
                body.last_Omega[:] = body.Omega
                body.last_velocity[:] = body.velocity
                body.last_Q[:] = body.Q
                body.last_h[:] = body.h
                body.last_u[:] = body.u
                body.last_mom[:] = body.mom
                body.last_t = body.t
            body.last_a[:] = body.a
            # translate and rotate
            body.h -= h_start
            body.last_h[:] = 0.0
            body.last_position[:] = body.position
            #body.rotation_matrix[:] = np.linalg.solve(Q_start,body.Q)
            body.rotation_matrix[:] = np.matmul(np.linalg.inv(body.Q),Q_start)
            body.rotation_euler[2] -= math.asin(body.rotation_matrix[1,0])#cek hack
            body.Shape.translate(body.h[:nd])
            body.barycenter[:] = body.Shape.barycenter
            body.position[:] = body.Shape.barycenter
            #logEvent("6DOF its = " + `its` + " residual = "+`r`)
            logEvent("6DOF time "+`body.t`)
            logEvent("6DOF DT "+`DT`)
            logEvent("6DOF n "+`n`)
            logEvent("6DOF force "+`body.FT[1]`)
            logEvent("displacement, h = "+`body.h`)
            logEvent("rotation, Q = "+`body.Q`)
            logEvent("velocity, v = "+`body.velocity`)
            logEvent("angular acceleration matrix, Omega = "+`body.Omega`)

        body.step = step
        #body.scheme = 'Forward_Euler'

# ----- SHAPES ----- #
tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)


# ----- BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setTank()  # sliding mesh nodes
tank.BC['y-'].setTank()  #sliding mesh nodes




domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.genMesh
he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.use_gmsh = opts.use_gmsh
geofile='mesh'+str(opts.he)
domain.geofile=geofile


# MESH REFINEMENT

if opts.use_gmsh:
    import py2gmsh
    from MeshRefinement import geometry_to_gmsh
    mesh = geometry_to_gmsh(domain)
    grading = opts.refinement_grading
    he = opts.he
    he_max = 10.
    ecH = 3.
    field_list = []

    def mesh_grading(start, he, grading):
        return '{he}*{grading}^(1+log((-1/{grading}*(abs({start})-{he})+abs({start}))/{he})/log({grading}))'.format(he=he, start=start, grading=grading)

    me01 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist = '(abs({zcoord}-y))'.format(zcoord=water_level)
    me01.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me01]

    me3 = py2gmsh.Fields.MathEval(mesh=mesh)
    dist_z = '(abs(abs({z_p}-y)+abs(y-{z_n})-({z_p}-{z_n}))/2.)'.format(z_p=max(caisson.vertices[:,1]), z_n=min(caisson.vertices[:,1]))
    dist_x = '(abs(abs({z_p}-x)+abs(x-{z_n})-({z_p}-{z_n}))/2.)'.format(z_p=max(caisson.vertices[:,0]), z_n=min(caisson.vertices[:,0]))
    me3.F = '{he}*{grading}^(Sqrt({dist_x}^2+{dist_z}^2)/{he})'.format(he=he, grading=grading, dist_x=dist_x, dist_z=dist_z)
    me3.F = mesh_grading(he=he, start=dist, grading=grading)
    field_list += [me3]

    # background field
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)

    # max element size
    mesh.Options.Mesh.CharacteristicLengthMax = he_max

    mesh.writeGeo(geofile+'.geo')



# passed in added_mass_p.py coefficients
max_flag = 0
max_flag = max(domain.vertexFlags)
max_flag = max(domain.segmentFlags+[max_flag])
max_flag = max(domain.facetFlags+[max_flag])
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
if opts.use_chrono:
    for s in system.subcomponents:
        if type(s) is crb.ProtChBody:
            for i in range(s.i_start, s.i_end):
                flags_rigidbody[i] = 1
else:
    flags_rigidbody[1:5] = 1


##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = opts.weak_factor/nu_0#Re
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
if opts.sc == 0.25:
    sc = 0.25 # default: 0.5. Test: 0.25
    sc_beta = 1.5 # default: 1.5. Test: 1.
    epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
elif opts.sc == 0.5:
    sc = 0.5
    sc_beta = 1.5
    epsFact_consrv_diffusion = 10.0 # default: 1.0. Test: 0.1. Safe: 10.
else:
    import sys
    sys.quit()
ns_forceStrongDirichlet = opts.strong_dir
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
    rd_shockCapturingFactor  = 0.75
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
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

tolfac = 0.001
mesh_tol = 0.001
ns_nl_atol_res = max(1.0e-8,tolfac*he**2)
vof_nl_atol_res = max(1.0e-8,tolfac*he**2)
ls_nl_atol_res = max(1.0e-8,tolfac*he**2)
mcorr_nl_atol_res = max(1.0e-8,tolfac*he**2)
rd_nl_atol_res = max(1.0e-8,tolfac*he)
kappa_nl_atol_res = max(1.0e-8,tolfac*he**2)
dissipation_nl_atol_res = max(1.0e-8,tolfac*he**2)
mesh_nl_atol_res = max(1.0e-8,mesh_tol*he**2)
am_nl_atol_res = 0.001#max(1.0e-8,mesh_tol*he**2)

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
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
