from math import *
import proteus.MeshTools
from proteus import *
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent


from proteus import Context
#===============================================================================
# Context
#===============================================================================
ct = Context.Options([
    ("T",                   1.0,"Time interval [0, T]"),
    ("Refinement",          5, "Specify initial mesh size by giving number of cells in each direction"),
    ("spaceOrder",          1,"FE space for velocity"),
    ("parallel",            False,"Use parallel or not"),
    ("dt_fixed",            0.001,"fixed time step"),
    ("nDTout",              100,"output number"),
    #########################################
    ("use_supg",            1,"Use supg or not"),
    ("nonconservative",     0,"0=conservative"),
    ("forceStrongDirichlet",False,"strong or weak"),
    ("use_sbm",             False,"use sbm instead of imb"),
    #########################################
    ("hk",1.0,"composite quadrature rule"),
    ("structured",          True,"Use structured mesh"),
    ("use_repulsive_force", False,"use repulsive force from Glowinski paper")
], mutable=True)

#===============================================================================
# use Balls
#===============================================================================
use_ball_as_particle = 1

use_supg = ct.use_supg##########

#===============================================================================
# Use SBM or not
#===============================================================================
if ct.use_sbm:
    USE_SBM=1
else:
    USE_SBM=0


#
#===============================================================================
# Chrono parameters
#===============================================================================
g_chrono = [0.0, -981, 0.0]
g = [0.0, -981, 0.0]

# Domain and mesh
L = [1.0, 1.0, 0.6] #1.0 10.0 0.6
# Initial condition
container_dim=[L[0],L[1],L[2]] #Dimensions of the container (Height/Width/depth)
container_cent=[L[0]/2,L[1]/2,0.0] #Position of the center of the container"

particle_diameter=0.25

particle_density=1.5#times 6 to make it like a cylinder. since it is a ball. 

dT_Chrono=ct.dt_fixed/1000.0
#===============================================================================
# Use balls as particles
#===============================================================================
ball_center = np.array([[0.25,0.6,0.0],
                        ],'d') #0.25,9.6,0.0
nParticle = ball_center.shape[0]
ball_radius = np.array([particle_diameter*0.5]*nParticle,'d')
ball_velocity = np.zeros((nParticle, 3),'d')
ball_angular_velocity = np.zeros((nParticle, 3),'d')

#===============================================================================
# Fluid parameter
#===============================================================================
# Water
rho_0 = 1.0
nu_0 = 0.05


# Air
rho_1 = rho_0
nu_1 = nu_0

# Sediment

rho_s = rho_0
nu_s = nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Initial condition
waterLine_x = 2*L[0]
waterLine_z = 2*L[1]

wall_x = 1.5
Um = 0.0


#===============================================================================
# #  Discretization -- input options
#===============================================================================
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = ct.Refinement
sedimentDynamics=False
genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = True
parallel = ct.parallel
if parallel:
   usePETSc = True
   useSuperlu=False
else:
   usePETSc = False
   
#===============================================================================
# Numerical parameters
#===============================================================================
timeDiscretization = 'be'#vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = ct.spaceOrder
pspaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 0.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega
openTop=True
fl_H = L[1]
# Input checks
if spaceOrder not in [1, 2]:
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

#===============================================================================
# FE space
#===============================================================================
if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        comp_quad = Quadrature.CompositeTriangle(elementQuadrature, ct.hk)
        elementQuadrature = comp_quad
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 5)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 5)

if pspaceOrder == 1:
    if useHex:
        pbasis = C0_AffineLinearOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineLinearOnSimplexWithNodalBasis
elif pspaceOrder == 2:
    if useHex:
        pbasis = C0_AffineLagrangeOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineQuadraticOnSimplexWithNodalBasis

#===============================================================================
# # Domain and mesh
#===============================================================================
he = L[0]/float(4*Refinement-1)
he*=0.5
he*=0.5



weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType =proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0




# new style PointGauges
# pointGauges = PointGauges(gauges = ((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))), (('p',), ((0.5, 0.5, 0),))),
#                           activeTime=(0, 0.5),
#                           sampleRate=0,
#                           fileName='combined_gauge_0_0.5_sample_all.csv')

# lineGauges = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825': ((0.495, 0.0, 0.0), (0.495, 1.8, 0.0))}, linePoints=20)
# #'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
# lineGauges_phi = LineGauges_phi(lineGauges.endpoints, linePoints=20)

structured = ct.structured


if useHex:
    nnx = 4 * Refinement + 1
    nny = 2 * Refinement + 1
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        N = 2**Refinement
        nnx = int(L[0])*N
        nny = int(L[1])*N
        triangleFlag=1
        domain = Domain.RectangularDomain(L=[L[0],L[1]])
    else:
        vertices = [[0.0, 0.0],  #0
                    [L[0], 0.0],  #1
                    [L[0], L[1]],  #2
                    [0.0, L[1]],  #3
#                     [0.2-0.16,L[1]*0.2],
#                     [0.2-0.16,L[1]*0.8],
#                     [0.2+0.3,L[1]*0.8],
#                     [0.2+0.3,L[1]*0.2],
#                     # the following are set for refining the mesh
#                     [0.2-0.06,0.2-0.06],
#                     [0.2-0.06,0.2+0.06],
#                     [0.2+0.06,0.2+0.06],
#                     [0.2+0.06,0.2-0.06],
                    ]

                    
                    
        vertexFlags = [boundaryTags['bottom'],
                       boundaryTags['bottom'],
                       boundaryTags['top'],
                       boundaryTags['top'],
                       # the interior vertices should be flaged to 0
#                        0, 0, 0, 0,
#                        0, 0, 0, 0,
                       ]

        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0],
                    #Interior segments
#                     [4, 5],
#                     [5, 6],
#                     [6, 7],
#                     [7,4],
#                     [8,9],
#                     [9,10],
#                     [10,11],
#                     [11,8],
                    ]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left'],
#                         0,
#                         0,
#                         0,
#                         0,
#                         0,
#                         0,
#                         0,
#                         0,
                        ]

        regions = [[0.95*L[0], 0.2*L[1]],
#                    [0.2-0.15,0.2],[0.2,0.2],
                   ]
        regionFlags = [1,
#                        2,3,
                       ]
        regionConstraints=[0.5*he**2,
#                            0.5*(he/2.0)**2,0.5*(he/6.0)**2,
                           ]
        #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
        #            vertices.append(gaugeCoordinates)
        #            vertexFlags.append(pointGauges.flags[gaugeName])

        # for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
        #     for gaugeCoordinates in gaugeLines:
        #         vertices.append(gaugeCoordinates)
        #         vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags,
                                                      regionConstraints=regionConstraints)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        #triangleOptions = "VApq30Dena%8.8f" % ((he ** 2) / 2.0,)
        triangleOptions = "VApq30Dena"
        logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))

#===============================================================================
# Time stepping
#===============================================================================
T=ct.T
dt_fixed = ct.dt_fixed#0.03
dt_init = 0.5*dt_fixed#########HUGH EFFECT FOR EXTRAPOLATION
runCFL=0.3
nDTout = ct.nDTout
#nDTout = int(T/dt_fixed)
dt_output = T/nDTout
if dt_init < dt_fixed:
    tnList = [0.0,dt_init]+[i*dt_output for i in range(1,nDTout+1)]
else:
    tnList = [0.0]+[i*dt_output for i in range(1,nDTout+1)]

#===============================================================================
# Numerical parameters
#===============================================================================
ns_forceStrongDirichlet = True
ns_sed_forceStrongDirichlet = False
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
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.9
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor = 0.9
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_vos = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

ns_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
ns_sed_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vos_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.05 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
phi_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)

#turbulence
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
def velRamp(t):
    if t < 0.0:
        return Um*exp(min(0.0,100*(t-0.25)/0.25))
    else:
        return Um



#===============================================================================
# Glowinski parameters
#===============================================================================
force_range = 1.5*he
stiffness_particle = 1.0e-5
stiffness_wall = 0.5e-5
#===============================================================================
# AuxiliaryVariables
#===============================================================================
import Chrono
class ChronoModel(AuxiliaryVariables.AV_base):
    def __init__(self,timeStep=1e-3,
                 m_container_center=container_cent,
                 m_container_dims=container_dim,
                 m_gravity=g_chrono,
                 m_particles_diameter=particle_diameter,
                 m_particles_density=particle_density,
                 dt_init=dt_init):
        self.mtime=0
        self.dt_init=dt_init
        self.chmodel = Chrono.MBDModel(
            m_timeStep=timeStep,
            m_container_center=np.array(m_container_center,dtype="d"),
            m_container_dims=np.array(m_container_dims,dtype="d"),
            m_particles_diameter=particle_diameter,
            m_particles_density=particle_density,
            m_gravity=np.array(m_gravity,dtype="d"),
            nParticle = nParticle,
            ball_center=ball_center,
            ball_radius=ball_radius)
        self.nnodes = self.chmodel.getNumParticles()
        self.solidPosition = np.zeros((self.nnodes,3), 'd')
        self.solidVelocity = np.zeros((self.nnodes,3), 'd')
        self.solidAngularVelocity = np.zeros((self.nnodes,3), 'd')
        self.solidPosition_old = np.zeros((self.nnodes,3), 'd')
        self.solidVelocity_old = np.zeros((self.nnodes,3), 'd')
        self.solidAngularVelocity_old = np.zeros((self.nnodes,3), 'd')

        for i in range(self.nnodes):
            self.solidPosition[i,:] = self.chmodel.get_Pos(i)
            self.solidPosition_old[i,:] = self.chmodel.get_Pos(i)
            #initial velocity and initial angular velocity are 0
        self.proteus_dt = self.dt_init

        self.solidForces = np.zeros((nParticle,3),'d')

        #rkpm
        #now add the Python code that generates the list of rkpm nodes here
        #(OUTPUT : n_rkpm_nodes)

#        particle_diameter=0.25
        r=particle_diameter/2
        #nodes_th=36
        #nodes_r=10
        #n_rkpm_nodes = nodes_th*nodes_r+1
        #self.n_rkpm_nodes = n_rkpm_nodes
        #rkpm_nodes = np.zeros((n_rkpm_nodes,2)) #node_x=[x;y]
        
        
        nodes_r = 3
        dr = r / nodes_r
        dR = np.linspace(dr, r, nodes_r)
        nR = len(dR)
        nT = np.zeros(nR)

        for i in range (nR):
            nT[i] = int(10*(i+1))
        
        n_rkpm_nodes = int(sum(nT)+1)
        self.n_rkpm_nodes = n_rkpm_nodes
        
        rkpm_nodes = np.zeros((n_rkpm_nodes,2)) #node_x=[x;y]
        self.rkpm_force = np.zeros((n_rkpm_nodes,3))
        #rkpm_nodes[0][0] = center_x
        #rkpm_nodes[1][0] = center_y
        
        #PI = np.pi
        #dtheta = np.linspace(0.0, 360, nodes_th)
        #dR=np.linspace(dr, r, nodes_r)

        ###### need center of the particle
        rkpm_nodes[0][0] = center_x = ball_center[0,0]
        rkpm_nodes[0][1] = center_y = ball_center[0,1]

        def rtpairs(r, n):
            
            for i in range(len(r)):
                for j in range(int(n[i])):    
                    yield r[i], j*(2 * np.pi / n[i])
        #------------------------------------------------------------------
        #RKPM shape function (INPUT: X1(Node),GCOO(Points); OUTPUT: SHP)
        
        c = 1
        for r, t in rtpairs(dR, nT):
            rkpm_nodes[c][0]=r * np.cos(t) + center_x #X
            rkpm_nodes[c][1]=r * np.sin(t) + center_y #Y
            c+=1
#        for drr in dR:
#            for theta in dtheta:
#                rkpm_nodes[c][0] = drr * np.cos(theta/180*PI) + center_x #X
#                rkpm_nodes[c][1] = drr * np.sin(theta/180*PI) + center_y #Y
#                c += 1
        import matplotlib.pyplot as plt
        plt.scatter(rkpm_nodes[:,0],rkpm_nodes[:,1])
        self.rkpm_nodes = rkpm_nodes
        plt.savefig("rkpm_nodes.png")
        

    def RKPM_shape(self, X1,GCOO):
        #DEG = 1 #linear basis
        MSIZE = 3 # linear basis function [1,x,y]
        CONT = 3 # C3 conti
        support = 3
        h = 0.25/2/10 #dr #nodal spacing
        dilation = h * support

        XMXI_OA = (X1[0] - GCOO[0]) / dilation
        YMXI_OA = (X1[1] - GCOO[1]) / dilation
        #            XMXI_OA[I] = (X1[0] - GCOO[0][I]) / dilation
        #            YMXI_OA[I] = (X1[1] - GCOO[1][I]) / dilation
        
        # initialization
        H_FULL = np.zeros(MSIZE)
        #XMXI_OA = np.zeros(N)
        #YMXI_OA = np.zeros(N)
        M_FULL = np.zeros((MSIZE,MSIZE))
        PHI = np.zeros(N)
        H0 = np.zeros(MSIZE)
        SPH = np.zeros(N)
        B = np.zeros(MSIZE)

        H_FULL[0] = 1.0
        H_FULL[1] = XMXI_OA#[I]
        H_FULL[2] = YMXI_OA#[I]

        # kernel function
        (PHIX) = self.MLS_KERNEL(abs(XMXI_OA),dilation,CONT)
        (PHIY) = self.MLS_KERNEL(abs(YMXI_OA),dilation,CONT)
        PHI = PHIX * PHIY
        #PHI[I] = PHIX * PHIY

        for J in range(MSIZE):
            for K in range(MSIZE):
                M_FULL[J,K] += H_FULL[J] * H_FULL[K] * PHI

        MINV = np.linalg.pinv(M_FULL)
        H0[0] = 1.0
        B = np.matmul(H0,MINV)

        #            for I in range (N): #loop over nodes
        #                H_FULL[0] = 1
        #                H_FULL[1] = XMXI_OA#[I]
        #                H_FULL[2] = YMXI_OA#[I]
        C =np.matmul(B,H_FULL)
        SHP = C * PHI #RK shape function
        #SHP[I] = C * PHI[I] #RK shape function
        return (SHP)

    def MLS_KERNEL(self,XSA,AJ,ISPLINE):
        if XSA <0.5:
           PHI = 2.0/3 - 4*XSA**2 - 4*XSA**3
        elif XSA <1:
           PHI = 4.0/3 -4*XSA + 4*XSA**2 -4.0/3*XSA**3
        else:
           PHI = 0

        return (PHI)

#    def neighbor_list(self,XCOO,YCOO, xg, yg, dm):
#        # get a neighbor list and a Length ID
#        #XCOO,YCOO: Nodes
#        #xg,yg: Points
#        #dm: dilation
#        
#        np_x = len[xg]
#        np_y = len[yg]
#        nnodes_x = len[XCOO]
#        nnodes_y = len[YCOO]
#        nnodes=nnodes_x*nnodes_y
#        np=np_x*np_y
#        neighor_list=np.zeros[nnodes,50]
#        LengthID=np.zeros[1,nnodes]
#
#        for j in range(nnodes):
#            m=0
#            for i in range(np):
#                rx  = abs(XCOO[j]-xg[i])/dm
#                ry  = abs(YCOO[j]-yg[i])/dm
#                if (rx < 1.0 & ry < 1.0 ):
#                    m = m+1
#                    neighor_list[j][m]= i
#                else:
#                    pass
#            LengthID[j]= m
#        return (neighor_list, LengthID)
        #self.get_rkpm_v(node_x, x) -> v_n(x)

    def attachModel(self,model,ar):
        self.chmodel.attachModel(model,ar)
        self.model=model
        self.rkpm_test = np.zeros((self.n_rkpm_nodes,)+ self.model.levelModelList[-1].q['x'].shape[:2],'d')
        self.model.levelModelList[-1].coefficients.rkpm_test = self.rkpm_test
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin >= 0)
        assert(flagMax <= 7)
        self.nForces=flagMax+1
        assert(self.nForces <= 8)
        #initialize RKPM test functions
        #m.q['x'] is a nElements x nQuad x nd array of points
        #we want rkpm_test[i,eN,k] to be value of i-th test function at
        #quadrature point eN,k with physical location m.q['x'][eN,k]
        import matplotlib.pyplot as plt
        for i in range(self.n_rkpm_nodes):
            for eN in range(m.q['x'].shape[0]):
                for k in range(m.q['x'].shape[1]):
                    #RKPM_shape(X1,GCOO)
#                    self.rkpm_test[i,eN,k] = self.get_rkpm_v(self.rkpm_nodes[:][i], 
#                                  m.q['x'][eN, k])
                   self.rkpm_test[i,eN,k] = self.RKPM_shape(self.rkpm_nodes[i][:], m.q['x'][eN, k])
            fig = plt.figure()
            plt.scatter(self.rkpm_nodes[:,0],self.rkpm_nodes[:,1],marker='x')
            #if self.rkpm_test[i]>0.0:
            plt.scatter(m.q['x'][...,0],m.q['x'][...,1], c=self.rkpm_test[i,...],marker='o',s=10)
            plt.colorbar()  # show color scale
            #plt.scatter(m.q['x'][...,0],m.q['x'][...,1], c=self.rkpm_test[i,...],marker='o')
            plt.savefig("rkpm_test_i_{0}.png".format(i))

        #check consistency
        #assert abs(sum(self.rkpm_test[i])-1.0) > 1.0e-5, 'The sum of shape functions is not the 1'
        
        return self

    def get_u(self):
        return 0.
    def get_v(self):
        return 0.
    def get_w(self):
        return 0.
    def calculate_init(self):
        self.last_F = None
        self.calculate()

    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        self.solidForces[:] = self.model.levelModelList[-1].coefficients.particle_netForces[:nParticle,:]
        ######## Test 4: implement repulsive force (19) and (20) in Glowinski's paper; Use `friction_const=0` in Chrono.h
        if ct.use_repulsive_force:
            import RepulsiveForce as RF
            for i in range(nParticle):
                for j in range(nParticle):
                    self.solidForces[i,:] += RF.get_repulsive_force_glowinski(self.solidPosition[i,:], 0.5*particle_diameter,
                                                                              self.solidPosition[j,:], 0.5*particle_diameter,
                                                                              force_range, stiffness_particle)
                self.solidForces[i,:] += RF.get_repulsive_force_from_vetical_wall_glowinski(self.solidPosition[i,:], 0.5*particle_diameter,
                                                                              0.0,force_range, stiffness_wall)
                self.solidForces[i,:] += RF.get_repulsive_force_from_vetical_wall_glowinski(self.solidPosition[i,:], 0.5*particle_diameter,
                                                                              L[0],force_range, stiffness_wall)
                self.solidForces[i,:] += RF.get_repulsive_force_from_horizontal_wall_glowinski(self.solidPosition[i,:], 0.5*particle_diameter,
                                                                              0.0,force_range, stiffness_wall)

        ######## Test 2: Turn on rotation in solid solver
        self.solidMoments=self.model.levelModelList[-1].coefficients.particle_netMoments
#         self.solidMoments = 0.0*self.model.levelModelList[-1].coefficients.particle_netMoments
        try:
            self.proteus_dt = self.model.levelModelList[-1].dt_last
            self.mtime = self.model.stepController.t_model_last
        except:
            self.proteus_dt = self.dt_init
            self.mtime = None
            return

        print("time/dt before chrono calculation:" , self.mtime , self.proteus_dt)
        
        ######### Test 1: use 0 force to make the solid fall down straightly
#         self.solidForces[:] = 0.0
#         self.solidForces[0][1] = 0.0

        self.chmodel.calculate(self.solidForces,self.solidMoments,self.proteus_dt)
            # self.chplate.calcNodalInfo()
        
        self.solidPosition_old[:] = self.solidPosition
        self.solidVelocity_old[:] = self.solidVelocity
        self.solidAngularVelocity_old[:] = self.solidAngularVelocity
        for i in range(self.nnodes):
            self.solidPosition[i,:] = self.chmodel.get_Pos(i)
            self.solidVelocity[i,:] = self.chmodel.get_Vel(i)
            self.solidAngularVelocity[i,:] = self.chmodel.get_Angular_Vel(i)
        
        ######## update proteus
        self.model.levelModelList[-1].coefficients.ball_center[:] = self.solidPosition
        self.model.levelModelList[-1].coefficients.ball_velocity[:] = self.solidVelocity
#         self.model.levelModelList[-1].coefficients.ball_velocity[:] = (self.solidPosition-self.solidPosition_old)/myChModel.proteus_dt
        ######## Test 3: Turn on rotation in fluid solver
        self.model.levelModelList[-1].coefficients.ball_angular_velocity[:] = self.solidAngularVelocity ##better result
#         self.model.levelModelList[-1].coefficients.ball_angular_velocity = 0.5*(self.solidAngularVelocity+self.solidAngularVelocity_old)
        
        # for writeTime in tnList:
        #     if (self.mtime>0.0001):
        #         if (abs(self.mtime-writeTime)<0.0001):
        #             self.chmodel.writeFrame()
        #
        #write Python code to evaluate RKPM test functions at quadrature points 
        #
myChModel = ChronoModel(timeStep=dT_Chrono,
                        m_container_center=container_cent,
                        m_container_dims=container_dim,
                        m_particles_diameter=particle_diameter,
                        m_particles_density=particle_density,
                        m_gravity=(g_chrono[0],g_chrono[1],g_chrono[2]))

numPars=myChModel.chmodel.getNumParticles()
print(numPars)
