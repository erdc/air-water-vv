"""
Bar floating in half-filled tank
"""
from math import *
from proteus import *
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


from proteus import Context
opts=Context.Options([
    ("bar_dim", (0.33, 0.33, 0.2), "Dimensions of the bar"),
    ("scale" , 21, "Scaling Factor"),
    # ("tank_dim", (1.0,1.0,1.0), "Dimensions of the tank"),
    ("topTank", (1.0,1.0), "Top Rectangle of the tank"),
    ("heightTank", 1.0, "Height of the tank"),
    ("baseTank", (0.5,0.5), "Bottom Rectangle of the tank"),
    ("culvertDim",(0.2,0.25,0.3),"Dimensions of the culvert"),
    ("water_surface_height",0.3,"Height of free surface above bottom"),
    ("speed",0.0,"Speed of current if non-zero"),
    ("bar_height",0.5,"Initial height of bar center above bottom"),
    ("bar_rotation",(0.01,0,0),"Initial rotation about x,y,z axes"),
    ("refinement_level",0,"Set maximum element diameter to he/2**refinement_level"),
    ("gen_mesh",True,"Generate new mesh"),
    ("he",0.1,"Mesh size"),
    ("T",10.0,"Simulation time"),
    ("dt_init",1e-6,"Initial time step"),
    ("cfl",0.33,"Target cfl"),
    ("nsave",100,"Number of time steps to  save"),
    ("parallel",True,"Run in parallel"),
    ("free_x",(0.0,0.0,1.0),"Free translations"),
    ("free_r",(1.0,1.0,0.0),"Free rotations")])


#----------------------------------------------------
# Bridge properties
#----------------------------------------------------


#----------------------------------------------------
# Physical properties
#----------------------------------------------------
rho_0=998.2
nu_0 =1.004e-6

rho_1=1.205
nu_1 =1.500e-5

sigma_01=0.0

g=[0.0,0.0,-9.81]

#----------------------------------------------------
# Domain - mesh - quadrature
#----------------------------------------------------
nd = 3


waterLevel   =  opts.water_surface_height



nLevels = 1

he = opts.he
he *=(0.5)**opts.refinement_level
genMesh=opts.gen_mesh

boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':7, 'sponge':8}

bTdim=opts.baseTank
tTdim=opts.topTank
hTank=opts.heightTank

L=(1.0,1.0,1.0)   # Temporal Fix, needs to be eliminated and all dependents changed.

oD = (0.5,0.5,0.0)          #Origin of the Domain

cDim = opts.culvertDim
ct_inter=0.5*(cDim[2]/hTank)*(tTdim[1]-bTdim[1])

#tank
vertices=[#Bottom
            [oD[0]-bTdim[0]/2,oD[1]-bTdim[1]/2,oD[2]],              #0,(0,0,0)
            [oD[0]+bTdim[0]/2,oD[1]-bTdim[1]/2,oD[2]],              #1,(1,0,0)
            [oD[0]+bTdim[0]/2,oD[1]+bTdim[1]/2,oD[2]],              #2,(1,1,0)
            [oD[0]-bTdim[0]/2,oD[1]+bTdim[1]/2,oD[2]],              #3,(0,1,0)
            #Top
            [oD[0]-tTdim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank],        #4,(0,0,1)
            [oD[0]+tTdim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank],        #5,(1,0,1)
            [oD[0]+tTdim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank],        #6,(1,1,1)
            [oD[0]-tTdim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank],        #7,(0,1,1)
            #Culvert Source (-Y, Back) Top
            [oD[0]-cDim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank],         #8,(0,0,1)
            [oD[0]+cDim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank],         #9,(1,0,1)
            [oD[0]+cDim[0]/2,oD[1]-tTdim[1]/2+cDim[1],oD[2]+hTank], #10,(1,1,1)
            [oD[0]-cDim[0]/2,oD[1]-tTdim[1]/2+cDim[1],oD[2]+hTank], #11,(0,1,1)
            # Culvert Sink (+Y, Front) Top
            [oD[0]-cDim[0]/2,oD[1]+tTdim[1]/2-cDim[1],oD[2]+hTank], #12,(0,0,1)
            [oD[0]+cDim[0]/2,oD[1]+tTdim[1]/2-cDim[1],oD[2]+hTank], #13,(1,0,1)
            [oD[0]+cDim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank],         #14,(1,1,1)
            [oD[0]-cDim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank],         #15,(0,1,1),
            # Culvert Source (-Y, BAck) Bottom
            [oD[0]-cDim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank-cDim[2]],      #16, (0,0,0)
            [oD[0]+cDim[0]/2,oD[1]-tTdim[1]/2,oD[2]+hTank-cDim[2]],      #17, (1,0,0)
            [oD[0]+cDim[0]/2,oD[1]-tTdim[1]/2+cDim[1],oD[2]+hTank-cDim[2]],  #18, (1,1,0)
            [oD[0]-cDim[0]/2,oD[1]-tTdim[1]/2+cDim[1],oD[2]+hTank-cDim[2]],  #19, (0,1,0)
            # Culvert Sink (+Y, Front) Bottom
            [oD[0]-cDim[0]/2,oD[1]+tTdim[1]/2-cDim[1],oD[2]+hTank-cDim[2]],          #20, (0,0,0)
            [oD[0]+cDim[0]/2,oD[1]+tTdim[1]/2-cDim[1],oD[2]+hTank-cDim[2]],          #21, (1,0,0)
            [oD[0]+cDim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank-cDim[2]],                 #22, (1,1,0)
            [oD[0]-cDim[0]/2,oD[1]+tTdim[1]/2,oD[2]+hTank-cDim[2]],                 #23, (0,1,0)
            # Intersection Points (-Y, back)
            [oD[0]-cDim[0]/2, oD[1]-tTdim[1]/2+ct_inter, oD[2]+hTank-cDim[2]],  # 24, (1,0.5,0)
            [oD[0]+cDim[0]/2, oD[1]-tTdim[1]/2+ct_inter, oD[2]+hTank-cDim[2]],  # 25, (0,0.5,0)
            # Intersection Points (+Y, Front)
            [oD[0]+cDim[0]/2, oD[1]+tTdim[1]/2-ct_inter, oD[2]+hTank-cDim[2]],  # 26, (1,0.5,0)
            [oD[0]-cDim[0]/2, oD[1]+tTdim[1]/2-ct_inter, oD[2]+hTank-cDim[2]],  # 27, (0,0.5,0)
    
]


vertexFlags=[boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left']]

nVert=len(vertices)
# add sponge tag to the rest of the vertices.
for i in range(nVert-len(vertexFlags)):
    vertexFlags.append(boundaryTags['sponge'])

facets=[[[0,1,2,3]],                        #0   Bottom
        [[4,8,11,10,9,5,6,14,13,12,15,7]],  #1   Top H-shape
        [[0,3,7,4]],                        #2   Left
        [[1,2,6,5]],                        #3   Right
        [[4,24,8]],                         #4   Back top-left triangle
        [[4,0,24]],                         #5   Back top-bottom triangle
        [[0,1,25,24]],                      #6   Back center square
        [[1,5,25]],                         #7   Back bottom-right triangle
        [[5,9,25]],                         #8   Back top-right     triangle
        [[6,26,14]],                        #9   Front top-left triangle
        [[6,2,26]],                         #10  Front top-bottom triangle
        [[2,3,27,26]],                      #11  Front center square
        [[3,7,27]],                         #12  Front bottom-right triangle
        [[7,15,27]],                        #13  Front top-right     triangle
        [[8,9,10,11]],                      #14  Source Culvert Top
        [[16,17,9,8]],                      #15  Source Culvert Out Wall
        [[19,18,10,11]],                    #16  Source Culvert In Wall
        [[16,24,8]],                        #17  Source Culvert Left-Out
        [[24,19,9,8]],                      #18  Source Culvert Left-In
        [[17,25,9]],                        #19  Source Culvert Right-Out
        [[25,18,10,9]],                     #20  Source Culvert Right-In
        [[16,17,25,24]],                    #21  Source Culvert Bottom-Out
        [[25,18,19,24]],                    #22  Source Culvert Bottom-In
        [[12,13,14,15]],                    #23  Sink Culvert Top
        [[20,21,13,12]],                    #24  Sink Culvert In Wall
        [[23,22,14,15]],                    #25  Sink Culvert Out Wall
        [[20,27,15,12]],                    #26  Sink Culvert Left-In
        [[27,23,15]],                       #27  Sink Culvert Left-Out
        [[21,26,14,13]],                    #28  Sink Culvert Right-In
        [[26,22,14]],                       #29  Sink Culvert Right-Out
        [[20,27,26,21]],                    #30  Sink Culvert Bottom-In
        [[27,23,22,26]]                     #31  Sink Culvert Bottom-Out
        ]

facetFlags=[boundaryTags['bottom'],
            boundaryTags['top'],
            boundaryTags['left'],
            boundaryTags['right'],
            ]
for i in range(5):
    facetFlags.append(boundaryTags['back'])

for i in range(5):
    facetFlags.append(boundaryTags['front'])

for i in range(9):
    facetFlags.append(boundaryTags['sponge'])       #Sorce Sponge

for i in range(9):
    facetFlags.append(boundaryTags['sponge'])       #Sink Sponge

regions=[[oD[0],oD[1]-tTdim[1]/2+cDim[1]/2,oD[2]+hTank-cDim[2]/2]]
regions.append([oD[0],oD[1],oD[2]+0.5*hTank])
regions.append([oD[0],oD[1]+tTdim[1]/2-cDim[1]/2,oD[2]+hTank-cDim[2]/2])
regionFlags=[1.0,2.0,3.0]
holes=[]


#################################################
######              Caisson            #############
####################################################


# BAR Properties
he_min = 0.01
scale = 1.0 / opts.scale
caisson_width = 0.7*18.5625*25*0.0254  # 1:25 (in) -> 1:1 (in) -> 1:1 (mts)

#bar_center=[0.5*L[0],0.5*L[1],waterLevel]
bar_center = (oD[0],oD[1],oD[2]+0.5*opts.bar_height)
#Calculate Density
bar_mass = 6350     #Kg taken from datasheet
bar_dim=(6.9,8.5,0.86)    #Dimensions based on SolidWorks model, only used for density calculation
bar_volume = bar_dim[0]*bar_dim[1]*bar_dim[2] #m^3 scaled
bar_density=bar_mass/bar_volume

#Scaled Dimensions for ODE simulation
bar_dim=(round(6.9*scale,3),round(8.5*scale,3),round(0.86*scale,4))


#print opts.bar_dim, bar_dim

speed = opts.speed

#set up barycenters for force calculation
barycenters = numpy.zeros((8,3),'d')
barycenters[7,:] = bar_center

bar_inertia = [[(L[1]**2+L[2]**2)/12.0, 0.0                    , 0.0                   ],
               [0.0                   , (L[0]**2+L[2]**2)/12.0 , 0.0                   ],
               [0.0                   , 0.0                    , (L[0]**2+L[1]**2)/12.0]]

RBR_linCons  = [1,1,0]
RBR_angCons  = [1,0,1]



import pandas as pd

xl = pd.ExcelFile("Ribbon Bridge Section Line Segment_short.xlsx")
df = xl.parse(0)
y = np.asarray(df['x'].tolist())
z = np.asarray(df['y'].tolist())
# Offset to (0,0)

y -= 0.5 * (max(y) + min(y))
z -= 0.5 * (max(z) + min(z))

# Real dimension
prescale = 25 * 0.0254  # 1:25 in -> 1:1 in -> 1:1 mts
y = y * prescale
z = z * prescale
# Scale to experimental model
yp = y * scale      # 1:1 -> 1:scale
zp = z * scale
caisson_width *= scale



#dim = [caisson_width, max(yp) - min(yp), max(zp) - min(zp)]

# Interpolate points to he_caisson
yd = np.diff(yp)
zd = np.diff(zp)
dist = np.sqrt(yd ** 2 + zd ** 2)
u = np.cumsum(dist)
u = np.hstack([[0], u])

t = np.linspace(0, u.max(), int(u.max() / he_min))
yBase = np.interp(t, u, yp)
zBase = np.interp(t, u, zp)

yTop = np.linspace(zp.min() + he_min, zp.max(), int((zp.max() - zp.min()) / he_min), endpoint=False)
zTop = np.full(len(yTop), zp[-1])
y1 = np.hstack([yBase, yTop])
z1 = np.hstack([zBase, zTop])

n1 = len(y1)

# Extrude for 3d object

x1 = np.empty(n1)
x1.fill(caisson_width / 2)

x2 = -x1
y2 = y1
z2 = z1

x = np.hstack([x1, x2])
y = np.hstack([y1, y2])
z = np.hstack([z1, z2])


x += bar_center[0]
y += bar_center[1]
z += bar_center[2]

# print "original vertices"
# print ('\n'.join('{}:, {}'.format(*ver) for ver in enumerate(vertices)))

nStart = len(vertices)
for i in range(len(x)):
    vertices.append([x[i],y[i],z[i]])
    vertexFlags.append(boundaryTags['obstacle'])

# print "added vertices"
# print ('\n'.join('{}:, {}'.format(*ver) for ver in enumerate(vertices[nStart:], start=nStart)))


facet1=[]
for i in range(len(x1)):
    facet1 += [nStart + i]
facets.append([facet1])
facetFlags.append(boundaryTags['obstacle'])
facet2=[]
for i in range(len(x2)):
    facet2 += [nStart + n1 + i]
facets.append([facet2])
facetFlags.append(boundaryTags['obstacle'])


for i in range(n1-1):
    facets.append([[nStart+i , nStart+i+1, nStart+n1+i+1, nStart+n1+i]])
    facetFlags.append(boundaryTags['obstacle'])
facets.append([[nStart+i+1 , nStart, nStart+n1, nStart+n1+i+1]])
facetFlags.append(boundaryTags['obstacle'])
#todo, add initial rotation of bar
# print "facets"
# print ('\n'.join('{}:, {}'.format(*ver) for ver in enumerate(facets)))

holes.append(bar_center)

domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             regions=regions,
                                             regionFlags=regionFlags,
                                             holes=holes)


#go ahead and add a boundary tags member
domain.boundaryTags = boundaryTags
from proteus import Comm
comm = Comm.get()
if comm.isMaster():
    domain.writePoly("mesh")
else:
    domain.polyfile="mesh"
comm.barrier()
triangleOptions="VApq1.35q12feena%21.16e" % ((he**3)/6.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

quad_order = 3

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = False
openSides = False
openEnd = True
smoothBottom = False
smoothObstacle = False
movingDomain=True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init=opts.dt_init
T = opts.T
nDTout=opts.nsave
dt_out =  (T-dt_init)/nDTout
runCFL = opts.cfl

#----------------------------------------------------
water_depth  = waterLevel-oD[2]

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not opts.parallel
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
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
nd = 3
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
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
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

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*he**2)
rd_nl_atol_res = max(1.0e-12,0.01*he)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x,t):
    p_L = 0.0
    phi_L = L[2] - waterLevel
    phi = x[2] - waterLevel
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

import ode

def near_callback(args, geom1, geom2):
    """Callback function for the collide() method.

    This function checks if the given geoms do collide and
    creates contact joints if they do.
    """

    # Check if the objects do collide
    contacts = ode.collide(geom1, geom2)

    # Create contact joints
    world,contactgroup = args
    for c in contacts:
        c.setBounce(0.2)
        c.setMu(5000)
        j = ode.ContactJoint(world, contactgroup, c)
        j.attach(geom1.getBody(), geom2.getBody())

class RigidBar(AuxiliaryVariables.AV_base):
    def __init__(self,density=1.0,bar_center=(0.0,0.0,0.0),bar_dim=(1.0,1.0,1.0),barycenters=None,he=1.0,cfl_target=0.9,dt_init=0.0001):
        self.dt_init = dt_init
        self.he=he
        self.cfl_target=cfl_target
        self.world = ode.World()
        #self.world.setERP(0.8)
        #self.world.setCFM(1E-5)
        self.world.setGravity(g)
        self.g = np.array([g[0],g[1],g[2]])
        self.space = ode.Space()
        eps_x = L[0]- 0.75*L[0]
        eps_y = L[1]- 0.75*L[1]

        self.M = ode.Mass()
        self.totalMass = density*bar_dim[0]*bar_dim[1]*bar_dim[2]
        self.M.setBox(density,bar_dim[0],bar_dim[1],bar_dim[2])
        #bar body
        self.body = ode.Body(self.world)
        self.body.setMass(self.M)
        self.body.setFiniteRotationMode(1)
        #bar geometry
        self.bar = ode.GeomBox(self.space,bar_dim)
        self.bar.setBody(self.body)
        self.bar.setPosition(bar_center)
        self.boxsize = (bar_dim[0],bar_dim[1],bar_dim[2])
        #contact joints
        self.contactgroup = ode.JointGroup()
        self.last_position=bar_center
        self.position=bar_center
        self.last_velocity=(0.0,0.0,0.0)
        self.velocity=(0.0,0.0,0.0)
        self.h=(0.0,0.0,0.0)
        self.rotation = np.eye(3)
        self.last_rotation = np.eye(3)
        self.last_rotation_inv = np.eye(3)
        self.barycenters=barycenters
        self.init=True
        self.bar_dim = bar_dim
        self.last_F = np.zeros(3,'d')
        self.last_M = np.zeros(3,'d')

    def attachModel(self,model,ar):
        self.model=model
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin >= 0)
        assert(flagMax <= 8)
        self.nForces=flagMax+1
        assert(self.nForces <= 9)
        return self
    def get_u(self):
        return self.last_velocity[0]
    def get_v(self):
        return self.last_velocity[1]
    def get_w(self):
        return self.last_velocity[2]
    def calculate_init(self):
        self.last_F = None
        self.calculate()
    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        F = self.model.levelModelList[-1].coefficients.netForces_p[7,:] + self.model.levelModelList[-1].coefficients.netForces_v[7,:];
        M = self.model.levelModelList[-1].coefficients.netMoments[7,:]
        logEvent("x Force " +`self.model.stepController.t_model_last`+" "+`F[0]`)
        logEvent("y Force " +`self.model.stepController.t_model_last`+" "+`F[1]`)
        logEvent("z Force " +`self.model.stepController.t_model_last`+" "+`F[2]`)
        logEvent("x Moment " +`self.model.stepController.t_model_last`+" "+`M[0]`)
        logEvent("y Moment " +`self.model.stepController.t_model_last`+" "+`M[1]`)
        logEvent("z Moment " +`self.model.stepController.t_model_last`+" "+`M[2]`)
        logEvent("dt " +`dt`)
        scriptMotion=False
        linearOnly=False
        if self.last_F == None:
            self.last_F = F.copy()
        if scriptMotion:
            velocity = np.array((0.0,0.3/1.0,0.0))
            logEvent("script pos="+`(np.array(self.position)+velocity*dt).tolist()`)
            self.body.setPosition((np.array(self.position)+velocity*dt).tolist())
            self.body.setLinearVel(velocity)
        else:
            if linearOnly:
                Fstar = 0.5*(F+self.last_F) + np.array(self.world.getGravity())
                velocity_last = np.array(self.velocity)
                velocity = velocity_last + Fstar*(dt/self.totalMass)
                velocity[0] = 0.0
                vmax = self.he*self.cfl_target/dt
                vnorm = np.linalg.norm(velocity,ord=2)
                if vnorm > vmax:
                    velocity *= vmax/vnorm
                    logEvent("Warning: limiting rigid body velocity from "+`vnorm`+" to "+`vmax`)
                position_last = np.array(self.position)
                position = position_last + 0.5*(velocity_last + velocity)*dt
                self.body.setPosition(position.tolist())
                self.body.setLinearVel(velocity.tolist())
                msg = """
Fstar         = {0}
F             = {1}
F_last        = {2}
dt            = {3:f}
velocity      = {4}
velocity_last = {5}
position      = {6}
position_last = {7}""".format(Fstar,F,self.last_F,dt,velocity,velocity_last,position,position_last)
                logEvent(msg)
            else:
                nSteps=10
                # vnorm = np.linalg.norm((F+self.g)*dt)
                # vmax = self.he*self.cfl_target/dt
                # if vnorm > vmax:
                #     F *= vmax/vnorm
                #     logEvent("Warning: limiting rigid body velocity from "+`vnorm`+" to "+`vmax`)

                Fstar=F#0.5*(F+self.last_F)
                Mstar=M#0.5*(M+self.last_M)
                for i in range(nSteps):
                    self.body.setForce((Fstar[0]*opts.free_x[0],
                                        Fstar[1]*opts.free_x[1],
                                        Fstar[2]*opts.free_x[2]))

                    self.body.setTorque((Mstar[0] * opts.free_r[0],
                                         Mstar[1] * opts.free_r[1],
                                         Mstar[2] * opts.free_r[2]))
                    #self.space.collide((self.world,self.contactgroup), near_callback)
                    self.world.step(dt/float(nSteps))
        #self.contactgroup.empty()
        self.last_F[:] = F
        self.last_M[:] = M
        x,y,z = self.body.getPosition()
        u,v,w = self.body.getLinearVel()
        self.barycenters[7,0]=x
        self.barycenters[7,1]=y
        self.barycenters[7,2]=z
        self.last_velocity=copy.deepcopy(self.velocity)
        self.last_position=copy.deepcopy(self.position)
        self.last_rotation=self.rotation.copy()
        self.last_rotation_inv = inv(self.last_rotation)
        self.position=(x,y,z)
        self.velocity=(u,v,w)
        self.rotation=np.array(self.body.getRotation()).reshape(3,3)
        self.h = (self.position[0]-self.last_position[0],
                  self.position[1]-self.last_position[1],
                  self.position[2]-self.last_position[2])
        logEvent("%1.2fsec: pos=(%21.16e, %21.16e, %21.16e) vel=(%21.16e, %21.16e, %21.16e) h=(%21.16e, %21.16e, %21.16e)" % (self.model.stepController.t_model_last,
                                                                                    self.position[0],
                                                                                    self.position[1],
                                                                                    self.position[2],
                                                                                    self.velocity[0],
                                                                                    self.velocity[1],
                                                                                                            self.velocity[2],
                                                                                                            self.h[0],
                                                                                                            self.h[1],
                                                                                                            self.h[2]))


bar = RigidBar(density=100,bar_center=bar_center,bar_dim=opts.bar_dim,barycenters=barycenters,he=he,cfl_target=0.9*opts.cfl,dt_init=opts.dt_init)
