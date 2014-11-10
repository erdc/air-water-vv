from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
   
#  Discretization -- input options  

Refinement = 300 #he=0.005
genMesh=True
useOldPETSc=False
useSuperlu=False#True
timeDiscretization='be'#'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
openTop = True
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

L = (6.0 , 0.7)

#Obstacle or weir dimensions/position
obst_portions = (0.5,0.401) #(width,height)
obst_x_start = 4.5  # start x coordinate of the obstacle; caution to be in the domain's range
obst_x_end = obst_x_start + obst_portions[0] # end x coordinate of the obstacle; caution to be in the domain's range
obst = (obst_x_start,obst_portions[1],obst_x_end) #coordinates of the obstacle to be used to define the boundary

#Background refinement
he = L[0]/float(4*Refinement-1)


# Refinement parameters
x_refine = (2.5 , 4.0 , 5.5) #end of zone 1, end of zone 2, end of zone 3 (zone 4 is up to the right wall)
refinementLevel = (4 , 2) #refinemnt level for zone 1 and zone 2 and 4 respectively zone 3 has the basic refinement level

#Left boundary imposed velocity 
inflow_velocity = 0.047


weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

class PointGauges(AV_base):
    def  __init__(self,gaugeLocations={'pressure_1':(0.5,0.5,0.0)}):
        AV_base.__init__(self)
        self.locations=gaugeLocations
        self.flags={}
        self.files={}#will be opened  later
        pointFlag=100
        for name,point in self.locations.iteritems():
            self.flags[name] = pointFlag
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.p = model.levelModelList[-1].u[0].dof
        self.u = model.levelModelList[-1].u[1].dof
        self.v = model.levelModelList[-1].u[2].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self    
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'.txt','w')
                self.files[name].write('%22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (self.vertices[vnMask,0],self.vertices[vnMask,1],self.p[vnMask],self.u[vnMask],self.v[vnMask]))

class LineGauges(AV_base):
    def  __init__(self,gaugeEndpoints={'pressure_1':((0.5,0.5,0.0),(0.5,1.8,0.0))},linePoints=10):
        import numpy as  np
        AV_base.__init__(self)
        self.endpoints=gaugeEndpoints
        self.flags={}
        self.linepoints={}
        self.files={}#while open later
        pointFlag=1000
        for name,(pStart,pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name]=[]
            for scale in np.linspace(0.0,1.0,linePoints):
                self.linepoints[name].append(p0 + scale*direction)
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.p = model.levelModelList[-1].u[0].dof
        self.u = model.levelModelList[-1].u[1].dof
        self.v = model.levelModelList[-1].u[2].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self    
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'.txt','w')
                for x,y,p,u,v in zip(self.vertices[vnMask,0],self.vertices[vnMask,1],self.p[vnMask],self.u[vnMask],self.v[vnMask]):
                    self.files[name].write('%22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (x,y,p,u,v))

class LineGauges_phi(AV_base):
    def  __init__(self,gaugeEndpoints={'pressure_1':((0.5,0.5,0.0),(0.5,1.8,0.0))},linePoints=10): 
        import numpy as  np
        AV_base.__init__(self)
        self.endpoints=gaugeEndpoints
        self.flags={}
        self.linepoints={}
        self.files={}#while open later
        pointFlag=1000
        for name,(pStart,pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name]=[]
            for scale in np.linspace(0.0,1.0,linePoints):
                self.linepoints[name].append(p0 + scale*direction)
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.phi = model.levelModelList[-1].u[0].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self    
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'_phi.txt','w')
                for x,y,phi in zip(self.vertices[vnMask,0],self.vertices[vnMask,1],self.phi[vnMask]):
                    self.files[name].write('%22.16e %22.16e %22.16e\n' % (x,y,phi))


pointGauges = PointGauges(gaugeLocations={'pointGauge1':(0.05,0.65,0.0)})
lineGauges  = LineGauges(gaugeEndpoints={'lineGauge1':((3.4,0.0,0.0),(3.4,0.1,0.0))},linePoints=2)
lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=2)


structured=False

if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
        domain = Domain.RectangularDomain(L)
    else:
        
        vertices=[
                  [0.0,0.0],#0
                  [obst[0],0.0], #1
                  [obst[0],obst[1]], #2  
                  [obst[2],obst[1]], #3
                  [obst[2],0.0],#4 
                  [L[0],0.0],#5
                  [L[0],L[1]],#6
                  [0.0,L[1]], #7
                  #Mesh refinement region: Points #2, #3, #4 plus the following
                  [x_refine[0],0.0], # 8 
                  [x_refine[0],L[1]], #9
                  [x_refine[1],0.0], #10
                  [x_refine[1],L[1]], #11
                  [x_refine[2],0.0], #12
                  [x_refine[2],L[1]] #13
                  ] 
             #

        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'], 
                     boundaryTags['top'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top']                     
                     ]
                   


        segments=[[0,8], #0
                  [8,10], #1
                  [10,1], #2
                  [1,2], #3
                  [2,3],#4
                  [3,4],#5
                  [4,12],#6
                  [12,5],#7
                  [5,6],#8
                  [6,13],#9
                  [13,11],#10
                  [11,9],#11
                  [9,7],#12
                  [7,0],#13
                  [8,9],#14 
                  [10,11],#15
                  [12,13],#16
                  ]
                  

        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'], 
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['left'],
                      boundaryTags['empty'],
                      boundaryTags['empty'],
                      boundaryTags['empty']
                      ]

 

        regions=[
            [x_refine[1]+0.00001, 0.00001], # Background refinement Zone 3
            [0.00001,0.00001],# Zone 1
            [x_refine[0]+0.00001, 0.00001], # Zone 2
            [x_refine[2]+0.00001, 0.00001] # Zone 4
            ]

        regionFlags=[1,2,3,4]
        trigArea = 0.5*he**2
        refinementArea =(refinementLevel[0]**2 , refinementLevel[1]**2)
        regionConstraints = [trigArea, trigArea*float(refinementArea[0]), trigArea*float(refinementArea[1]), trigArea*float(refinementArea[1])]


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
        triangleOptions="VApq30Dena"

        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))


# Time stepping
T=10.0
dt_fixed = 0.02
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.9
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.5
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

ns_nl_atol_res = max(1.0e-10,0.001*he**2)
vof_nl_atol_res = max(1.0e-10,0.001*he**2)
ls_nl_atol_res = max(1.0e-10,0.001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

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
waterLine_x =obst_x_start+0.1
waterLine_z = 0.462

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[1]-waterLine_z 
    if phi_x <= 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

outflowHeight=-L[1]
def outflowPressure(x,t):
    return (L[1]-x[1])*rho_1*abs(g[1])
    #p_L = L[1]*rho_1*g[1]
    #phi_L = L[1] - outflowHeight
    #phi = x[1] - outflowHeight
    #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
    #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
