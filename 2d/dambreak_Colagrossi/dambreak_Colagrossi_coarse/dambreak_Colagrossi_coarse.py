from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent

#  Discretization -- input options  
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = 16
genMesh=True
movingDomain=False
applyRedistancing=True
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
#L = (0.584,0.350)
L = (3.22 , 1.8)
he = L[0]/float(4*Refinement-1)
#he*=0.5
#he*=0.5
#he*=0.5
#he*=0.5
weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False

class PointGauges(AV_base):
    def  __init__(self, gauges=((('u','v'), ((0.5, 0.5, 0), (1, 0.5, 0))),
                                (('p',),    ((0.5, 0.5, 0),))),
                  activeTime = (0, 0.5),
                  sampleRate = 0,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt'):

        AV_base.__init__(self)
        self.gauges=gauges
        self.measured_fields = set()

        # build up dictionary of location information from gauges
        # dictionary of dictionaries, outer dictionary is keyed by location (3-tuple)
        # inner dictionaries contain monitored fields, and closest node
        # closest_node is None if this process does not own the node
        self.locations = {}
        for gauge in gauges:
            fields, locations = gauge
            self.measured_fields.update(fields)
            for location in locations:
                # initialize new dictionary of information at this location
                if location not in self.locations:
                    l_d = {}
                    l_d['fields'] = set()
                    self.locations[location] = l_d
                # add any currently unmonitored fields
                self.locations[location]['fields'].update(fields)

        self.activeTime = activeTime
        self.sampleRate = sampleRate
        self.fileName = fileName
        self.file = open(self.fileName, 'w')
        self.flags={}
        self.files={}


    def findNearestNode(self, location):
        """Given a gauge location, attempts to locate the most suitable process for monitoring information about
        this location, as well as the node on the process closest to the location.
        """
        from proteus import flcbdfWrappers
        from proteus import Comm
        import numpy as np

        comm = Comm.get()

        node_distances = np.linalg.norm(self.vertices-location, axis=1)
        nearest_node = np.argmin(node_distances)
        nearest_node_distance = node_distances[nearest_node]

        # mpi4py refactor will simplify the below code
        global_min_distance = flcbdfWrappers.globalMin(nearest_node_distance)
        if nearest_node_distance > global_min_distance:
            tie_breaking_nearest_node_distance = nearest_node_distance + comm.size()
        elif nearest_node_distance < global_min_distance:
            raise FloatingPointError('Unexpected result in findNearestNode')
        else: # resolve ties with rank
            tie_breaking_nearest_node_distance = nearest_node_distance + comm.rank()
        tie_breaking_global_min_distance = flcbdfWrappers.globalMin(tie_breaking_nearest_node_distance)
        if tie_breaking_global_min_distance == tie_breaking_nearest_node_distance:
            return nearest_node
        else:
            return None

    def buildQuantityRow(self, m, quantity_id, quantity):
        """ Builds up gauge operator, currently just sets operator to use value of nearest node.
        TODO: get correct element from neighboring node and use element assembly coefficients
        """

        location, node = quantity
        m[quantity_id, node] = 1


    def attachModel(self,model,ar):
        import numpy as np
        from petsc4py import PETSc
        from collections import defaultdict

        self.model=model
        self.fieldNames = model.levelModelList[-1].coefficients.variableNames
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.measured_quantities = defaultdict(list)
        self.m = {}

        for location, l_d in self.locations.iteritems():
            l_d['nearest_node'] = self.findNearestNode(location)
            if l_d['nearest_node'] is not None:
                for field in l_d['fields']:
                    self.measured_quantities[field].append((location, l_d['nearest_node']))

        num_owned_nodes = model.levelModelList[-1].mesh.nNodes_global

        self.m = []
        self.field_ids = []
        self.dofsVecs = []
        self.gaugesVecs = []

        for field in self.measured_fields:
            m = PETSc.Mat().create(PETSc.COMM_SELF)
            m.setSizes([len(self.measured_quantities[field]), num_owned_nodes])
            m.setType('aij')
            m.setUp()
            # matrices are a list in same order as fields
            self.m.append(m)
            field_id = self.fieldNames.index(field)
            self.field_ids.append(field_id)
            # dofs are a list in same order as fields as well
            dofs = self.model.levelModelList[-1].u[field_id].dof
            dofsVec = PETSc.Vec().createWithArray(dofs, comm=PETSc.COMM_SELF)
            self.dofsVecs.append(dofsVec)


        for field, m in zip(self.measured_fields, self.m):
            for quantity_id, quantity in enumerate(self.measured_quantities[field]):
                self.buildQuantityRow(m, quantity_id, quantity)
            gaugesVec = PETSc.Vec().create(comm=PETSc.COMM_SELF)
            gaugesVec.setSizes(len(self.measured_quantities[field]))
            gaugesVec.setUp()
            self.gaugesVecs.append(gaugesVec)

        for m in self.m:
            m.assemble()

        self.outputHeader()
        # this is currently broken for initial time, need to fix initial model time
        # or enforce that calculate is called as soon as possible
        # after model time is set up
        # time = self.get_time()
        time = 0
        self.outputRow(time)
        self.last_output = time

        return self

    def get_time(self):
        """ Returns the current model time"""
        return self.model.levelModelList[-1].timeIntegration.t

    def attachAuxiliaryVariables(self,avDict):
        return self    

    def outputHeader(self):
        """ Outputs a single header for a CSV style file to self.file"""
        self.file.write("time ")
        for field in self.measured_fields:
            for quantity in self.measured_quantities[field]:
                location, node = quantity
                self.file.write(", %s [%22.16e %22.16e %22.16e]" % (field, location[0], location[1], location[2]))
        self.file.write('\n')

    def outputRow(self, time):
        """ Outputs a single row of currently calculated gauge data to self.file"""
        self.file.write("%22.16e " % time)
        for gaugesVec in self.gaugesVecs:
            for quantity in gaugesVec.getArray():
                self.file.write(", %22.16e" % (quantity,))
        self.file.write('\n')
        # disable this for better performance, but risk of data loss on crashes
        self.file.flush()
        self.last_output = time

    def calculate(self):
        """ Computes current gauge values, updates open output files
        """

        time = self.get_time()

        if self.activeTime[0] <= time <= self.activeTime[1] and time >= self.last_output + self.sampleRate:
            for field, m, dofsVec, gaugesVec in zip(self.measured_fields, self.m, self.dofsVecs, self.gaugesVecs):
                m.mult(dofsVec, gaugesVec)
            self.outputRow(time)

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
        self.tt=model.levelModelList[-1].timeIntegration.t
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
                    self.files[name].write('%22.16e %22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (self.tt,x,y,p,u,v))

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
        self.tt= model.levelModelList[-1].timeIntegration.t
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
                    self.files[name].write('%22.16e %22.16e %22.16e %22.16e\n' % (self.tt,x,y,phi))

pointGauges = PointGauges()
lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825':((0.495,0.0,0.0),(0.495,1.8,0.0))},linePoints=20)
#'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)


if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        domain = Domain.RectangularDomain(L)
        nnx=4*Refinement
        nny=2*Refinement
    else:
        vertices=[[0.0,0.0],#0
                  [L[0],0.0],#1
                  [L[0],L[1]],#2
                  [0.0,L[1]]]#3
        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['top']]
        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,0]]
        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['left']]
        regions=[[1.2 ,0.6]]
        regionFlags=[1]

#        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
#            vertices.append(gaugeCoordinates)
#            vertexFlags.append(pointGauges.flags[gaugeName])

        for gaugeName,gaugeLines in lineGauges.linepoints.iteritems():
            for gaugeCoordinates in gaugeLines:
                vertices.append(gaugeCoordinates)
                vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
# Time stepping
T=0.2
dt_fixed = 0.01
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False#True
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
    epsFact_consrv_diffusion = 0.1
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

ns_nl_atol_res = max(1.0e-8,0.001*he**2)
vof_nl_atol_res = max(1.0e-8,0.001*he**2)
ls_nl_atol_res = max(1.0e-8,0.001*he**2)
rd_nl_atol_res = max(1.0e-8,0.005*he)
mcorr_nl_atol_res = max(1.0e-8,0.001*he**2)
kappa_nl_atol_res = max(1.0e-8,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-8,0.001*he**2)

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
waterLine_x = 1.2
waterLine_z = 0.6

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

