from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
   
#  Discretization -- input options  

Refinement = 8
genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=False
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
L = (1.6 ,0.61, 0.75)
obst_portions = (0.12,0.12,0.75) #(width_x,width_y,height)
obst_x_start = 0.9 # start x coordinate of the obstacle; caution to be in the domain's range
obst_x_end = obst_x_start + obst_portions[0] # end x coordinate of the obstacle; caution to be in the domain's range
obst_y_start =0.25  # start y coordinate of the obstacle; caution to be in the domain's range
obst_y_end=obst_y_start+obst_portions[1] # end y coordinate of the obstacle; caution to be in the domain's range
obst = (obst_x_start,obst_x_end,obst_y_start,obst_y_end,obst_portions[2]) #coordinates of the obstacle to be used to define the boundary


he = L[0]/float(4*Refinement-1)
nLevels = 1
weak_bc_penalty_constant = 100.0
quasi2D=False
if quasi2D:#make tank one element wide
    L = (L[0],he,L[2])
from proteus.MeshAdaptPUMI  import MeshAdaptPUMI 
hmin = he
hmax = 10.0*he
adaptMesh = True
adaptMesh_nSteps = 5
adaptMesh_numIter = 2
MeshAdaptMesh=MeshAdaptPUMI.MeshAdaptPUMI(hmax=hmax, hmin=hmin, numIter=adaptMesh_numIter,sfConfig="isotropic",maType="isotropic")
useModel=False
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

structured=False  


if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back','obst']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
        domain = Domain.RectangularDomain(L)
    else:
                  #plane 1
        vertices=[[0.0,0.0,0.0], #0
                  [obst[0],0.0,0.0], #1
                  [obst[0],0.0,L[2]], #2
                  [0.0,0.0,L[2]], #3
                  #plane 2
                  [0.0,obst[2],0.0], #4
                  [obst[0],obst[2],0.0], #5
                  [obst[0],obst[2],L[2]], #6
                  [0.0,obst[2],L[2]], #7
                  #plane 3                  
                  [0.0,obst[3],0.0], #8
                  [obst[0],obst[3],0.0], #9
                  [obst[0],obst[3],L[2]], #10
                  [0.0,obst[3],L[2]], #11
                  #plane 4
                  [0.0,L[1],0.0], #12
                  [obst[0],L[1],0.0], #13
                  [obst[0],L[1],L[2]], #14
                  [0.0,L[1],L[2]], #15
                  #plane 1 extension 1
                  [obst[1],0.0,0.0], #16
                  [obst[1],0.0,L[2]], #17
                  #plane 2 extension 1
                  [obst[1],obst[2],0.0], #18
                  [obst[1],obst[2],L[2]], #19
                  #plane 3 extension 1
                  [obst[1],obst[3],0.0], #20
                  [obst[1],obst[3],L[2]], #21
                  #plane 4 extension 1                 
                  [obst[1],L[1],0.0], #22
                  [obst[1],L[1],L[2]], #23
                  #plane 1 extension 2
                  [L[0],0.0,0.0], #24
                  [L[0],0.0,L[2]], #25
                  #plane 2 extension 2
                  [L[0],obst[2],0.0], #26
                  [L[0],obst[2],L[2]], #27
                  #plane 3 extension 2
                  [L[0],obst[3],0.0], #28
                  [L[0],obst[3],L[2]], #29
                  #plane 4 extension 2
                  [L[0],L[1],0.0], #30
                  [L[0],L[1],L[2]] #31
                  ]                

        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'], 
                     boundaryTags['top'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'], 
                     boundaryTags['top'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['bottom'],                 
                     boundaryTags['top'],                
                     boundaryTags['bottom'],                 
                     boundaryTags['top'],                
                     boundaryTags['bottom'],                   
                     boundaryTags['top'],                   
                     boundaryTags['bottom'],                 
                     boundaryTags['top'],                  
                     boundaryTags['bottom'],                   
                     boundaryTags['top']                 
                     ]
     #    for v,vF in zip(vertices,vertexFlags):
             
     #        vertices.append([v[0],L[1],v[2]])
     #        vertexFlags.append(vF)
     

        segments=[[0,1],
                  [1,16],
                  [16,24],
                  [24,25],
                  [25,17],
                  [17,2],
                  [2,3],
                  [3,0],
                  [12,13],
                  [13,22],
                  [22,30],
                  [30,31],
                  [31,23],
                  [23,14],
                  [14,15],
                  [15,12],
                  [0,4],
                  [4,8],
                  [8,12],
                  [24,26],
                  [26,28],
                  [28,30],
                  [25,27],
                  [27,29],
                  [29,31],
                  [3,7],
                  [7,11],
                  [11,15],
                  [5,18],
                  [18,19],
                  [19,6],
                  [6,5],
                  [9,20],
                  [20,21],
                  [21,10],
                  [10,9],
                  [5,9],
                  [18,20],
                  [19,21],
                  [6,10]]
  


        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['left'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'], 
                      boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['left'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['top'],
                      boundaryTags['bottom'],
                      boundaryTags['obst'],
                      boundaryTags['top'],
                      boundaryTags['obst'],
                      boundaryTags['bottom'],
                      boundaryTags['obst'],
                      boundaryTags['top'],
                      boundaryTags['obst'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['top'],
                      boundaryTags['top']]


        facets=[[[0,4,8,12,15,11,7,3]],#left
                [[24,26,28,30,31,29,27,25]],#right
                [[5,18,19,6]],  #front facet of obstacle
                [[9,20,21,10]],#back facet of obstacle
                [[5,9,10,6]], #left facet of obstacle
                [[18,20,21,19]], #right facet of obstacle
                [[0,1,5,4]],#bottom
                [[4,5,9,8]],#bottom
                [[8,9,13,12]],#bottom
                [[1,16,18,5]],#bottom
                [[18,16,24,26]],#bottom
                [[18,26,28,20]],#bottom
                [[20,28,30,22]],#bottom
                [[20,22,13,9]],#bottom
                [[3,2,6,7]],  #top
                [[7,6,10,11]],#top
                [[11,10,14,15]], #top
                [[2,17,19,6]],#top
                [[19,17,25,27]],#top
                [[19,27,29,21]],#top
                [[21,29,31,23]],#top
                [[10,21,23,14]],#top
                [[0,1,16,24,25,17,2,3]], #front
                [[12,13,22,30,31,23,14,15]]    #back        
                ]
        

        facetFlags=[boundaryTags['left'],
                    boundaryTags['right'],
                    boundaryTags['obst'],
                    boundaryTags['obst'],
                    boundaryTags['obst'],
                    boundaryTags['obst'],
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['bottom'], 
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['top'],
                    boundaryTags['front'],
                    boundaryTags['back']
                    ]

 
    
    #    for s,sF in zip(segments,segmentFlags):
    #     facets.append([[s[0],s[1],s[1]+4,s[0]+4]])
    #     facetFlags.append(sF)
    #     if s[0]<4:
    #         frontFacet.append(s[0])
    #         backFacet.append(s[0]+4)
    #    facets.append([frontFacet])
    #    facetFlags.append(boundaryTags['front'])
    #    facets.append([backFacet])
    #    facetFlags.append(boundaryTags['back'])
        regions=[]
        regionFlags=[]
        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
# Time stepping
T=0.8
dt_fixed = 0.05
dt_init = min(0.1*dt_fixed,0.1*he)
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
g = [0.0,0.0,-9.8]

# Initial condition
waterLine_x =0.4
waterLine_z = 0.30

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[2]-waterLine_z 
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

