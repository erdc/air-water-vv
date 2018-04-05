from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges
from proteus import WaveTools as WT
import vegZoneVelocityInterp as vegZoneInterp


from proteus import Context
opts=Context.Options([
    ("wave_type", 'linear', "type of waves generated: 'linear', 'Nonlinear', 'single-peaked', 'double-peaked'"),
    ("depth", 0.457, "water depth at leading edge of vegetation (not at wave generator)[m]"),
    ("wave_height", 0.192, "wave height at leading edget of vegetation [m]"),
    ("peak_period", 2.0, "Peak period [s]"),
    ("peak_period2", 1.5, "Second peak period (only used in double-peaked case)[s]"),
    ("peak_wavelength",3.91,"Peak wavelength in [m]"),
    ("parallel", False, "Run in parallel"),
    ("gauges", True, "Enable gauges"),
    ("tank_height", 1.0, "height of wave tank")])

#wave generator
windVelocity = (0.0,0.0, 0.0)
veg_platform_height = 17.2/44.0 + 6.1/20.0
depth = opts.depth #veg_platform_height + opts.depth
inflowHeightMean = depth
inflowVelocityMean = (0.0,0.0,0.0)
period = opts.peak_period
omega = 2.0*math.pi/period
waveheight = opts.wave_height
amplitude = waveheight/ 2.0
wavelength = opts.peak_wavelength
k = 2.0*math.pi/wavelength
waveDir = numpy.array([1,0,0])
g = numpy.array([0,-9.81,0])
if opts.wave_type == 'linear':
    waves = WT.MonochromaticWaves(period = period, # Peak period
                                  waveHeight = waveheight, # Height
                                  depth = depth, # Depth
                                  mwl = inflowHeightMean, # Sea water level
                                  waveDir = waveDir, # waveDirection
                                  g = g, # Gravity vector, defines the vertical
                                  waveType="Linear")
elif opts.wave_type == 'Nonlinear':
    waves = WT.MonochromaticWaves(period = period, # Peak period
                                  waveHeight = waveheight, # Height
                                  wavelength = wavelength,
                                  depth = depth, # Depth
                                  mwl = inflowHeightMean, # Sea water level
                                  waveDir = waveDir, # waveDirection
                                  g = g, # Gravity vector, defines the vertical
                                  waveType="Fenton",
                                  Ycoeff = [0.04160592, #Surface elevation Fourier coefficients for non-dimensionalised solution
                                       0.00555874,
                                       0.00065892,
                                       0.00008144,
                                       0.00001078,
                                       0.00000151,
                                       0.00000023,
                                       0.00000007],
                                  Bcoeff = [0.05395079,
                                       0.00357780,
                                       0.00020506,
                                       0.00000719,
                                       -0.00000016,
                                       -0.00000005,
                                       0.00000000,
                                       0.00000000])
elif opts.wave_type == 'single-peaked':
    waves = WT.RandomWaves( Tp = period, # Peak period
                            Hs = waveheight, # Height
                            mwl = inflowHeightMean, # Sea water level
                            depth = depth, # Depth
                            #fp = 1./period, #peak Frequency
                            bandFactor = 2.0, #fmin=fp/Bandfactor, fmax = Bandfactor * fp
                            N = 101, #No of frequencies for signal reconstruction
                            #mwl = inflowHeightMean, # Sea water level
                            waveDir = waveDir, # waveDirection
                            g = g,
                            spectName = "JONSWAP") # Gravity vector, defines the vertical
                            #gamma=3.3,
                            #spec_fun = JONSWAP)
elif opts.wave_type == 'double-peaked':
    waves = WT.DoublePeakedRandomWaves( #Tp = period, # Peak period
                                        Hs = waveheight, # Height
                                        d = depth, # Depth
                                        fp = 1./period, #peak Frequency
                                        bandFactor = 2.0, #fmin=fp/Bandfactor, fmax = Bandfactor * fp
                                        N = 101, #No of frequencies for signal reconstruction
                                        mwl = inflowHeightMean, # Sea water level
                                        waveDir = waveDir, # waveDirection
                                        g = g, # Gravity vector, defines the vertical
                                        gamma=10.0,
                                        spec_fun = WT.JONSWAP,
                                        Tp_2 = opts.peak_period2)

gauges=opts.gauges
#  Discretization -- input options
genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=not opts.parallel
timeDiscretization='be'#'be','vbdf','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = True #False #True #False #True #False
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

#for debugging, make the tank short
#L = (45.4,opts.tank_height)
he = 0.055 #0.025 #float(wavelength)/130.0 #100.0 #50.0#0.0#100
L = (12.2, opts.tank_height, 0.4 ) #0.5)
GenerationZoneLength = wavelength
AbsorptionZoneLength= 45.4-37.9-28.7
spongeLayer = True
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
epsFact_solid = xSponge/2.0
#zone 2
xSponge_2 = 37.9-28.7
xRelaxCenter_2 = 0.5*(37.9+45.4-28.7)
epsFact_solid_2 = AbsorptionZoneLength/2.0

weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False


gauge_dx=0.075
PGL=[]
gauge_top = L[1] #19.5/44.0 + 1.5
veg_gauge_bottom_y = 0.0 #17.2/44.0 + 6.1/20.0
LGL=[#[(6.1, (6.1 - 5.4)/44.0, 0.5*L[2]), (6.1,gauge_top,0.5*L[2])],#1 Goda
     #[#(6.4, (6.4 - 5.4)/44.0, 0.5*L[2]), (6.4,gauge_top,0.5*L[2])],#2 Goda
     #[#(7.0, (7.0 - 5.4)/44.0, 0.5*L[2]), (7.0,gauge_top,0.5*L[2])],#3 Goda
     #[#(26.0, veg_gauge_bottom_y, 0.5*L[2]), (26.0,gauge_top,0.5*L[2])],#4 veg
     [(1.2, veg_gauge_bottom_y, 0.5*L[2]), (1.2,gauge_top,0.5*L[2])],#5
     [(1.7, veg_gauge_bottom_y, 0.5*L[2]), (1.7,gauge_top,0.5*L[2])],#6
     [(2.2, veg_gauge_bottom_y, 0.5*L[2]), (2.2,gauge_top,0.5*L[2])],#7
     [(2.8, veg_gauge_bottom_y, 0.5*L[2]), (2.8,gauge_top,0.5*L[2])],#8
     [(3.8, veg_gauge_bottom_y, 0.5*L[2]), (3.8,gauge_top,0.5*L[2])],#9
     [(5.3, veg_gauge_bottom_y, 0.5*L[2]), (5.3,gauge_top,0.5*L[2])],#10
     [(7.0, veg_gauge_bottom_y, 0.5*L[2]), (7.0,gauge_top,0.5*L[2])],#11
     [(8.7, veg_gauge_bottom_y, 0.5*L[2]), (8.7,gauge_top,0.5*L[2])],#12
     [(10.5, veg_gauge_bottom_y, 0.5*L[2]), (10.5,gauge_top,0.5*L[2])]]#13
# for i in range(0,int(L[0]/gauge_dx+1)): #+1 only if gauge_dx is an exact
#  PGL.append([gauge_dx*i,0.5,0])
#  LGL.append([(gauge_dx*i,0.0,0),(gauge_dx*i,L[1],0)])


gaugeLocations=tuple(map(tuple,PGL))
columnLines=tuple(map(tuple,LGL))


#pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
#                                (('p',),    gaugeLocations)),
#                  activeTime = (0, 1000.0),
#                  sampleRate = 0,
#                  fileName = 'combined_gauge_0_0.5_sample_all.txt')

#print gaugeLocations
#print columnLines

fields = ('vof',)

columnGauge = LineIntegralGauges(gauges=((fields, columnLines),),
                                 fileName='column_gauge3D.csv')

#v_resolution = max(he,0.05)
#linePoints = int((gauge_top - veg_gauge_bottom_y)/v_resolution)
# lineGauges  = LineGauges(gauges=((('u','v'),#fields in gauge set
#                                   (#lines for these fields
#                                       ((26.0, veg_gauge_bottom_y, 0.5*L[2]),(26.0, gauge_top, 0.5*L[2])),
#                                   ),#end  lines
#                               ),#end gauge set
#                              ),#end gauges
#                          fileName="vegZoneVelocity.csv")

#lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)


# if useHex:
#     nnx=ceil(L[0]/he)+1
#     nny=ceil(L[1]/he)+1
#     hex=True
#     domain = Domain.RectangularDomain(L)
# else:
#     boundaries=['left','right','bottom','top','front','back']
#     boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
#     if structured:
#         nnx=ceil(L[0]/he)+1
#         nny=ceil(L[1]/he)+1
#     elif spongeLayer:
#         #tp = opts.tank_height
#         bp = 20.0*(opts.tank_height-0.2)
#         vertices=[[0.0,                                                   0.0                 ],#0 begin wave paddle bottom
#                   [5.4,                                                   0.0                 ],#1 end wave paddle bottom, begin incline 1
#                   [5.4 + 17.2,                                            17.2/44.0           ],#2 end incline 1, begin incline 2
#                   [5.4 + 17.2 + 6.1,                                      17.2/44.0 + 6.1/20.0  ],#3 end incline 2, begin pre veg platform
#                   [5.4 + 17.2 + 6.1 + 1.2,                                17.2/44.0 + 6.1/20.0  ],#4 end pre veg platform, begin veg zone
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8,                          17.2/44.0 + 6.1/20.0  ],#5 end veg zone, begin post veg platform
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2,                    17.2/44.0 + 6.1/20.0  ],#6 -- sponge, end post veg platorm, begin slope bottom
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2 + bp,            17.2/44.0 + 6.1/20.0 + bp/20.0],#7 end slope bottom
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2 + bp,            17.2/44.0 + 6.1/20.0 + bp/20.0+0.2],#8 end slope top
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2,                    19.5/44.0 + 1.5],#9 -- sponge begin sponge top
#                   [0.0,                                                   19.5/44.0 + 1.5]]#10 begin wave paddle top

#         vertexFlags=[boundaryTags['bottom'],#0
#                      boundaryTags['bottom'],#1
#                      boundaryTags['bottom'],#2
#                      boundaryTags['bottom'],#3
#                      boundaryTags['bottom'],#4
#                      boundaryTags['bottom'],#5
#                      boundaryTags['bottom'],#6
#                      boundaryTags['bottom'],#7
#                      boundaryTags['top'],#8
#                      boundaryTags['top'],#9
#                      boundaryTags['top']]#10
#         segments=[[0,1],#0
#                   [1,2],#1
#                   [2,3],#2
#                   [3,4],#3
#                   [4,5],#4
#                   [5,6],#5
#                   [6,7],#6
#                   [7,8],#7
#                   [8,9],#8
#                   [9,10],#8
#                   [10,0],#9
#                   [6,9]]#10

#         segmentFlags=[boundaryTags['bottom'],#0
#                       boundaryTags['bottom'],#1
#                       boundaryTags['bottom'],#2
#                       boundaryTags['bottom'],#3
#                       boundaryTags['bottom'],#4
#                       boundaryTags['bottom'],#5
#                       boundaryTags['bottom'],#6
#                       boundaryTags['right'],#7
#                       boundaryTags['top'],#8
#                       boundaryTags['top'],#9
#                       boundaryTags['left'],#10
#                       0]#11

#         regions=[[0.5,0.5],
#                   [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2+1.0,                    17.2/44.0 + 6.1/20.0 + 1.0]]
#         regionFlags=[1,2]
#         domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
#                                                       vertexFlags=vertexFlags,
#                                                       segments=segments,
#                                                       segmentFlags=segmentFlags,
#                                                       regions=regions,
#                                                       regionFlags=regionFlags)
#         #go ahead and add a boundary tags member
#         domain.boundaryTags = boundaryTags
#         domain.writePoly("mesh")
#         domain.writePLY("mesh")
#         domain.writeAsymptote("mesh")
#         triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
#         print triangleOptions
#         logEvent("""Mesh generated using: triangle -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
#         porosityTypes      = numpy.array([1.0,
#                                           1.0,
#                                           1.0])
#         dragAlphaTypes = numpy.array([0.0,
#                                       0.0,
#                                       0.5/1.004e-6])
#         dragBetaTypes = numpy.array([0.0,0.0,0.0])

#         epsFact_solidTypes = np.array([0.0,0.0,epsFact_solid_2])

#     else:
#         vertices=[[0.0,0.0],#0
#                   [L[0],0.0],#1
#                   [L[0],L[1]],#2
#                   [0.0,L[1]]]#3

#         vertexFlags=[boundaryTags['bottom'],
#                      boundaryTags['bottom'],
#                      boundaryTags['top'],
#                      boundaryTags['top']]
#         segments=[[0,1],
#                   [1,2],
#                   [2,3],
#                   [3,0]
#                   ]
#         segmentFlags=[boundaryTags['bottom'],
#                       boundaryTags['right'],
#                       boundaryTags['top'],
#                       boundaryTags['left']]

#         regions=[ [ 0.1*L[0] , 0.1*L[1] ],
#                   [0.95*L[0] , 0.95*L[1] ] ]
#         regionFlags=[1,2]
#         domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
#                                                       vertexFlags=vertexFlags,
#                                                       segments=segments,
#                                                       segmentFlags=segmentFlags,
#                                                       regions=regions,
#                                                       regionFlags=regionFlags)
#         #go ahead and add a boundary tags member
#         domain.boundaryTags = boundaryTags
#         domain.writePoly("mesh")
#         domain.writePLY("mesh")
#         domain.writeAsymptote("mesh")
#         triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)

#         logEvent("""Mesh generated using: triangle -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
        domain = Domain.RectangularDomain(L)
    elif spongeLayer:
        # vertices=[[0.0,0.0,0.0],#0
        #           [xSponge,0.0,0.0],#1
        #           [xSponge_2,0.0,0.0],#2 
        #           [L[0],0.0,L[2]-0.2],#3
        #           [L[0], 0.0, L[2]],#4
        #           [xSponge_2, 0.0, L[2]],#5
        #           [xSponge,0.0,L[2]],#6
        #           [0.0,0.0,L[2]]]#7
        
               
        # vertexFlags=[boundaryTags['front'],
        #              boundaryTags['front'],
        #              boundaryTags['front'],
        #              boundaryTags['front'],
        #              boundaryTags['front'],
        #              boundaryTags['front'],            
        #              boundaryTags['front'],       
        #              boundaryTags['front']]


        # for v,vf in zip(vertices,vertexFlags):
        #     vertices.append([v[0],L[1], v[2]])
        #     vertexFlags.append(boundaryTags['back'])

        # print vertices
        # print vertexFlags

        # segments=[[0,1],
        #           [1,2],
        #           [2,3],
        #           [3,4],
        #           [4,5],
        #           [5,6],
        #           [6,7],
        #           [7,0],
        #           [1,6],
        #           [2,5]]
                 
        # segmentFlags=[boundaryTags['bottom'],
        #              boundaryTags['bottom'],
        #              boundaryTags['bottom'],                   
        #              boundaryTags['right'],
        #              boundaryTags['top'],
        #              boundaryTags['top'],
        #              boundaryTags['top'],
        #              boundaryTags['left'],
        #              boundaryTags['empty'],
        #              boundaryTags['empty'] ]
        

        bp = 20.0*(opts.tank_height-0.2)
        vertices=[##[0.0,                                                   0.0 ,0.0                ],#0 begin wave paddle bottom
                  #[5.4,                                                   0.0     ,0.0            ],#1 end wave paddle bottom, begin incline 1
                  #[5.4 + 17.2,                                            17.2/44.0   ,0.0        ],#2 end incline 1, begin incline 2
                  [5.4 + 17.2 + 6.1,                                      17.2/44.0 + 6.1/20.0 ,0.0 ],#3 end incline 2, begin pre veg platform
                  [5.4 + 17.2 + 6.1 + 1.2,                                17.2/44.0 + 6.1/20.0 ,0.0 ],#4 end pre veg platform, begin veg zone
                  [5.4 + 17.2 + 6.1 + 1.2 + 9.8,                          17.2/44.0 + 6.1/20.0 ,0.0 ],#5 end veg zone, begin post veg platform
                  [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2,                    17.2/44.0 + 6.1/20.0 ,0.0 ],#6 -- sponge, end post veg platorm, begin slope bottom
                  [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2 + bp,            17.2/44.0 + 6.1/20.0 + bp/20.0 ,0.0],#7 end slope bottom
                  [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2 + bp,            17.2/44.0 + 6.1/20.0 + bp/20.0+0.2 ,0.0],#8 end slope top
                  [5.4 + 17.2 + 6.1 + 1.2 + 9.8 + 1.2,                    19.5/44.0 + 1.5 ,0.0],#9 -- sponge begin sponge top
                  [28.7,                                                   19.5/44.0 + 1.5 ,0.0]]#10 begin wave paddle top
        vertices = np.array(vertices)
        vertices[:,0] -= 28.7
        vertices[:,1] -= 17.2/44.0 + 6.1/20.0
        vertices = vertices.tolist()
        

        vertexFlags=[#boundaryTags['bottom'],#0
                     #boundaryTags['bottom'],#1
                     #boundaryTags['bottom'],#2
                     boundaryTags['bottom'],#3
                     boundaryTags['bottom'],#4
                     boundaryTags['bottom'],#5
                     boundaryTags['bottom'],#6
                     boundaryTags['bottom'],#7
                     boundaryTags['top'],#8
                     boundaryTags['top'],#9
                     boundaryTags['top']]#10

        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1], L[2]])
            #vertexFlags.append(boundaryTags['back'])
        vertexFlags += vertexFlags
        vertices[0][0] = 0.0
        vertices[8][0] = 0.0
        segments=[[0,1],#0 bottom
                  [1,2],#1 bottom
                  [2,3],#2 bottom
                  [3,4],#3 bottom
                  [4,5],#4 bottom
                  [5,6],#5 bottom
                  [6,7],#6 botton
                  [7,8],#7 right
                  [8,9],#8 top 
                  [9,10],#9 top
                  [10,0],#10 left
                  [6,9]]#11 interior
        segments=[[0,1], #0 bottom
                  [1,2], #1 bottom
                  [2,3], #2 bottom
                  [3,4], #3 bottom
                  [4,5], #4 right
                  [5,6], #5 top
                  [6,7], #6 top
                  [7,0], #7 left,
                  [3,6]] #8 interior
                  
                  

        segmentFlags=[# boundaryTags['bottom'],#0
                      # boundaryTags['bottom'],#1
                      # boundaryTags['bottom'],#2
                      boundaryTags['bottom'],#3
                      boundaryTags['bottom'],#4
                      boundaryTags['bottom'],#5
                      boundaryTags['bottom'],#6
                      boundaryTags['right'],#7
                      boundaryTags['top'],#8
                      boundaryTags['top'],#9
                      boundaryTags['left'],#10
                      0]#11

        facets=[]
        facetFlags=[]

        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+8,s[0]+8]])
            facetFlags.append(sF)

        #bf=[[0,1,6,7],[1,2,5,6],[2,3,4,5]] #
        #bf=[[0,1,2,3,4,5,6,9,10], [6,7,8,9]]
        bf=[[0,1,2,3,6,7], [3,4,5,6]]
        #bf=[[0,1,6,7], [1,2
        tf=[]
        for i in range(0,2):
         facets.append([bf[i]])
         facetFlags.append(boundaryTags['front'])
         tf=[ss + 8 for ss in bf[i]]
         facets.append([tf])
         facetFlags.append(boundaryTags['back'])


        # for i in range(0,2):
        #  facetFlags.append(boundaryTags['front'])
        #  facetFlags.append(boundaryTags['back'])
        # bf=[[6,9,6+11,9+11]]
        # for i in range(0,1):
        #  facets.append([bf[i]])
        #  facetFlags.append(boundaryTags['empty'])
        

        print facets
        print facetFlags

        regions=[ [ 0.1*L[0] , 0.1*L[1] , 0.5*L[2]],
                  [0.95*L[0] , 0.95*L[1] , 0.5*L[2]] ]
        # regions=[[xRelaxCenter, 0.5*L[1],0.5*L[2]],
        #          [xRelaxCenter_2, 0.5*L[1], 0.5*L[2]],
        #          [0.5*L[0],0.5*L[1], 0.5*L[2]]]
        regionFlags=[1,2]
        # import pdb
        # pdb.set_trace()
        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags,
                                                     )
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

        porosityTypes      = numpy.array([1.0,
                                          1.0,
                                          1.0# ,
                                          # 1.0
                                      ])
        dragAlphaTypes = numpy.array([0.0,
                                      0.0,
                                      0.5/1.004e-6])

        dragBetaTypes = numpy.array([0.0,0.0,0.0]) #,0.0])
        
        epsFact_solidTypes = np.array([0.0,0.0,0.0]) #,0.0])

    else:             
        vertices=[[0.0,0.0,0.0],#0
                  [L[0],0.0,0.0],#1
                  [L[0],L[1],0.0],#2       
                  [0.0,L[1],0.0]]#3
        
               
        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom']]


        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1],L[2]])
            vertexFlags.append(boundaryTags['top'])

        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,0]]
                 
        segmentFlags=[boundaryTags['front'],                   
                     boundaryTags['right'],
                     boundaryTags['back'],
                     boundaryTags['left']]

        facets=[]
        facetFlags=[]

        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+4,s[0]+4]])
            facetFlags.append(sF)

        bf=[[0,1,2,3]]
        tf=[]
        for i in range(0,1):
         facets.append([bf[i]])
         tf=[ss + 4 for ss in bf[i]]
         facets.append([tf])

        for i in range(0,1):
         facetFlags.append(boundaryTags['bottom'])
         facetFlags.append(boundaryTags['top'])

        for s,sF in zip(segments,segmentFlags):
            segments.append([s[1]+4,s[0]+4])
            segmentFlags.append(sF)
        

        regions=[[0.5*L[0],0.5*L[1], 0.0]]
        regionFlags=[1]

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
T=200.0 #480.0  #480.0 #40*period
dt_fixed = period/30.0#2.0*0.5/20.0#T/2.0#period/21.0
dt_init = min(0.001*dt_fixed,0.001)
runCFL=5.90
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.0
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
    epsFact_redistance = 1.5
    epsFact_consrv_diffusion = 10.0
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
rd_nl_atol_res = max(1.0e-10,0.05*he)
mcorr_nl_atol_res = max(1.0e-10,0.001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
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
waterLine_x = 4*L[0]
waterLine_z = inflowHeightMean


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


def theta(x,t):
    return k*x[0] - omega*t + pi/2.0

def z(x):
    return x[1] - inflowHeightMean

#sigma = omega - k*inflowVelocityMean[0]
h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0

residence_time = 2.0*period
def ramp(t):
    if t < residence_time:
        return 1.0-exp(-5.0*t/residence_time)
    else:
        return 1.0
#waveData
def waveHeight(x,t):
    # import pdb
    # pdb.set_trace()
    return  float(vegZoneInterp.interp_phi.__call__(t%vegZoneInterp.time[-1])) 
    # return inflowHeightMean + ramp(t)*waves.eta(x,t)
# def waveVelocity_u(x,t):
#     return waves.u(x,t)[0]*ramp(t)
# def waveVelocity_v(x,t):
#     return waves.u(x,t)[1]*ramp(t)
# def waveVelocity_w(x,t):
#     return waves.u(x,t)[2]*ramp(t)

def waveHeight(x,t):
    return  float(vegZoneInterp.interp_phi.__call__(t%vegZoneInterp.time[-1])) 
def waveVelocity_u(x,t):
    return vegZoneInterp.interpU.__call__(t%vegZoneInterp.time[-1],x[1])[0][0] 
def waveVelocity_w(x,t):
    return 0.0
def waveVelocity_v(x,t):
    return vegZoneInterp.interpW.__call__(t%vegZoneInterp.time[-1],x[1])[0][0]

# import pdb
# pdb.set_trace()
#solution variables

def wavePhi(x,t):
    return x[1] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = waveVelocity_v(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowVelocity_w(x,t):
    waterspeed = waveVelocity_w(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[2]+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

outflowHeight=inflowHeightMean

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)

def outflowPhi(x,t):
    return x[1] - outflowHeight

def outflowPressure(x,t):
  if x[1]>inflowHeightMean:
    return (L[1]-x[1])*rho_1*abs(g[1])
  else:
    return (L[1]-inflowHeightMean)*rho_1*abs(g[1])+(inflowHeightMean-x[1])*rho_0*abs(g[1])


    #p_L = L[1]*rho_1*g[1]
    #phi_L = L[1] - outflowHeight
    #phi = x[1] - outflowHeight
    #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
    #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

# def twpflowVelocity_w(x,t):
#     return 0.0

def zeroVel(x,t):
    return 0.0

from collections import  namedtuple


def zeroVel(x,t):
    return 0.0

from collections import  namedtuple

RelaxationZone = namedtuple("RelaxationZone","center_x sign u v w")

class RelaxationZoneWaveGenerator(AV_base):
    """ Prescribe a velocity penalty scaling in a material zone via a Darcy-Forchheimer penalty

    :param zones: A dictionary mapping integer material types to Zones, where a Zone is a named tuple
    specifying the x coordinate of the zone center and the velocity components
    """
    def __init__(self,zones):
        assert isinstance(zones,dict)
        self.zones = zones
    def calculate(self):
        for l,m in enumerate(self.model.levelModelList):
            for eN in range(m.coefficients.q_phi.shape[0]):
                mType = m.mesh.elementMaterialTypes[eN]
                if self.zones.has_key(mType):
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN,k]
                        m.coefficients.q_phi_solid[eN,k] = self.zones[mType].sign*(self.zones[mType].center_x - x[0])
                        m.coefficients.q_velocity_solid[eN,k,0] = self.zones[mType].u(x,t)
                        m.coefficients.q_velocity_solid[eN,k,1] = self.zones[mType].v(x,t)
                        #m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
        m.q['phi_solid'] = m.coefficients.q_phi_solid
        m.q['velocity_solid'] = m.coefficients.q_velocity_solid

rzWaveGenerator = RelaxationZoneWaveGenerator(zones={
                                                    # 1:RelaxationZone(xRelaxCenter,
                                                    #                  1.0,
                                                    #                  twpflowVelocity_u,
                                                    #                  twpflowVelocity_v,
                                                    #                  twpflowVelocity_w),
                                                    2:RelaxationZone(xRelaxCenter_2,
                                                                     -1.0, #currently Hs=1-exp_function
                                                                     zeroVel,
                                                                     zeroVel,
                                                                     zeroVel)})

#beam info
beam_quadOrder=3
beam_useSparse=False
beamFilename="wavetankBeams"
beamLocation=[]
beamLength=[]
beamRadius=[]
EI=[]
GJ=[]
lam = 0.05 #0.07 #0.05 #3.0*2.54/100.0 #57.4e-3
lamx = 3.0**0.5*lam
xs = 1.2
ys = 0.0
xList=[]
yList = []
while xs <= 11.0:
    xList.append(xs)
    xs += lam
while ys<= L[2]:
    yList.append(ys)
    ys+=lamx
for i in xList:
    for j in yList:
        beamLocation.append((i,j))
        beamLength.append(0.415)
        beamRadius.append(0.0032) #32)
        EI.append(0.0142) # needs to be fixed
        GJ.append(2.0*0.0142) # needs to be fixed

xs = 1.2+0.5*lam
ys = 0.5*lamx
xList=[]
yList = []
while xs <= 11.0:
    xList.append(xs)
    xs += lam

while ys<= L[2]:
    yList.append(ys)
    ys+=lamx

for i in xList:
    for j in yList:
        beamLocation.append((i,j))
        beamLength.append(0.415)
        beamRadius.append(0.0032) #32)
        EI.append(0.0142) # needs to be fixed
        GJ.append(2.0*0.0142) # needs to be fixed
nBeamElements = int(beamLength[0]/he*0.5)
nBeamElements=max(nBeamElements,3)
print nBeamElements
