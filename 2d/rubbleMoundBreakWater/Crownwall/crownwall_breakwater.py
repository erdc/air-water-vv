from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
import numpy as np
from proteus.mprans import BodyDynamics as bd
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import Gauges as ga



# Options (Geometry, Physical Properties, Waves, Numerical Options, Obstacle Dimentions)

opts=Context.Options([
    
# Geometry

("tank_height",1.5,"Vertical Dimention of the tank"),
("waterline_x", 10000, "used in levelset"),

# Physical Properties
 ("rho_0", 998.2, "Water density"),
 ("nu_0", 1.004e-6,"Water viscosity"),
 ("rho_1", 1.205, "Air density"),
 ("nu_1",1.5e-5, "Air viscosity"),
 ("sigma_01", 0.,"surface tension"),
 ("g", np.array([0., -9.805, 0.]), "gravity"),
 ("n_mound", 0.39, "porosity of rock mound"),
 ("n_armour", 0.5, "porosity of the armour layer"),



 # alpha  Darcy - porous media

("alpha_armour",0.1 ,"Darcy Coef Ergun 1952"),
("alpha_filter",0.9 ,"Darcy Coef Ergun 1952"),
("alpha_core",0.63 ,"Darcy Coef by Polidoro et al 2018 "),
("alpha_baselayer",0.62 ,"Darcy Coef by Polidoro et al 2018"),

 # beta  Forchheimer - porous media

("beta_armour",105.8 ,"Forchheimer Coef Ergun 1952 "),
("beta_filter", 705.3, "Forchheimer Coef Ergun 1952 "),
("beta_core", 59.63, "Forchheimer Coef  by Polidoro et al 2018 "),
("beta_baselayer",59.49, "Forchheimer Coef  by Polidoro et al 2018"),


# Waves
("Tstart", 0, "Start time"),
("Ntotalwaves", 200.,"totalnumber of waves"),
("fract", 4000., "fraction of duration (Tend = Ntotalwaves*Tp/1.1./fract"),
("x0", np.array([0.,0.,0.]), "Position vector for the tinme series"),
("Tp", 3.5, "Peak wave period"),
("Hs", 0.096, "Significant wave height"),
("mwl", 0.4, "Mean Water level"),
("depth", 0.3 , "Water depth"),
("waveDir", np.array([1.,0.,0.]),"Direction of the waves"),
("N", 2000, "Number of frequency components"),
("bandFactor", 2.0 ,"Spectal Band Factor"),
("spectName", "JONSWAP","Name of Spectral Distribution"),
("spectral_params",{"gamma": 3.3, "TMA":False,"depth": 0.4} ,"Spectral Distribution Characteristics"),
("seed", 420,"Seed for random phases"),
("Lgen",None , "Length of the generation zone"),
("Nwaves", 15, "Number of waves per window"),
("Nfreq", 32 , "Number of fourier components per window"),
("dx", 0.01, "vertical spacing in probes [m] "),
    



# Numerical Options
("refinement_level", 100.,"he=wavelength/refinement_level"),
("cfl", 0.5,"Target cfl"),
("ecH", 1.5,"Smoothing Coefficient"),
("Np", 200.," Output points per period Tp/Np" ),
("dt_init", 0.0001 , "initial time step" )
])
    
    



# --- DOMAIN
domain = Domain.PlanarStraightLineGraphDomain()




################	Random Waves Class	 ################

np.random.seed(opts.seed)
phi = 2*np.pi*np.random.rand(opts.N)
Tend=opts.Ntotalwaves*opts.Tp/1.1

wave = wt.RandomWavesFast(Tstart=opts.Tstart,
                         Tend=Tend,
                         x0=opts.x0,
                         Tp=opts.Tp,
                         Hs=opts.Hs,
                         mwl=opts.mwl,
                         depth=opts.depth,
                         waveDir=opts.waveDir,
                         g=opts.g,
                         N=opts.N,
                         bandFactor=opts.bandFactor,
                         spectName=opts.spectName,
                         spectral_params=opts.spectral_params,
                         phi=phi,
                         Lgen=opts.Lgen,
                         Nwaves=opts.Nwaves,
                         Nfreq=opts.Nfreq,
                         checkAcc=True,
                         fast=True)


# Wave length calculation
wave_length=wave.wavelength           



################	Boundary Declaration	################




boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'porousLayer': None,
                        'moving_porousLayer': None,
                       }

boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'sponge' : 5,

		}


################	TANK OUTLINE GEOMETRY	################ 


                     
vertices=[[0.0,0.0], #0
        [0.39, 0.],   #1
	[3.57, 0.159 ],#2
	[5.42, 0.198], #3
	[6.7, 0.224],   #4
	[7.97, 0.251], #5
	[9.25, 0.278], #6
	[10.52, 0.304], #7
	[11.80, 0.331], #8
	[13.07, 0.357], #9
	[14.17, 0.36], #10
	#[14.74, 0.351],
	#[15.31, 0.362],
	[15.89, 0.36],#11
	[15.94, 0.414],#12
	#[18.54, 0.414],
	#[20.54, 0.414], 
	[15.94+wave_length, 0.414], #13
	[15.94+wave_length, 1.5], #14
	[15.94, 1.5], #15
	[0.,1.5], #16
	[-wave_length,1.5], #17
	[-wave_length, 0 ] #18
]

              


vertexFlags=np.array([1, #0 
                        1, #1 
                        1, #2
                        1, #3
                        1, #4 
                        1, #6 
                        1, #7 
                        1, #8 
                        1, #9
                        1, #10 
                        1,#11
			1,#12
			1,#13
			3,#14
			3,#15
			3,#16
			4,#17
			4,#18
			
                        ])           

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,8],
          [8,9],
          [9,10],
          [10,11],
          [11,12],
	        [12,13],
	        [13,14],
	        [14,15],
	        [15,16],
	        [16,17],
          [17,18],
          [18,0]
      ]

segmentFlags=np.array([ 
	           1,#[0,1]
          1,#[1,2]
          1,#[2,3]
          1,#[3,4]
          1,#[4,5],
          1,#[5,6],
          1,#[6,7],
          1,#[7,8],
          1,#[8,9],
          1,#[9,10],
          1,#[10,11],
          1,#[11,12],
	        1,#[12,13],
	        2,#[13,14],
	        3,#[14,15],
	        3,#[15,16],
	        3,#[16,17],
	        4,#[17,18],
	        1,#[18,0],
 	        5,#[0,16],
	        5,#[12,15],



])



regions=[
[5,1.],
[-0.5*wave_length,0.3],
[15.94+0.5*wave_length,0.4]
] 
        
regionFlags =np.array([1,2,3])  



################# Impermeable  ##############################
                    
obs_boundaryOrientations = {'obstacle': None}

obs_boundaryTags = {'obstacle' : 1,}

obs_vertices=[

            [15.39,0.454],#0
            [15.67,0.454],#1
            [15.67,0.674],#2
            [15.43,0.674],#3
            [15.41,0.835],#4
            [15.39,0.835],#5
            
	    
            ] 

obs_vertexFlags=np.array([1, #0 
                        1, #1
                        1, #2
                        1, #3
                        1, #4
                        1 #5
                         
                        ])           

obs_segments=[[0,1],
              [1,2],
              [2,3],
              [3,4],
              [4,5],
              [5,0]
	      
              ]

obs_segmentFlags=np.array([ 1, #0 
                            1, #1
                            1, #2
                            1, #3
                            1, #4
                            1 #5
                            
                            ])

obs_regions=[[15.4, 0.5]]         
obs_regionFlags =np.array([1])        

obstacle = st.CustomShape(domain, vertices=obs_vertices, vertexFlags=obs_vertexFlags,
                      segments=obs_segments, segmentFlags=obs_segmentFlags,
                      regions=obs_regions, regionFlags=obs_regionFlags,
                      boundaryTags=obs_boundaryTags, 
                      boundaryOrientations=obs_boundaryOrientations)


obstacle.setHoles([[15.4, 0.5]])




##############################################################
#################        porous zones        #################


############################################################
####### layer_1 acropodes


layer1_boundaryOrientations = {'layer_1': None, 'solid_wall':None}





layer1_boundaryTags = {'layer_1' : 1,'solid_wall':2}

layer1_vertices=np.array([(14.459,0.36),#0
			(14.579, 0.36),#1
			(15.124, 0.763),#2
			(15.39, 0.763),#3
			(15.39, 0.835), #4
           		(15.101, 0.835),#5
               
               

])


layer1_vertexFlags=np.array([
2, #0 
2, #1
1, #2
2, #3
2, #4
1 #5
])           



layer1_segments=[[0,1],
[1,2],
[2,3],
[3,4],
[4,5],
[5,0]
]

layer1_segmentFlags=np.array([ 2, #0,1 
                            1, #1,2
                            1, #2,3
                            2, #3,4
                            1, #4,5
                            1 #5,0
                            
                            ])




layer1_regions=[[14.5, 0.37]]         
layer1_regionFlags =np.array([1])        

layer_1 = st.CustomShape(domain,
 vertices=layer1_vertices,
 vertexFlags=layer1_vertexFlags,
 segments=layer1_segments,
 segmentFlags=layer1_segmentFlags,
 regions=layer1_regions,
 regionFlags=layer1_regionFlags,
 boundaryTags=layer1_boundaryTags, 
 boundaryOrientations=layer1_boundaryOrientations)


#########################################
####### layer_2 filter zone


layer2_boundaryOrientations = {'layer_2': None, 'solid_wall':None}

layer2_boundaryTags = {'layer_2' : 1,'solid_wall':2}



layer2_vertices= np.array([(14.579, 0.36),#0
		   (14.669,0.36),
 		   (15.144, 0.713),
                   (15.39, 0.713),
		   (15.39, 0.763),
                   (15.124, 0.763)#5
                   ])
		   

layer2_vertexFlags=np.array([
2, #0 
2, #1
1, #2
2, #3
2, #4
1,#5
])           

layer2_segments=[[0,1], #0
[1,2],#1
[2,3],#2
[3,4],#3
[4,5],#4
[5,0]#5
]

layer2_segmentFlags=np.array([ 2, #0 ,1
                            1, #1,2
                            1, #2,3
                            2, #3,4
                            1, #4,5
                            1 #5,6
                            
                            ])

layer2_regions=[[14.6, 0.37]]         
layer2_regionFlags =np.array([1])        

layer_2 = st.CustomShape(domain,
 vertices=layer2_vertices,
 vertexFlags=layer2_vertexFlags,
 segments=layer2_segments,
 segmentFlags=layer2_segmentFlags,
 regions=layer2_regions,
 regionFlags=layer2_regionFlags,
 boundaryTags=layer2_boundaryTags, 
 boundaryOrientations=layer2_boundaryOrientations)

###########################################
####### layer_3 core


layer3_boundaryOrientations = {'layer_3': None,'solid_wall':None}

layer3_boundaryTags = {'layer_3' : 1,'solid_wall':2}


layer3_vertices=np.array([  (14.669,0.36),
			(15.39, 0.36),
			(15.39, 0.454), #extra point when the interface changes
			(15.39, 0.713),
      (15.144, 0.713),

]) 



 

layer3_vertexFlags=np.array([
2, #0 
2, #1
1, #2 new point when interface changes
1, #3
1, #4


])           

layer3_segments=[[0,1],#0 
[1,2],#1 
[2,3],#2 
[3,4], #3 
[4,0], #4 

]

layer3_segmentFlags=np.array([ 2, #0 0,1
                            1, #1    1,2
                            1, #2    2,3 CHANGED
                            1, #3    3,4
                            1, #4    4,0
                            ])



layer3_regions=[[15., 0.4]]         

layer3_regionFlags =np.array([1])        



layer_3 = st.CustomShape(domain,
 vertices=layer3_vertices,
 vertexFlags=layer3_vertexFlags,
 segments=layer3_segments,
 segmentFlags=layer3_segmentFlags,
 regions=layer3_regions,
 regionFlags=layer3_regionFlags,
 boundaryTags=layer3_boundaryTags, 
 boundaryOrientations=layer3_boundaryOrientations)


#########################################
###### base layer


base_layer_boundaryOrientations = {'base_layer': None,'solid_wall':None}


base_layer_boundaryTags = {'base_layer' : 1, 'solid_wall':2}

base_layer_vertices=np.array([
            (15.39,0.454),
            (15.39,0.36),
            (15.89,0.36),
            (15.94,0.414),
            (15.67,0.454),
])


base_layer_vertexFlags=np.array([
1, #0 
2, #1
2, #2
1, #3
1, #4

])           

base_layer_segments=[[0,1],
[1,2],
[2,3],
[3,4],
[4,0]
]

base_layer_segmentFlags=np.array([ 1, #0 0,1
                            2, #1 1,2
                            1, #2  2,3
                            1, #3  3,4
                            1, #4   4,0
                            
                            ])

base_layer_regions=[[15.41, 0.4]]         
base_layer_regionFlags =np.array([1])        

base_layer = st.CustomShape(domain,
 vertices=base_layer_vertices,
 vertexFlags=base_layer_vertexFlags,
 segments=base_layer_segments,
 segmentFlags=base_layer_segmentFlags,
 regions=base_layer_regions,
 regionFlags=base_layer_regionFlags,
 boundaryTags=base_layer_boundaryTags, 
 boundaryOrientations=base_layer_boundaryOrientations)

#################################################
##### setting the porous media && drag forces ##


n_mound=opts.n_mound
n_armour=opts.n_armour


dragAlpha_accropodes=float(opts.alpha_armour*n_armour**2/opts.nu_0)
dragBeta_accropodes=float(opts.beta_armour*n_armour**3/opts.nu_0)


dragAlpha_filter=float(opts.alpha_filter*n_mound**2/opts.nu_0)
dragBeta_filter=float(opts.beta_filter*n_mound**3/opts.nu_0)


dragAlpha_core=float(opts.alpha_core*n_mound**2/opts.nu_0)
dragBeta_core=float(opts.beta_core*n_mound**3/opts.nu_0)

dragAlpha_baselayer=float(opts.alpha_baselayer*n_mound**2/opts.nu_0)
dragBeta_baselayer=float(opts.beta_baselayer*n_mound**3/opts.nu_0)


########################################
########### porous zones ###############


layer_1.setPorousZones(flags=1, dragAlpha=dragAlpha_accropodes, dragBeta=dragBeta_accropodes, porosity=np.array([n_armour]))

layer_2.setPorousZones(flags=1, dragAlpha=dragAlpha_filter, dragBeta=dragBeta_filter, porosity=np.array([n_armour]))

layer_3.setPorousZones(flags=1, dragAlpha=dragAlpha_core, dragBeta=dragBeta_core, porosity=np.array([n_mound]))


base_layer.setPorousZones(flags=1, dragAlpha=dragAlpha_baselayer,dragBeta=dragBeta_baselayer, porosity=np.array([n_mound]))



################	Mesh Refinement		################

#Characteristic Cell size
he=wave_length/opts.refinement_level
ecH=opts.ecH
smoothing=ecH*he




#######################################################
################	Tank Setup	################      


tank = st.CustomShape(domain,
 			vertices=vertices,
 			vertexFlags=vertexFlags,
			segments=segments,
 			segmentFlags=segmentFlags,
   			regions=regions, 
			regionFlags=regionFlags,
                        boundaryTags=boundaryTags,
 			boundaryOrientations=boundaryOrientations)



tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

tank.BC['sponge'].setNonMaterial()
obstacle.BC['obstacle'].setFreeSlip()




##########################################
##########################################


## Setting internal boundary conditions## 




##########################################
## outer layer accropodes
##########################################

layer_1.BC['layer_1'].setNonMaterial()

layer_1.BC['solid_wall'].setFreeSlip()



##########################################
## filter layer
##########################################

layer_2.BC['layer_2'].setNonMaterial()
#layer_2.BC['solid_wall'].setNoSlip()
layer_2.BC['solid_wall'].setFreeSlip()



##########################################
## core of the breakwater
##########################################


layer_3.BC['layer_3'].setNonMaterial()

layer_3.BC['solid_wall'].setFreeSlip()



##########################################
## base layer below impermeale structure
##########################################


base_layer.BC['base_layer'].setNonMaterial()

base_layer.BC['solid_wall'].setFreeSlip()




################################################################################
################	Generation and Absorption zones		################





dragAlpha = 5*(2*np.pi/opts.Tp)/1e-6


tank.setGenerationZones(flags=2,
			 epsFact_solid=wave_length/2.,
			 center=(-wave_length/2,0.35),
			 orientation=(1.,0.,0.), 
			 waves=wave,
			 dragAlpha=dragAlpha)

tank.setAbsorptionZones(flags=3,
			 epsFact_solid=wave_length/2.,
			 center=(15.94+wave_length*0.5,0.957),
 		         orientation=(-1.,0.,0.),
                         dragAlpha=dragAlpha)




waterLine_x=opts.waterline_x
waterLine_z = opts.mwl



def signedDistance(x):
    phi_x = x[0]- waterLine_x
    phi_y = x[1] - opts.mwl
    if phi_x < 0.0:
        if phi_y < 0.0:
            return max(phi_x, phi_y)
        else:
            return phi_y
    else:
        if phi_y < 0.0:
            return phi_x
        else:
            return np.sqrt(phi_x ** 2 + phi_y ** 2)


class P_IC:
    def __init__(self):
        self.waterLevel=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_height - opts.mwl)*opts.rho_1*opts.g[1] - (opts.mwl - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_height - opts.mwl)*opts.rho_1*opts.g[1]
class AtRest:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest()}
class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()


#################################################################
#################	Two Phase Flow Set up	 ################

Duration= Tend/opts.fract
dt_output = opts.Tp/opts.Np

outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                       dt_init=opts.dt_init,
                                       # cfl=cfl,
                                       dt_output=dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None, # set with SpatialTools,
                                             )


################	Numerical Parameters	################

params = myTpFlowProblem.Parameters

myTpFlowProblem.useSuperLu=False#True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma_01

m = params.Models

#m.rans2p.n.conservativeFlux = {0:'pwl-bdm-opt'}
#m.rans2p.p.coefficients.useVF=0.
#m.rans2p.p.coefficients.weak_bc_penalty_constant=10.
#myTpFlowProblem.movingDomain = True


#m.rans2p.n.maxNonlinearIts=100
#m.rdls.p.coefficients.epsFact=0.75
m.moveMeshElastic.index = 0
m.rans2p.index = 1
m.vof.index = 2
m.ncls.index = 3
m.rdls.index = 4
m.mcorr.index = 5


myTpFlowProblem.movingDomain = True
myTpFlowProblem.useSuperLu=False#True



################	Gauges		################ 


#free surface elevation at the crest

LG=[((15.4,0.836,0.),
     (15.4,1.5,0.))]

x1=15.4


probes=np.linspace(0.836,opts.tank_height,opts.dx)

for i in probes:
	LG.append((x1,i,0.),)



myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables+= [ga.LineGauges(gauges=((('u',), LG),),
					activeTime=(0.,Duration),
					sampleRate=0.,
					fileName='velocity.csv')
					]


myTpFlowProblem.Parameters.Models.vof.auxiliaryVariables+= [ga.LineGauges(gauges=(((('vof'),), LG),),
			activeTime=(0.,Duration),
			sampleRate=0.,
			fileName='vof.csv')
			]


# probe array

LG2=[((0.0,0.,0.),(0.0,0.9,0.))]
probes1=np.linspace(0.,0.9,opts.dx)

for i in probes1:
	LG2.append((0.,i,0.),)

   	
myTpFlowProblem.Parameters.Models.vof.auxiliaryVariables+=[ga.LineIntegralGauges(gauges=(((('vof'),), LG2),),
			activeTime=(0.,Duration),
			sampleRate=0.,
			fileName='line_integral_gauge_generation.csv')]


#Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']





##################################################################################





