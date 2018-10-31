from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem
from proteus import WaveTools as wt

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time",1.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",4,"Refinement level"),
    ("grav", [0, -9.81, 0], "Gravity vector in m/s^2"),
    # waves
    ("water_level", 1., "Water level from y=0"),
    ("wave_period", 1.0, "Period of the waves in s"),
    ("wave_height", 0.025, "Height of the waves in m"),
    ("depth", 1., "Wave depth in m"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wavelength", 3., "Wavelength in m"),
    ])

###################
# DOMAIN AND MESH #
###################
L=(25.0,1.5)
coord = [0,0,0]
refinement = opts.refinement
domain = RectangularDomain(L=L,x=coord)
inlet_x_position = coord[0]
outlet_x_position = L[0]


# CREATE REFINEMENT #
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1

he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#############################################
###    Monochromiatic Waves (linear)      ###
#############################################

# general options
waterLevel = opts.water_level

# waves
period = opts.wave_period
height = opts.wave_height
mwl = opts.water_level
depth = opts.depth
direction = opts.wave_dir
# define wave as the Monochromiatic wave from WaveTools
wave = wt.MonochromaticWaves(period, height, mwl, depth, np.array(opts.grav), direction)

#####################################
# DEFINE SOME NEEDED FUNCTIONS HERE #
#####################################
def bathymetry(X):
    x = X[0]
    # assume flat bathymetry for now
    return 0.0*x

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        # here we initialize with water height = depth
        return depth
class momU_at_t0(object):
    def uOfXT(self,X,t):
        return 0.0
class momV_at_t0(object):
    def uOfXT(self,X,t):
        return 0.0

###############################
##### BOUNDARY CONDITIONS #####
###############################
def water_height_DBC(X,flag):
    if X[0]==outlet_x_position:
        # here we set outlet to have same depth as initial depth
        return lambda x,t: depth
def x_mom_DBC(X,flag):
    if X[0]==inlet_x_position:
        # here we set momentum (waterHeight * vel_x) flow at the inlet
        return lambda X,t: depth * wave.u(X,t)[0]
def y_mom_DBC(X,flag):
        # no vel_y flow for now
        return lambda x,t: 0.0



# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': momU_at_t0(),
                     'y_mom': momV_at_t0()}
boundaryConditions = {'water_height': water_height_DBC,
                      'x_mom': x_mom_DBC,
                      'y_mom': y_mom_DBC}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=0,
                                              cfl=0.33,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=0
mySWFlowProblem.physical_parameters['mannings']=0
