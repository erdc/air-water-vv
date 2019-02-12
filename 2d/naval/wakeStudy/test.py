from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np
from proteus.mprans import BodyDynamics as bd
opts=Context.Options([
    # predefined test cases
    ("water_level", 5., "Height of free surface above bottom"),
    # Geometry
    ("tank_dim", (180., 7.5,), "Dimensions of the tank"),
    ("tank_sponge", (1., 0.), "Length of relaxation zones (front/back, left/right)"),
    ("tank_BC", 'FreeSlip', "tank boundary conditions: NoSlip or FreeSlip"),
    # waves
    ('wave', True, 'Enable  generation'),
    ("period", 3.5 , "wave period"),
    ("period_long", 3.5 , "wave period"),
    ("wave_height", 0.5, "wave height"), 
    ("wave_dir", np.array([1.,0.,0.]),"Direction of the waves"),
    # numerical options
    ("refinement_level", 0.0,"he=walength/refinement_level"),
    ("he", 0.2,"he=walength/refinement_level"),
    ("cfl", 0.4,"Target cfl"),
    ("duration", 60., "Durarion of the simulation"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("useVF", 1.0, "For density and viscosity smoothing"),
    ('movingDomain', not True, "Moving domain and mesh option"),
    ('conservativeFlux', True,'Fix post-processing velocity bug for porous interface'),
    # obstacle dimensions (tweak to set up the geometry of the sloping beach)
    ("slope_start", 100, "x coordinate where the beach slope is starting within the tank"), # approx. 3 wave lengths away form the generation
    ("beach_slope", 10, "1/slope"),
    ("beach_crest", 7, "elevation of the highest point of the beach")
    ])


# --- DOMAIN
domain = Domain.PlanarStraightLineGraphDomain()


# --- Phisical constants
rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g =np.array([0.,-9.81,0.])
gAbs=sqrt(sum(g**2))
waterLevel = opts.water_level

# --- WAVE input
period = opts.period
height = opts.wave_height
mwl = opts.water_level
depth = opts.water_level
direction = opts.wave_dir
N = 32
Nwaves = 16
overl = 0.7
cutoff = 0.1
#wave = wt.MonochromaticWaves(period, height, mwl, depth, g, direction)



wave1 = wt.TimeSeries(
    timeSeriesFile="Time_series_short.txt",
    skiprows=0,
    timeSeriesPosition=np.array([0.,0.,0.]),
    depth = opts.water_level,
    N = N,
    mwl=opts.water_level,
    waveDir = opts.wave_dir,
    g =g,
    rec_direct=False,
    window_params = {"Nwaves":Nwaves,"Tm":opts.period,"Window":"costap","Overlap":overl,"Cutoff":cutoff},
    Lgen = np.array([1,0,0])
    )
    
wave2 = wt.TimeSeries(
    timeSeriesFile="Time_series_long.txt",
    skiprows=0,
    timeSeriesPosition=np.array([0.,0.,0.]),
    depth = opts.water_level,
    N = N,
    mwl=opts.water_level,
    waveDir = opts.wave_dir,
    g =g,
    rec_direct=True,
    Lgen = np.array([1,0,0])
    )
    
wave = wave1
time = np.loadtxt("Time_series_short.txt")[:,0]
time = time - time[0]
data = time.copy()
x0 = np.zeros(3)
for i,t in enumerate(time):
    data[i] = wave.eta(x0,t)

np.savetxt("out_short.txt",zip(time,data))

wave = wave2

time = np.loadtxt("Time_series_long.txt")[:,0]
time = time - time[0]
data = time.copy()
x0 = np.zeros(3)
for i,t in enumerate(time):
    data[i] = wave.eta(x0,t)

np.savetxt("out_long.txt",zip(time,data))


wave = wt.CombineWaves([wave1,wave2])

time = np.loadtxt("Time_series_long.txt")[:,0]
time = time - time[0]
data = time.copy()
x0 = np.zeros(3)
for i,t in enumerate(time):
    data[i] = wave.eta(x0,t)

np.savetxt("out_all.txt",zip(time,data))

