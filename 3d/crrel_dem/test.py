from DEM import DEM
import numpy as np
dem=DEM("two_particles.yaml")
F=np.zeros((2,3),'d')
M=np.zeros((2,3),'d')
dt=1.0
dem.step(F,M,dt)
dem.step(F,M,dt)
import pdb
pdb.set_trace()
