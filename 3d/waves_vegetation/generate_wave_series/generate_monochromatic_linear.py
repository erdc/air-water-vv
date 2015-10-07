from proteus import WaveTools
import numpy as np
import math

mwl = 0.5
mw = WaveTools.MonochromaticWaves(period = 1.0,
                                  waveHeight = 0.1,
                                  mwl = mwl,
                                  depth = 1.0,
                                  g = np.array([0.0, 0.0, -9.81]),
                                  waveDir = np.array([1.0, 0.0, 0.0]),
                                  wavelength=None,
                                  waveType="Linear",
                                  Ycoeff = None, 
                                  Bcoeff =None, 
                                  meanVelocity = 0.0,
                                  phi0 = math.pi/2.0)

f = open('monochromatic_linear_time_series.txt', 'w')
num_points = 10000
tList = np.linspace(0.0,100.0,num_points)
etaList = np.zeros(tList.shape)
for i in range(num_points):
    etaList[i] = mw.eta(0.0,0.0,0.0,tList[i]) + mwl

combined = np.vstack((tList,etaList)).transpose()
np.savetxt("monochromatic_linear_time_series.txt", combined)


    
