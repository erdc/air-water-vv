import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

def getSplines(filepath):

    vals = np.loadtxt(filepath+"/vegZoneVelocity.csv",
                      skiprows = 1,
                      delimiter=',')

    f = open(filepath + "/vegZoneVelocity.csv", 'r')
    first = f.readline()
    f.close()
    g = first.split(',')
    zListU=[]
    zListW=[]
    for i in range(1,len(g)):
        h = g[i].split()
        if h[0] == 'u':
            zListU.append(float(h[3]))
        elif h[0] == 'v':
            zListW.append(float(h[3]))

    time = vals[:,0]
    zU = np.array(zListU)-zListU[0]
    zW = np.array(zListW)-zListW[0]

    Uarray = vals[:,1:len(zListU)+1]
    Warray = vals[:,len(zListU)+1::]

    interpU = interp.RectBivariateSpline(time,zU,Uarray, kx=1,ky=1)
    interpW = interp.RectBivariateSpline(time,zW, Warray, kx=1, ky=1)
    #plt.contourf(tt,zzU,Uarray)
    #plt.imshow(Uarray)
    #plt.show()
    vals = np.loadtxt(filepath+ "/column_gauge.csv",
                      skiprows=1,
                      delimiter=',')


    time = vals[:,0]
    water_height = 1.2472727272727273-vals[:,4]

    interp_phi=interp.interp1d(time,water_height)

    return (interp_time,interp_phi, interpU, interpW)

