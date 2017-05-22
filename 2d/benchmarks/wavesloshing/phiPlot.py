import numpy as np
from scipy import *
from pylab import *
import collections as cll
import csv

# Put relative path below
filename='pointGauge_levelset.csv'

# Reading file
with open (filename, 'rb') as csvfile:
    data=csv.reader(csvfile, delimiter=",")
    a=[]
    time=[]
    probes=[]
    nRows=0

    for row in data:
# Time steps
        if nRows!=0:
            time.append(float(row[0]))

# Probes location ######################
        if nRows==0:                   #
            for i in row:              #
                if i!= '      time':   #
                    i=float(i[14:24])  #
                    probes.append(i)   #
########################################

        row2=[]
        for j in row:
            if j!= '      time' and nRows>0.:
                j=float(j)
                row2.append(j)
        a.append(row2)
        nRows+=1

#####################################################################################
   
# Choose which probes to plot  
    print('Number of probes : '+ str(len(probes)))
    x = int(raw_input('Enter which probes to plot (range from 1 to number of probes): ')) 
    phi=[]
    for k in range(1,nRows):
        phi.append(a[k][x]+0.05)    
# Plot phi in time
    import matplotlib.pyplot as plt
    plt.plot(time,phi)
    plt.xlabel('time [sec]')    
    plt.ylabel('phi [m]')
    plt.suptitle('Position of the interface at the left boundary plotted against time')
    plt.xlim((0,2.5))
    plt.ylim((0.04,0.06))
    plt.grid(True)
    plt.show()
    savefig('phi_in_time.png')
    
#####################################################################################
    
    k = 31.41592653589793
    h = 0.05
    eps = 0.15707963267948966
    g = (0, -9.81, 0)
    
    def w0(h):
        w_0 = np.sqrt(np.tanh(h))
        return w_0
    
    def w2(w0):
        w_2 = 1./32.*(9.*w0**(-7)-12.*w0**(-3)-3*w0-2*w0**5)
        return w_2
    
    def omega(h, eps):
        w_0 = w0(h)
        w_2 = w2(w_0)
        w = w_0+0.5*eps**2*w_2
        return w

    def eta0(x, t):
        eta_0 = np.sin(t)*np.cos(x)
        return eta_0
    
    def eta1(x, t, w0):
        eta_1 = 1./8.*((w0**2-w0**(-2))+(w0**(-2)-3*w0**(-6))*np.cos(2.*t))*np.cos(2.*x)
        return eta_1
    
    def eta2(x, t, w0):
        b11 = 1./32.*(3.*w0**(-8)+6.*w0**(-4)-5.+2.*w0**4)
        b13 = 3./128.*(9.*w0**(-8)+27.*w0**(-4)-15.+w0**4+2*w0**8)
        b31 = 1./128.*(3.*w0**(-8)+18.*w0**(-4)-5.)
        b33 = 3./128.*(-9.*w0**(-12)+3.*w0**(-8)-3.*w0**(-4)+1)
        eta_2 = b11*np.sin(t)*np.cos(x)+b13*np.sin(t)*np.cos(3*x)+b31*np.sin(3*t)*np.cos(x)+b33*np.sin(3*t)*np.cos(3*x)
        return eta_2
    
    def eps_eta(x, t, h, eps):
        w_0 = w0(h)
        eta_0 = eta0(x, t)
        eta_1 = eta1(x, t, w_0)
        eta_2 = eta2(x, t, w_0)
        epseta = eps*eta_0+eps**2*eta_1+0.5*eps**3*eta_2
        return epseta

    def eta(x, t):
        x_ = x*k
        h_ = h*k
        w_ = omega(h_, eps)
        t_ = (t+2*np.pi/(w_*np.sqrt(k*(-g[1])))*0.25)*(w_*np.sqrt(k*(-g[1])))
        eta = eps_eta(x_, t_, h_, eps)/k
        return eta
    
# Print an output file to validate the results
    Phi_f_Cal = phi[-1] 
    Phi_f_Ana = -eta(0.0,time[len(time)-1])+0.05 
    err = 100*abs(Phi_f_Ana-Phi_f_Cal)/Phi_f_Ana
    val = open('validation_phi.txt', 'w')
    val.write('Gauges taken at the left boundary'+'\n')
    val.write('phi at t=3.0s:'+'\n')
    val.write('Analitical'+'\t'+'Simulation'+'\t'+'Error'+'\n')
    val.write(str(Phi_f_Ana)+'\t'+str(Phi_f_Cal)+'\t'+str(err))
    val.close()
    
#####################################################################################

# Print an output file
    info = open('probes.txt','w')
    string1=str(probes)
    string2=string1.replace('[',' ')
    string3=string2.replace(']',' ')   
    string4=string3.replace(',','\t') 
    info.write(str('x')+'\t'+string4+'\n')
    for j in range(1,nRows):
        string5=str(a[j])
        string6=string5.replace('[',' ')
        string7=string6.replace(']',' ')   
        string8=string7.replace(',','\t')   
        info.write(string8+'\n')
    info.close()




