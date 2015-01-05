#from scipy import *
from pylab import *
from numpy import *
import collections as cll
from tank import *
#from math import pi


# Read data
NumberOfProbes=6
NumberOfLines=1052
x_loc=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
y_loc=-0.5 #y=0 in the mean free surface
t_0 = zeros(int(NumberOfProbes),float)
filename=[]
for i in range (0,6):
 filename.append("probes_data" + str(i) + ".csv")

n=NumberOfLines

for j in range (0,6):
  # t_0[j]=x_loc[j]/(float(wavelength)/float(period))
  y_p = zeros(int(n),float)
  y_phi = zeros(int(n),float)
  time = zeros(int(n),float)
  headerline = 1
  fid = open(filename[j],'r')
  fid.seek(0)
  D = fid.readlines()
  header = D[:headerline]
  n=len(D)
  print n
  a = numpy.array(zeros((int(n-1.0),19)))
  for i in range (1,int(n)):
    a[i-1,:]=D[i].split(',')
  a = numpy.array(a,dtype=float32)

  # Create arrays to be plotted

  n=n-1.0
  meanp = mean(a[2:,5])
  y_p[2:]=a[2:,5]-meanp
  y_p[2:] /= float(rho_0)*9.81*cosh(float(k)*(float(inflowHeightMean)+float(y_loc)))/cosh(float(k)*inflowHeightMean)
  y_phi[2:]=-a[2:,9]+float(y_loc)
  time[2:]= a[2:,15]
  
 
  y_theor=zeros(int(n),float)
  kk = 0
  for ts in time:
    y_theor[kk]=float(waveheight)/2.0*cos(2.0*3.14/float(wavelength)*x_loc[j]-2.0*3.14/float(period) * ts )
    kk+=1


  fig1 = figure()
  fig1.add_subplot(1,1,1)
  plot(time[:],y_p[:],"k-",lw=1.5,label='Proteus-p')
  plot(time[:],y_theor[:],"r-",lw=1,label='Analytical')
  legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
  ##xlim(80,90)
  ylim(-0.08,0.08)
  grid()
  title('Free surface elevation at x='+ str(x_loc[j]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
  xticks(fontsize = 14) 
  yticks(fontsize = 14) 
  ylabel(r"$\eta (m)$",fontsize=16)
  xlabel(r'$time (s)$', fontsize=16)
  savefig('probe_graph_p' + str(j) +'.png',dpi=100)

  fig2 = figure()
  fig2.add_subplot(1,1,1)
  plot(time[:],y_phi[:],"k-",lw=1.5,label='Proteus-phi')
  plot(time[:],y_theor[:],"r-",lw=1,label='Analytical')
  legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
  ylim(-0.08,0.08)
  grid()
  title('Free surface elevation at x='+ str(x_loc[j]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
  xticks(fontsize = 14) 
  yticks(fontsize = 14) 
  ylabel(r"$\eta (m)$",fontsize=16)
  xlabel(r'$time (s)$', fontsize=16)
  savefig('probe_graph_phi' + str(j) +'.png',dpi=100)
