#from scipy import *
from pylab import *
from numpy import *
import collections as cll
from tank import *
#from math import pi


# Read data
NumberOfProbes=1
NumberOfLines=44
x_loc=[0.5]
y_loc=-0.5 #y=0 in the mean free surface
t_0 = zeros(int(NumberOfProbes),float)
filename=[]
for i in range (0,int(NumberOfProbes)):
 filename.append("probes_data" + str(i) + ".csv")

n=NumberOfLines

for j in range (0,int(NumberOfProbes)):
  #read paraview data
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
  meanp = mean(a[:,5])
  y_p[:]=a[:,5]-meanp
  y_p[:] /= float(rho_0)*9.81*cosh(float(k)*(float(inflowHeightMean)+float(y_loc)))/cosh(float(k)*inflowHeightMean)
  y_phi[:]=-a[:,9]+float(y_loc)
  p=a[:,5]
  u=a[:,6]
  v=a[:,7]
  time = a[:,15]
   
  #read PROTEUS gauge's data
  n=50
  filename_PR= ['combined_gauge_0_0.5_sample_all.txt']
  p_PR = zeros(int(n),float)
  u_PR = zeros(int(n),float)
  v_PR = zeros(int(n),float)
  time_PR = zeros(int(n),float)
  headerline = 1
  fid = open(filename_PR[j],'r')
  fid.seek(0)
  D = fid.readlines()
  header = D[:headerline]
  n=len(D)
  print n
  b = numpy.array(zeros((int(n-1.0),6)))
  for i in range (1,int(n)):
    b[i-1,:]=D[i].split(',')
  b = numpy.array(b,dtype=float32)
  p_PR=b[:,1]
  u_PR=b[:,2]
  v_PR=b[:,4]
  time_PR=b[:,0]

  y_theor=zeros(int(n),float)
  kk = 0
  for ts in time:
    y_theor[kk]=float(waveheight)/2.0*cos(2.0*3.14/float(wavelength)*x_loc[j]-2.0*3.14/float(period) * ts )
    kk+=1


  fig1 = figure()
  fig1.add_subplot(1,1,1)
  plot(time_PR,p_PR,"k-",lw=1.5,label='p_GF')
  plot(time,p,"bo",lw=1,label='p_PV')
  legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
  ##xlim(80,90)
  #ylim(-0.08,0.08)
  grid()
  #title('Free surface elevation at x='+ str(x_loc[j]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
  xticks(fontsize = 14) 
  yticks(fontsize = 14) 
  ylabel(r"$p (N/m^2)$",fontsize=16)
  xlabel(r'$time (s)$', fontsize=16)
  savefig('probe_graph_p' + str(j) +'.png',dpi=100)

  fig2 = figure()
  fig2.add_subplot(1,1,1)
  plot(time_PR,u_PR,"k-",lw=1.5,label='u_GF')
  plot(time,u,"bo",lw=1,label='u_PV')
  plot(time_PR,v_PR,"k-",lw=1.5,label='v_GF')
  plot(time,v,"bo",lw=1,label='v_PV')
  legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
  #ylim(-0.08,0.08)
  grid()
  #title('Free surface elevation at x='+ str(x_loc[j]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
  xticks(fontsize = 14) 
  yticks(fontsize = 14) 
  ylabel(r"$velocity (m/s)$",fontsize=16)
  xlabel(r'$time (s)$', fontsize=16)
  savefig('probe_graph_vel' + str(j) +'.png',dpi=100)
