#from scipy import *
from pylab import *
from numpy import *
import collections as cll
from tank import *



# Read data
NumberOfProbes=1 
NumberOfLines=89 #number of lines in the output files
x_loc=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0] #x coordinate of the probes location 
y_loc=0.0 #y coordinate y=0 in the mean free surface
t_0 = zeros(int(NumberOfProbes),float)
filename=[]
for i in range (0,int(NumberOfProbes)):
 filename.append("prob0" + str(i) + ".csv") #names of the paraview probe data files 

n=NumberOfLines


for j in range (0,int(NumberOfProbes)):

  # Proteus results using paraview's probes ---------------------------------------------
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

  

  n=n-1.0
  y_phi[2:]=-a[2:,9]+float(y_loc)-inflowHeightMean
  time[2:]= a[2:,15]
  
  # Theoritical ---------------------------------------------------------- 
  y_theor=zeros(int(n),float)
  kk = 0
  loc=[x_loc[j],y_loc]
  Y =  [0.04160592,  #Surface elevation Fourier coefficients for non-dimensionalised solution corresponing to the wave input parameters 
         0.00555874,
         0.00065892,
         0.00008144,
         0.00001078,
         0.00000151,    
         0.00000023,     
         0.00000007] 

  order=len(Y)
  for ts in time:
    
    waterDepth = 0
    for i in range(0,int(order)):  
       waterDepth += Y[i]*cos((i+1)*theta(loc,ts))/k
    
    y_theor[kk]=waterDepth
    kk+=1
  y_theor-=mean(y_theor) 
  y_phi-=mean(y_phi)

  #plot---------------------------------------------------------------
  fig1 = figure()
  fig1.add_subplot(1,1,1)
  plot(time[2:],y_phi[2:],"k-",lw=1.5,label='Proteus-phi')
  plot(time[2:],y_theor[2:],"r-",lw=1,label='Analytical')
  legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
  ylim(-0.1,0.1)
  grid()
  title('Free surface elevation at x='+ str(x_loc[j]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
  xticks(fontsize = 14) 
  yticks(fontsize = 14) 
  ylabel(r"$\eta (m)$",fontsize=16)
  xlabel(r'$time (s)$', fontsize=16)
  savefig('probe_graph_phi' + str(j) +'.png',dpi=100)
