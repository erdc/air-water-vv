#from scipy import *
from pylab import *
from numpy import *
#import numpy 
import collections as cll
from tank import *
#from math import pi

# Read data

#Select type of graph
TypeOfGraph=1    #TimeS series=0 ; Profile snapshot=1 
NumberOfLines=101  #number of lines in the output files 101

if TypeOfGraph==0 :
   NumberOfTimeSeries=1
   NumberOfProfiles=1
   it=NumberOfTimeSeries
   filename=[]
   for i in range (0,int(it)):
     filename.append("probes" + str(i) + ".csv") #names of the paraview probe data files 
elif TypeOfGraph==1: 
   NumberOfProfiles=5
   it=NumberOfProfiles
   t=[3.2, 4.0, 8.0, 12.0, 16.0] #times of profiles 
   filename=[]
   for i in range (0,int(it)):
     filename.append("./outputs_solitary_239272/t="+ str(t[i])+ "sy0xline.csv") #names of the paraview probe data files 

#Read proteus data using paraview's probes
n=NumberOfLines
for j in range (0,int(it)):
  phi = zeros(int(n),float)
  z_PhiDimensionless = zeros(int(n),float)
  p = zeros(int(n),float)
  u = zeros(int(n),float)
  v = zeros(int(n),float)
  x = zeros(int(n),float)
  y = zeros(int(n),float)
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

  phi[:]=a[:,9]
  p[:]=a[:,5]
  u[:]=a[:,6]
  v[:]=a[:,7]
  x[:]= a[:,16]
  y[:]=a[:,17]
  y_loc=y[0]-inflowHeightMean #y coordinate y=0 in the mean free surface
  zz=[y_loc]
  z_PhiDimensionless=(-phi+float(y_loc))/waveheight
  n=n-1.0
  if TypeOfGraph==0:
    t = zeros(int(n),float)
    t[:]=a[:,15]
 
  eta=zeros(n,float)
  Uhorz=zeros(n,float)
  Uvert=zeros(n,float)

#Creat indipendant analytical graph-----------
#x=[0]
#z=[-1]
#for i in range(1,NumberOfLines):
#  x.append(x[i-1]+0.1)

#for i in range(1,51):
#  z.append(0.04+z[i-1])

#print z
#t=10.0
  
  jj=0

  h = inflowHeightMean
  H_w= waveheight
  celerity = sqrt(abs(g[1])*(waveheight+inflowHeightMean))

  if TypeOfGraph==0:
   NumberOfTimeIterations=n
   NumberOfXlocations=1
  elif TypeOfGraph==1: 
   NumberOfTimeIterations=NumberOfProfiles
   NumberOfXlocations=n

  
 # for ii in range(0,int(NumberOfTimeIterations)):
  kk=0
  x+=trans 
  for i in range(0,int(NumberOfXlocations)):
    
    try :
      a0=cosh(k*(celerity * t[j] - x[i]))
    except Exception:
      a0=1e300
    try :
      eta[kk]= waveheight/((a0)**2.0)
    except Exception:
      eta[kk]=1e-300
  
    try :
      a1=cosh(sqrt( 3.0 * H_w / h**3.0) * (celerity  * t[j] - x[i]))
    except Exception:
      a1=1e300

    try :
      a2=cosh(k*(celerity * t[j] - x[i]))
    except Exception:
      a2=1e300

    
    try :
      Uhorz[kk] =  1.0 /(4.0 * h**4 ) * sqrt(abs(g[1]) * h) *  H_w  * (             
               2.0 * h**3 + h**2 * H_w  + 12.0 * h * H_w * zz[jj] + 6.0 *  H_w * zz[jj]**2.0 +
              (2.0 * h**3 - h**2 * H_w - 6.0 * h * H_w * zz[jj] - 3.0 * H_w * zz[jj]**2 ) * a1)/(a2)**4
    except Exception:
      Uhorz[kk]=1e-300


    try:
      Uvert[kk] =  -( 1.0 / ( 4.0 * sqrt( abs(g[1])* h) ) * sqrt(3.0) * abs(g[1]) * (H_w / h**3.0)** 1.5  * (h + zz[jj])*(
                2.0 * h**3 - 7.0 * h**2.0 * H_w + 10.0 * h * H_w * zz[jj] + 5.0 * H_w * zz[jj]**2.0 +
                (2.0 * h**3.0 + h**2.0 * H_w - 2.0 * h * H_w * zz[jj] - H_w * zz[jj]**2.0)*
                cosh(sqrt( 3.0 * H_w / h**3.0) * (celerity * t[j] - x[i] )))/(
                cosh(sqrt( 3.0 * H_w / ( 4.0 * h**3.0))*
                (celerity * t[j] - x[i]))   )** 4.0*( 
                tanh( sqrt( 3.0 * H_w / ( 4.0 * h**3.0))*(celerity * t[j] - x[i]))))
    except Exception:
      Uvert[kk]=1e-300
    
    kk+=1

  ind_eta=eta/H_w
  AnalyticalPressure=rho_0*abs(g[1])*(eta+inflowHeightMean)+(L[1]-eta-inflowHeightMean)*rho_1*abs(g[1])
  

  fig1 = figure()
  fig1.add_subplot(1,1,1)

  if TypeOfGraph==0 :
      plot(t,ind_eta,"k-",lw=1.5,label='eta_analitical/Ho')
      plot(t,Uhorz,"b-",lw=1.5,label='u @z='+ str(zz[jj]))
      plot(t,Uvert,"r-",lw=1.5,label='v@z='+ str(zz[jj]))
      plot(t[:], z_PhiDimensionless[:],"k--",lw=1.5,label='eta_(phi)proteus/Ho')
      plot(t,u ,"g--",lw=1.5,label='u_proteus')
      plot(t,v ,"r--",lw=1.5,label='v_proteus')
      #ylabel(r"$\eta (m)$",fontsize=16)
      xlabel(r'$time (s) $', fontsize=16)
      legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
      ylim(-0.2,1.4)
      grid()
      #title('Free surface elevation at x='+ str(x_loc[jj]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
      xticks(fontsize = 14) 
      yticks(fontsize = 14) 
      #savefig('probe_graph_phi' + str(jj) +'.png',dpi=100)
      savefig('TimeseriesGraph @x=' + str(x[0]) + 'y=' + str(y[0]) +'.png')
  elif TypeOfGraph==1 :
      plot(x-trans,AnalyticalPressure/1000,"b--",lw=1.5,label='Pressure Analytical')
      plot(x-trans,p/1000,"k-",lw=1.5,label='Pressure Proteus')
      #plot(x-trans,ind_eta,"k-",lw=1.5,label='eta_analitical/Ho')
      #plot(x-trans,Uhorz,"b-",lw=1.5,label='u @z='+str(zz[jj]))
      #plot(x-trans,Uvert,"r-",lw=1.5,label='v@z='+str(zz[jj]))
      #plot(x[2:]-trans, z_PhiDimensionless[2:],"k--",lw=1.5,label='eta_(phi)proteus/Ho')
      ylabel(r"$p (KN/m^2)$",fontsize=16)
      xlabel(r'$x (m) $', fontsize=16)
      legend(bbox_to_anchor=[0.99,0.99],ncol=2,fontsize=13)
      ylim(9.75,10.35)
      grid()
      #title('Free surface elevation at x='+ str(x_loc[jj]) + ', y='+ str(y_loc) + '(m)',fontsize=16)
      xticks(fontsize = 14) 
      yticks(fontsize = 14) 
      #savefig('probe_graph_phi' + str(jj) +'.png',dpi=100)
      savefig('TimeSnaphotGraph @z=' + str(zz[jj]) + 't=' + str(t[j]) +'.png')

  


