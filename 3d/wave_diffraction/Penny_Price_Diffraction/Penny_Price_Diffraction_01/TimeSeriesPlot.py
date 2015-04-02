#from scipy import *
import pylab 
import numpy 
import collections as cll

# Read data
filename=[]
filename.append("integration_time_series_loc09_upto08s0.csv")
filename.append("integration_time_series_loc09_upto3s00.csv")
filename.append("integration_time_series_loc102_upto08s000.csv")
filename.append("integration_time_series_loc102_upto3s000.csv")
filename.append("data1.5s0.csv")
filename.append("data1.5s0back0.csv")



plotfilename=[]
plotfilename.append("force.png")


legendname=[]
legendname.append("Proteus")
legendname.append("Experiment")

expfilename="ExperimentalData.csv"
headerline = 1
fid = open(expfilename,'r')
fid.seek(0)
D0 = fid.readlines()
header0 = D0[:headerline]
n=len(D0)
print n
F = numpy.zeros((int(n-1.0),2),dtype=float)
for i in range (1,int(n)):
 F[i-1,:]=D0[i].split(',')


for j in range (0,6):
  headerline = 1
  fid = open(filename[j],'r')
  fid.seek(0)
  D = fid.readlines()
  header = D[:headerline]
  n=len(D)
  print n
  print filename[j]
  a = numpy.zeros((int(n-1.0),18),dtype=float)
  for i in range (1,int(n)):
    a[i-1,:]=D[i].split(',')
 # a = numpy.array(a,dtype=float32)

  # Create arrays to be plotted
  n=n-1.0
  if j==0:
   y_1 = numpy.zeros(int(n),dtype=float)
   time_1 = numpy.zeros(int(n),dtype=float)
   y_1[:] = a[:,5]
   time_1[:]= a[:,14]
  elif j==1:  
   y_2 = numpy.zeros(int(n),dtype=float)
   time_2 = numpy.zeros(int(n),dtype=float)
   y_2[:] = a[:,5]
   time_2[:]= a[:,14]
  elif j==2:
   y_3 = numpy.zeros(int(n),dtype=float)
   time_3 = numpy.zeros(int(n),dtype=float)
   y_3[:] = a[:,5]
   time_3[:]= a[:,14] 
  elif j==3: 
   y_4 = numpy.zeros(int(n),dtype=float)
   time_4 = numpy.zeros(int(n),dtype=float)
   y_4[:] = a[:,5]
   time_4[:]= a[:,14]
  elif j==4:
   y_5 = numpy.zeros(int(n),dtype=float)
   time_5 = numpy.zeros(int(n),dtype=float)
   y_5[:] = a[:,5]
   time_5[:]= a[:,14] 
  elif j==5: 
   y_6 = numpy.zeros(int(n),dtype=float)
   time_6 = numpy.zeros(int(n),dtype=float)
   y_6[:] = a[:,5]
   time_6[:]= a[:,14]
time_2[:]=time_2[:]+time_1[-1]
time_4[:]=time_4[:]+time_3[-1]


n1=int(len(y_1)+len(y_2))
n2=int(len(y_3)+len(y_4))
print n1
print n2
y_front = numpy.zeros(int(n1),dtype=float)
y_back = numpy.zeros(int(n2),dtype=float) 
time_front = numpy.zeros(int(n1),dtype=float)
time_back = numpy.zeros(int(n2),dtype=float) 
print y_5
y_front[:]=numpy.append(y_1,y_2)
y_back[:]=numpy.append(y_3,y_4)
time_front[:]=numpy.append(time_1,time_2)
time_back[:]=numpy.append(time_3,time_4)
print time_front
print time_back
netforce=y_front-y_back
  
pylab.fig1 = pylab.figure(figsize=(16, 10))
        
pylab.fig1.add_subplot(1,1,1)
#pylab.plot(time[:],y[:],"k-",lw=2)#,label=legendname[j])
pylab.plot(time_front[:],netforce[:],"k-",lw=2.5,label=legendname[0])
pylab.plot(F[:,0],F[:,1],"bo",lw=0.2,label=legendname[1])
pylab.plot(time_5,y_5-y_6,"k--")
#legend()
pylab.legend(ncol=2,loc=2,fontsize=20)


pylab.ylim(-25.0,55.0)

pylab.grid()
  #title('Surface elevation at $x=$' + str(probex[int(j)])+'m')
  #title('Free surface elevation at x=0.0, y=0.0')
pylab.xticks(fontsize = 18) 
pylab.yticks(fontsize = 18) 
pylab.xlabel(r'$Time (s)$', fontsize=26)
pylab.ylabel(r'$Force (N)$',fontsize=26)

pylab.show()
  
  #savefig(plotfilename[j])
pylab.savefig(plotfilename[0],dpi=200)
