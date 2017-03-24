#from scipy import *
from pylab import *
from numpy import *
import collections as cll

#Define parameters
L=3.141 #wave length
h=float(L)/2.0 #water depth
m=1.0 #order of the solution 
km_sym=(2.0*m-1.0)*3.141/float(L)
km_asym=2.0*m*3.141/float(L)
am=0.041*2.0 #amplitude
symmetric=True
if symmetric==True:
 km=km_sym
else:
 km=km_asym

w=sqrt(9.81*km*tanh(km*h)) #freq

# Read data

filename='probe_data00000.csv'
headerline = 1
fid = open(filename,'r')
fid.seek(0)
D = fid.readlines()
header = D[:headerline]
n=len(D)
a = array(zeros((int(n-1.0),19)))
for i in range (1,int(n)):
 a[i-1,:]=D[i].split(',')
a = array(a,dtype=float32)

# Create arrays to be plotted
n=n-1.0
y = zeros(int(n),float)
time = zeros(int(n),float)
y_theor=zeros(int(n),float)

y[:] = a[:,9]
time[:]= a[:,15]
x_loc=0.0 #define location
y_loc=h

print time
for i in range (int(n)):
 y_theor[i]=-y_loc-1.0/9.81*am*cos(km*x_loc)*w*cos(w*time[i])*cosh(km*h) 

#j=str(2)

fig1 = figure(1)
fig1.add_subplot(1,1,1)
plot(time,y,"b-",lw=2,label='Proteus')
plot(time,y_theor,"r-",lw=2,label='Theoritical')
legend()

##xlim(80,90)

ylim(-1.64,-1.46)

grid()
#title('Surface elevation at $x=$' + str(probex[int(j)])+'m')
title('Free surface elevation at x=0.0, y=0.0')
xlabel('Time (s)')
ylabel('phi (m)')

show()
savefig('plottrial.pdf')
savefig('plottrial.png',dpi=100)
