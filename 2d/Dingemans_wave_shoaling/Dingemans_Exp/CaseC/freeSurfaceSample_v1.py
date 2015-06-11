
import numpy
from scipy import *
from pylab import *
import collections as cll
import os

#Sampling locations

fide1 = open('expdata/bec01gh.dat','r')
fide2 = open('expdata/bec02gh.dat','r')

exper1 = loadtxt(fide1,skiprows=7)
#dim = len(exper1)
#exp = zeros((dim,12),float)

#'''
#data[:,0] = data[:,0]/sqrt(2)
#data[:,1:] = data[:,1:]/2.
#exper[:,0] = exper[:,0]/sqrt(2)
#exper[:,1:] = exper[:,1:]/(2)

#data2[:,0] = data2[:,0]/sqrt(2)
#data2[:,1:] = data2[:,1:]/2.
#'''



L = (54.25,1.26)
GenerationZoneLength = 4.0
AbsorptionZoneLength= 6.0
xSponge = GenerationZoneLength
xRelaxCenter = xSponge/2.0
xSponge_2 = L[0]-AbsorptionZoneLength
xRelaxCenter_2 = 0.5*(xSponge_2+L[0])

x_loc=[0.0,xRelaxCenter,xSponge,24.04,30.04,34.04,xSponge_2,xRelaxCenter_2,L[0]]



NumLocations=int(len(x_loc))
NumCol=3*NumLocations +1
y_loc=-0.2 #y=0 in the mean free surface

#Read data from PROTEUS
filename_VOF=[]
filename_VOF.append('../../non_linear_wave_shoaling_GAZ_65_EPS15_A/column_gauge.csv')

             
fid = open(filename_VOF[0],'r')
fid.seek(0)
headerline = 1
D = fid.readlines()
header = D[:headerline]
n=len(D)-1
y_VOF = zeros((int(n),int(NumLocations)),float)
for ss in range (0,int(NumLocations)):
    time_P = zeros(int(n),float)
    fid = open(filename_VOF[0],'r')
    fid.seek(0)
    headerline = 1
    D = fid.readlines()
    header = D[:headerline]
    n=len(D)
    print n
    b3 = numpy.array(zeros((int(n-1.0),int(NumLocations+1))))
    for i in range (1,int(n)):
      b3[i-1,:]=D[i].split(',')
    b3 = numpy.array(b3,dtype=float32)
    y_VOF[:,ss]=b3[:,ss+1]
    y_VOF[:,ss]-=mean(y_VOF[:,ss])
    y_VOF[:,ss]=-y_VOF[:,ss]
    time_P = b3[:,0]


#fid2 = open('setup02_take2/surfaceElevation/0.0029994/surfaceElevation.dat','r')
#data2=loadtxt(fid2,skiprows=4)




fig1 = figure(1,figsize=(7,21),dpi=100)

#Probe 1
fig1.add_subplot(3,1,1)
for t in gca().get_yticklabels():
    t.set_fontsize(14)
for t in gca().get_xticklabels():
    t.set_fontsize(14)
#line1= plot(data[:,0]+1.97,data[:,6]-0.9,'-b')
line= plot(exper1[:,0],exper1[:,3],'xk',markersize=2)
line2=plot(time_P[:]+10.1,y_VOF[:,3],'-k',lw=1) 
#line2= plot(data2[:,0]+1.86,data2[:,4]-0.9,'-k',lw=1)
xlim(55,61)
ylim(-0.07,0.15)
grid()

text(55.1,0.13,'x=20.04 m',size=15)
figl=figlegend((line,line2),('Dingemans (1994)','Numerical model'),'center right',
          bbox_to_anchor = (0.9, 0.95),ncol=2,prop={'size':14},frameon=False)

#Probe 2
fig1.add_subplot(3,1,2)
for t in gca().get_yticklabels():
    t.set_fontsize(14)
for t in gca().get_xticklabels():
    t.set_fontsize(14)
#line1= plot(data[:,0]+1.97,data[:,6]-0.9,'-b')
line= plot(exper1[:,0],exper1[:,4],'xk',markersize=2)
line2=plot(time_P[:]+10.1,y_VOF[:,4],'-k',lw=1)
#line2= plot(data2[:,0]+1.82,data2[:,6]-0.9,'-k',lw=1)
ylabel('$\eta$ (m)',size=20)
grid()
text(60.1,0.13,'x=26.04 m',size=15)

xlim(60,66)
ylim(-0.07,0.15)


#Probe 3
fig1.add_subplot(3,1,3)
for t in gca().get_yticklabels():
    t.set_fontsize(13)
for t in gca().get_xticklabels():
    t.set_fontsize(13)
#line1= plot(data[:,0]+1.97,data[:,8]-0.9,'-b')
line= plot(exper1[:,0],exper1[:,5],'xk',markersize=2)
line2=plot(time_P[:]+10.1,y_VOF[:,5],'-k',lw=1)
#line2= plot(data2[:,0]+1.84,data2[:,8]-0.9,'-k',lw=1)
text(63.1,0.08,'x=30.04 m',size=15)

xlim(63,69)
ylim(-0.1,0.1)


grid()
#title(r'Dynamic pressure at $x=$' + str(probex[int(j)])+'m')
xlabel('$t$ (s)',size=20)

#savefig('a.pdf')
subplots_adjust(top=0.9,bottom=0.15,right=0.95,left=0.2,wspace=0.5,hspace=0.5)
savefig('Dingemans.png',dpi=300)              
show()

