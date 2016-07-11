##!/usr/bin/env python
## -*- coding: cp1252 -*-
from pylab import *
import pylab as p 
import numpy

#from scipy.interpolate import interp1d
from scipy import interpolate


Str_1 = loadtxt('BroadWeir.txt')
x = Str_1[:,0]
y = Str_1[:,1]

NumCases = 6
CaseCode = ['18_R6_INF2', '18_R6_INF3','18_R6_INF4','18_R6_INF5','18_R5_INF3','18_R7_INF3' ]
NumCol = 15
NumProfiles=[4,5,4,4,5,12] #paraview files for each time snapshot 
HS_filename=[]

HS_runs=[8 , 9 , 10, 11, 9, 9] #experimental data
SnapTime=[5,10,15,20,25,30] #time snaphots 
NumOfSnapShots=len(SnapTime)
FieldVariable= 1 #1 for phi=0; 2 for vof=0.5

for i in range(0,int(NumCases)):
# read experimental data
 HS_filename.append('Profile_H_S_Run_' + str(int(HS_runs[i])) + '.txt') 
 values_1 = loadtxt(HS_filename[i])
 H1 = values_1[:,0]-3.0
 P1 = values_1[:,1]
 if i==5: 
  NumOfSnapShots=2

 kkk=0
 PR_filename=[]
# read proteus data
 for jj in range(0,int(NumOfSnapShots)):
  NN=0 

  a=[]
  b=[]
  for j in range(0,int(NumProfiles[i])):
   if FieldVariable==1:
    PR_filename.append('ParaviewProfiles_' + CaseCode[i]+'/t_' + str(int(SnapTime[jj])) + '_phi' + str(int(j)) +'.csv') 
   elif FieldVariable==2:
    PR_filename.append('ParaviewProfiles_' + CaseCode[i]+ '/t_' + str(int(SnapTime[jj])) + '_vof' + str(int(j)) +'.csv') 
   fid = open(PR_filename[kkk],'r')
   kkk+=1
   fid.seek(0)
   headerline = 1
   D = fid.readlines()
   header = D[:headerline]
   n=len(D)
   b1 = numpy.array(zeros((int(n-1.0),int(NumCol))))
   for ii in range (1,int(n)):
     b1[ii-1,:]=D[ii].split(',')
   b1 = numpy.array(b1,dtype=float32)
   for kk in range (0,int(n-1)):
     a.append(b1[kk,12])
     b.append(b1[kk,13])
   NN+=n-1
   #print NN
   #print 'Snapshot'
   #print jj
   #print 'Profile'
   #print j
   #print a
 
  #print b
  #print NN
  H2=zeros(int(NN),float)
  P2=zeros(int(NN),float)
  H2[:] = a 
  P2[:] = b
  #print H2
  #print P2
  fig1 = figure(1,figsize = (8,14),dpi = 25)


#------------------------------------------
#PROFILE 
#------------------------------------------
  
  
  fig1.add_subplot(6,1,int(jj+1))
  line1=plot(H1,P1,'k-',ms=1, lw=1.0, alpha=1.0, mfc='black')
  line2=plot(H2,P2,'kd',ms=2.0, lw=0.5, alpha=1.0, mfc='black')
  lineS=plot(x,y,'k-',ms=1.5, lw=2, alpha=1.0, mfc='black')

  p.axis([1.4,2.2,0.2,0.65])
  ax = p.gca()
  ax.set_autoscale_on(False)
  maxlx   = MultipleLocator(0.1)      #majorLocatorx
  minlx   = MultipleLocator(0.01)     #minorLocatorx  
  maxly   = MultipleLocator(0.1)      #majorLocatory
  minly   = MultipleLocator(0.01)     #minorLocatory 
  ax.xaxis.set_major_locator(maxlx)
  ax.xaxis.set_minor_locator(minlx)
  ax.yaxis.set_major_locator(maxly)
  ax.yaxis.set_minor_locator(minly)

  plt.legend(('H & S - Run'+ str(HS_runs[i]),'Proteus-'+ str(i+1)+'Pr'+ str(jj)),'upper right',fancybox=False)
  leg = plt.gca().get_legend()
  ltext  = leg.get_texts()
  llines = leg.get_lines()  
  frame  = leg.get_frame()

  #frame.set_facecolor('0.80')    # set the frame face color
  plt.setp(ltext, fontsize=10)    # the legend text fontsize
  plt.setp(llines, linewidth=1) # the legend linewidth
  leg.draw_frame(False)          # don't draw the legend frame
  plt.setp(ltext, fontsize=10)
  #grid(True)

  xlabel('Distance $x[m]$',fontsize=12)
  ylabel('Elevation $y[m]$',fontsize=12)


 subplots_adjust(top=0.9,bottom=0.1,right=0.95,left=0.1,wspace=0.3,hspace=0.3)
 if FieldVariable==1:
   nm='phi'
 elif FieldVariable==2:
  nm='vof'
 savefig('Prof_Comparison_Run' + str(CaseCode[i]) +'_'  + nm +'.pdf')
 plt.close(fig1)
 #show()


#------------------------------------------
#ERROR
#------------------------------------------
 H2new=[]
 P2new = []
 for j in range(len(H2)):
    if H2[j]> H1[0] and H2[j] < H1[len(H1)-1] and P2[j] > 0.400:
       H2new.append(H2[j])  
       P2new.append(P2[j])

 f1 = interpolate.interp1d(H1,P1,kind='cubic')
 Error1 = abs((f1(H2new)-P2new)/(f1(H2new)))*100
 
 AverageError=average(Error1)
 print 'Average Error (%) ' + CaseCode[i]
 print AverageError 

 fig2 = figure(1,figsize = (7,6),dpi = 25)

#-----------------------------------------------------------------------------------

 line1=plot(H2new,Error1,'ks',ms=3, lw=0.5, alpha=1.0, mfc='black',label='Proteus-1')

 #xscale('log')
 yscale('log')

 p.axis([1.0,2.2,0.001,100])
 ax = p.gca()
 #ax.set_autoscale_on(True)


 plt.legend(loc=1,fancybox=False)
 leg = plt.gca().get_legend()
 ltext  = leg.get_texts()
 llines = leg.get_lines()  
 frame  = leg.get_frame()




 ##frame.set_facecolor('0.80')    # set the frame face color
 plt.setp(ltext, fontsize=10)    # the legend text fontsize
 plt.setp(llines, linewidth=1) # the legend linewidth
 leg.draw_frame(False)          # don't draw the legend frame
 plt.setp(ltext, fontsize=10)#,size='large')


 grid(which='both')

 title('Surface elevation error\nRun 6a')
 xlabel('Distance $x[m]$',fontsize=14)
 ylabel('$e_y$ (%)',fontsize=14)

 savefig('ErrorRun_' + str(CaseCode[i]) +'Pr' + str(jj)+'.pdf')
 #show()
 plt.close(fig2)
