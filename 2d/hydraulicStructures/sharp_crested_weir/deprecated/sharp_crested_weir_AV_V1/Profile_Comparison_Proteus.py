from __future__ import print_function
from __future__ import division
##!/usr/bin/env python
## -*- coding: cp1252 -*-
from builtins import str
from builtins import range
from past.utils import old_div
from pylab import *
import pylab as p 
import numpy

#from scipy.interpolate import interp1d
from scipy import interpolate


Str_1 = loadtxt('SharpWeir.txt')
x = Str_1[:,0]
y = Str_1[:,1]

NumCases=1
CaseCode = ['18_R4_INF2']#, '18_R5_INF1']
NumCol=15
NumProfiles=[6]#,8]
HS1_filename=[]
HS2_filename=[]

HS_runs=[ 0.4 , 0.2]#0.3 ,0.4, 0.5, 0.6]
SnapTime=[10]#,10,15,20,25,30] #time snaphots 
NumOfSnapShots=len(SnapTime)
FieldVariable= 1 #1 for phi=0; 2 for vof=0.5



for i in range(0,int(NumCases)):
# read experimental data
 HS1_filename.append('UP_h_' + str(HS_runs[i]) + '_Theory_Mesh_1.txt') 
 values_1 = loadtxt(HS1_filename[i])
 H1 = values_1[:,0]
 P1 = values_1[:,1]

 HS2_filename.append('DOWN_h_' + str(HS_runs[i]) + '_Theory_Mesh_1.txt') 
 values_3 = loadtxt(HS2_filename[i])
 H3 = values_3[:,0]
 P3 = values_3[:,1]
 PR_filename=[]
 kkk=0
 for jj in range(0,int(NumOfSnapShots)):
  a=[]
  b=[]
  # read proteus data
  NN=0 
  for j in range(0,int(NumProfiles[i])):
   if FieldVariable==1:
    PR_filename.append('Paraview_Profiles_' + CaseCode[i]+'/t_' + str(int(SnapTime[jj])) + '_phi' + str(int(j)) +'.csv') 
   elif FieldVariable==2:
    PR_filename.append('Paraview_Profiles_' + CaseCode[i]+ '/t_' + str(int(SnapTime[jj])) + '_vof' + str(int(j)) +'.csv')

   fid = open(PR_filename[kkk],'r')#
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
     a.append(b1[kk,12]-2.5)
     b.append(b1[kk,13])
   NN+=n-1
   print(j)
  #print a
 
  #print b
  #print NN
  H2=zeros(int(NN),float)
  P2=zeros(int(NN),float)
  H2[:] = a 
  P2[:] = b
  #print H2
  #print P2
  fig1 = figure(1,figsize = (8,35),dpi = 25)


#------------------------------------------
#PROFILE 
#------------------------------------------

  fig1.add_subplot(6,1,int(jj+1))
  line1=plot(H1,P1,'k-',ms=1, lw=0.5, alpha=1.0, mfc='black')
  line2=plot(H2,P2,'ko',ms=5.5, lw=0.5, alpha=1.0, mfc='white')
  line3=plot(H3,P3,'k--',ms=1, lw=0.5, alpha=1.0, mfc='black')
  lineS=plot(x,y,'k-',ms=1.5, lw=3, alpha=1.0, mfc='black')

  p.axis([-0.15,0.8,0.5,1.5])
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

  plt.legend(('Theor. Upper WES','Proteus-'+ str(i+1) + 'Pr'+str(jj),'Theor. Lower Montes (1992)'),'upper right',fancybox=True)
  leg = plt.gca().get_legend()
  ltext  = leg.get_texts()
  llines = leg.get_lines()  
  frame  = leg.get_frame()

  #frame.set_facecolor('0.80')    # set the frame face color
  plt.setp(ltext, fontsize=12)    # the legend text fontsize
  plt.setp(llines, linewidth=1) # the legend linewidth
  leg.draw_frame(False)          # don't draw the legend frame
 
  #grid(True)

  xlabel('Distance $x[m]$',fontsize=13)
  ylabel('Elevation $y[m]$',fontsize=13)


 subplots_adjust(top=0.95,bottom=0.1,right=0.95,left=0.1,wspace=0.25,hspace=0.25)
 if FieldVariable==1:
   nm='phi'
 elif FieldVariable==2:
  nm='vof'
 savefig('Prof_Comparison_Run' + str(CaseCode[i]) +'_'  + nm +'.pdf')
 plt.close(fig1)

 #show()

#------------------------------------------
#ERROR LOWER
#------------------------------------------

 H2new=[]
 P2new = []
 deltaP = 0.05
 f1 = interpolate.interp1d(H3,P3,kind='cubic')

 for j in range(len(H2)):
    if H2[j]> H3[0] and H2[j] < H3[len(H3)-1] and P2[j] < (f1(H2[j]) + deltaP) and P2[j] > (f1(H2[j]) - deltaP):
       H2new.append(H2[j])  
       P2new.append(P2[j])

 Error1 = abs(old_div((f1(H2new)-P2new),(f1(H2new))))*100

 AverageError=average(Error1)
 print('Average Error Lower Profile (%)') 
 print(AverageError)
#-----------------------------------------------------
 fig2 = figure(2,figsize = (8,6),dpi = 25)
 line1=plot(H2new,Error1,'ks',ms=5, lw=0.5, alpha=1.0)


 #xscale('log')
 yscale('log')
 xlim(0,0.8)
 #p.axis([-0.5,0.4,0.0,1.30])
 ax = p.gca()
 ax.set_autoscale_on(True)


 grid(which='both')

 title('Test Case 1 - LOWER Profile',fontsize=14)
 xlabel('Distance $x[m]$',fontsize=14)
 ylabel('$E_y$[%]',fontsize=14)

 savefig('SW_ErrorLow_'+str(CaseCode[i]) +'Pr' + str(jj)+ '.png')
 close(fig2)
#------------------------------------------
#ERROR UPPER
#------------------------------------------

 H2new=[]
 P2new = []
 deltaP = 0.1
 f1 = interpolate.interp1d(H1,P1,kind='cubic')

 for j in range(len(H2)):
    if H2[j]> H1[0] and H2[j] < H1[len(H1)-1] and P2[j] < (f1(H2[j]) + deltaP) and P2[j] > (f1(H2[j]) - deltaP):
       H2new.append(H2[j])  
       P2new.append(P2[j])

 Error2 = abs(old_div((f1(H2new)-P2new),(f1(H2new))))*100
 AverageError=average(Error2)
 print('Average Error Upper Profile (%)') 
 print(AverageError)
#----------------------------------------------------
 fig3 = figure(3,figsize = (8,6),dpi = 25)
 line1=plot(H2new,Error2,'ks',ms=5, lw=0.5, alpha=1.0)


 #xscale('log')
 yscale('log')

 xlim(-0.15,0.401)

 ax = p.gca()
 ax.set_autoscale_on(True)


 grid(which='both')

 title('Test Case 1 - Upper Profile',fontsize=14)
 xlabel('Distance $x[m]$',fontsize=14)
 ylabel('$E_y$[%]',fontsize=14)

 savefig('SW_ErrorUp_'+str(CaseCode[i]) +'Pr' + str(jj)+ '.png')
 close(fig3)
