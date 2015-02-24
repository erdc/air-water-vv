from numpy import *
#from scipy import *
from pylab import *
import collections as cll
from TimeSeriesPlot_v2 import *
## Developped by Aggelos Dimakopoulos HR Wallingford Ltd
## Revised by Ian Chandler HR Wallingford Ltd

nprobes = 2

#Preparing data for zero crossing analysis

Time=zeros(int(len(time_PR[:])),float)
prdata=zeros((int(len(time_PR[:])),int(nprobes)),float)
Time[:]= time_PR[:]
ndt=int(len(time_PR[:]))

probex=[0,25.0,50.0]
prdata[:,0] = y_p_PR[:,1]
prdata[:,1] = y_p_PR[:,2]
print prdata

#   Start time of analysis
stime = raw_input('Enter start time for analysis (press Enter for default 0): ')
try:
    stime = float(stime)    
    ja = min(find(Time>=stime))
except:
    stime = 0
    ja = 0
print('Start time entered: '+str(stime))
#   end time of analysis. Preferable to put integer number of periods
etime = raw_input('Enter end time for analysis (preferably an integer number of periods, press Enter for default end time): ')
try:
    etime = float(etime)    
    jb = min(find(Time>=etime))
except:
    etime = Time[-1]
    jb = len(Time)-1
print('End time entered: '+str(etime))

# Min period accounted for zero crossing
minT = raw_input('Enter min period (press Enter for T=0): ')
try:
    minT = float(minT)    
except:
    minT=0
print('Min period entered: '+str(minT))


# Demeaning the time series and storing wave setup
meanh=zeros(nprobes,float)

for j in range(nprobes):
  meanh[j] = average(prdata[ja:jb,j])
for j in range(nprobes):
  prdata[:,j] =prdata[:,j]- meanh[j]


#Zero crossing analysis
#   definitions ncrosses - > integer array with the number of upcrossings for each probe
#               upc      - > list for indexing position of upcrossings
#               Tavg     - > array that will return the average zero crossing period per probe
ncrosses = zeros(nprobes,int)
Tavg =zeros(nprobes,float)
izeros =0.
upc = []
s = 0.

print('Conducting zero crossing analysis')

# Zero crossing loop

#   External looping for each probe position
for j in range(nprobes):
    upc.append([])
    tint = 0.
    #   Internal looping for scanning the probe time series
    for i in range(ndt-1):
        #   Start and end time condition
        if Time[i]>stime and Time[i]<etime:
            a1 = prdata[i,j] 
            a2 = prdata[i+1,j]
            # Upcrossing detection
            if  sign(a1*a2)<0 and sign(a1)<0 and Time[i]-tint >=minT:
                izeros +=1
                upc[j].append(int(i))
                if izeros==1:
                    s_i = (prdata[i+1,j]-prdata[i,j])/(Time[i+1]-Time[i])
                    s_i = -prdata[i,j]/s_i+Time[i]
                s = (prdata[i+1,j]-prdata[i,j])/(Time[i+1]-Time[i])
                s = -prdata[i,j]/s+Time[i]
                tint = Time[i]
    Tavg[j] = (s-s_i)/float(izeros-1)
    ncrosses[j] = izeros
    izeros = 0
    print j
    if abs(float(probex[j]))<1e-05:
         T1 = Tavg[j]  


indx = []

Havg = zeros(nprobes,float)
Hrms = zeros(nprobes,float)
H3 = zeros(nprobes,float)
H10 = zeros(nprobes,float)
Hmax = zeros(nprobes,float)

#  Average height calculation for each probe. Each period height is calculated between two successive upcrossings
for j in range(nprobes):
    Havg[j] = 0.
    indx = upc[j]
    H = zeros(len(upc[j]),float)
    for i in range(ncrosses[j]-1):
        i1 = indx[i]
        i2 = indx[i+1]
        #print(str(i)+' '+str(i1)+' '+str(i2)+' '+str(j))
        H[i] = (max(prdata[i1:i2,j])-min(prdata[i1:i2,j]))
        Havg[j] = Havg[j] + H[i]
        Hrms[j] = Hrms[j] + (H[i])**2
    Havg[j] = Havg[j]/float(ncrosses[j]-1)
    Hrms[j] = sqrt(Hrms[j]/float(ncrosses[j]-1))
    dum = sorted(H)
    dum1 = dum[int(2.*len(H)/3):]
    dum2 = dum[int(9.*len(H)/10):]
    H3[j] = mean(dum1)
    H10[j] = mean(dum2)
    Hmax[j] = max(H)
    if abs(float(probex[j]))<1e-05:
         H1 = Havg[j]  



#Writing analysis
info = open('FSEProbe_stats.txt','w')

info.write('Analysis performed between '+ str(stime) +' and ' + str(etime) + '\n')
info.write('Probe No'+'\t'+'x'+'\t\t'+'Period'+'\t\t'+'Mean Height'+'\t\t' + 'Hrms'+'\t'+'H1/3'+'\t'+'H1/10'+'\t'+'Hmax'+'\n')
#At FSEProbe_stats.txt it is printed: 1 Probe No / 2 Probe x/ 3 Probe y / 4 Avg Period / 5 Avg H / 6 Input H error / Input T error
for j in range(nprobes):
    info.write(str(j+1)+'\t\t'+str(probex[j])+'\t\t'+str(Tavg[j])+'\t'+str(Havg[j])+'\t' + str(Hrms[j])+'\t'+str(H3[j])+'\t'+str(H10[j])+'\t'+str(Hmax[j])+'\t'+'\n')

info.close()

print('Finished zero crossing analysis')


