from numpy import *
#from scipy import *
from pylab import *
import collections as cll
import scipy.signal as signal
## Developped by Aggelos Dimakopoulos HR Wallingford Ltd
## Revised by Ian Chandler HR Wallingford Ltd

# Put relative path below
filename='combined_gauge_0_0.5_sample_all.txt'



#Reading probe data file
print('Reading data')

# Put header number
headerline = 1

# Read header info 
fid = open(filename,'r')

prdata=loadtxt(fid,skiprows=headerline,delimiter=",")

fid.seek(0)
D = fid.readlines()
header = D[:headerline]

del D

ndt = len(prdata)


# Extracting probe info
probes=(header[0].rstrip("\n")).split(",")
probes=probes[1:]
tot_probes = len(probes)

nprobes=0
probes_p = []
probes_u = []
probes_v = []
#getting pressure probes
indx = []
ii=-1
for probe in probes:
    ii+=1
    if "p" in probe:
        indx.append(ii)
        nprobes+=1
        probes_p.append([])
        temp = probe[:-1].split()
        for s in temp:
            try:
                coord = float(s)
                probes_p[nprobes-1].append(coord)
            except:
                continue
probex = zeros(nprobes,float)
probey = zeros(nprobes,float)
probez = zeros(nprobes,float)
          
#Writing probe info to test file
print('Writing pressure probe information')

info = open('FSEprobe_info.txt','w')

info.write('Found '+ str(nprobes) +' probes in '+filename + '\n') 
info.write('Probe Locations:\n \t x \t y \t z \n' )


ii=-1
for coor_p in probes_p:
    ii+=1
    info.write('\t' + str(coor_p[0]) + '\t' + str(coor_p[1]) + '\t' + str(coor_p[2])+'\n')
    probex[ii] = coor_p[0]
    probey[ii] = coor_p[1]
    probez[ii] = coor_p[2]

ndt = len(prdata)
info.write('Found '+ str(ndt) +' timesteps ' + '\n')
info.close()





#Preparing data for zero crossing analysis
# Storing time
time = zeros(ndt,float)
time[:] = np.array(prdata[:,0])

# Storing fs elevation
prdata =prdata[:,1:]
p_data = zeros((ndt,nprobes),float)
ii = -1
for i in indx:
    ii+=1
    p_data[:,ii] = prdata[:,i]

del prdata    
prdata = p_data













#   Start time of analysis
stime = raw_input('Enter start time for analysis (press Enter for default 0): ')
try:
    stime = float(stime)    
    ja = min(find(time>=stime))
except:
    stime = 0
    ja = 0
print('Start time entered: '+str(stime))
#   end time of analysis. Preferable to put integer number of periods
etime = raw_input('Enter end time for analysis (preferably an integer number of periods, press Enter for default end time): ')
try:
    etime = float(etime)    
    jb = min(find(time>=etime))
except:
    etime = time[-1]
    jb = len(time)-1
print('End time entered: '+str(etime))

pmin = raw_input('Enter wave period for smoothing: ')
pmin = float(pmin)


# Demeaning the time series and storing wave setup
meanh=zeros(nprobes,float)

for j in range(nprobes):
    meanh[j] = average(prdata[ja:jb,j])
for j in range(nprobes):
    prdata[:,j] =prdata[:,j]- meanh[j]

#Filtering
data1 = zeros(shape(prdata),float)    
data1[:,:] = prdata[:,:]
    
for jj in range(nprobes):
    for kk in range(ja,jb):
    
        try:        
            plow = where(time[:]<time[kk]-pmin)[0][-1]
            phigh = where(time[:]>time[kk]+pmin)[0][0]
            prdata[kk,jj] = mean(data1[plow:phigh,jj])
        except:
            prdata[kk,jj] = data1[kk,jj]
           


plot(time,prdata[:,0],"k--",lw=3)
plot(time,data1[:,0])
show()
         


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
s_i = 1e-6
#   External looping for each probe position
for j in range(nprobes):
    upc.append([])
    tint = 0.
    #   Internal looping for scanning the probe time series
    for i in range(ndt-1):
        #   Start and end time condition
        if time[i]>stime and time[i]<etime:
            a1 = prdata[i,j] 
            a2 = prdata[i+1,j]
            # Upcrossing detection
            if  sign(a1*a2)<0 and sign(a1)<0:
                izeros +=1
                upc[j].append(int(i))
                if izeros==1:
                    s_i = (prdata[i+1,j]-prdata[i,j])/(time[i+1]-time[i])
                    s_i = -prdata[i,j]/s_i+time[i]
                s = (prdata[i+1,j]-prdata[i,j])/(time[i+1]-time[i])
                s = -prdata[i,j]/s+time[i]
                tint = time[i]
    Tavg[j] = (s-s_i)/float(izeros-1)
    ncrosses[j] = izeros
    izeros = 0
    if abs(float(probex[j]))<1e-05:
         T1 = Tavg[j]  

prdata[:,:] = data1[:,:]

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


# Put wavelenght here
Lo = float(raw_input("Enter wavelength: "))
kw = 2*pi/Lo
H0 = float(raw_input("Enter waveheight: "))
d0 = float(raw_input("Enter depth: "))
Ap = array(raw_input("Enter [z in (0, -depth), density] of water for pressure probes (press enter for vof probes): "))
try:
    Ap = cosh(kw(Ap[0]+depth))/cosh(kw*depth)/9.81/Ap[1]
except:
    Ap = 1.
H0 = 0.025
kw = 2*pi/Lo
Lp = Lo/6.
x_1 = float(probex[0])
Dx_1 = float(probex[1])-float(probex[0])
for kk in range(1,nprobes):
    x=float(probex[kk])    
    Dx = x-x_1
    if(round(Dx,3) != round(Dx_1,3)): print "PROBES ARE NOT EQUIDISTANT"
    x_1 = x
    Dx_1 = Dx
       
Np = int(round(Lp/Dx))
Np2 = 2*Np

D = kw*Np*Dx
print "Reflection phase: Target 60 Deg, Array %s Deg" %(180.*D/pi)




Lamda = zeros(nprobes,float)
Gamma = zeros(nprobes,float)
Hi = zeros(nprobes,float)
Hr = zeros(nprobes,float)
Amp = Havg/2./9810.

Ap = cosh(kw*(0.5))/cosh(kw*1)

for j in range(nprobes-Np2):
    
    A1 = Amp[j]*Amp[j]
    A2 = Amp[j+Np]*Amp[j+Np]
    A3 = Amp[j+Np2]*Amp[j+Np2]
    Lamda[j] = (A1 + A3 - 2.*A2*cos(2*D))/(4.*sin(D)*sin(D))
    Gamma[j] = 0.5*sqrt(
        ((2*A2-A1-A3)/(2.*sin(D)*sin(D)))**2+((A1-A3)/sin(2*D))**2)
    
Hi = sqrt(Lamda + Gamma) + sqrt(Lamda - Gamma)
Hr = sqrt(Lamda + Gamma) - sqrt(Lamda - Gamma)
Hi[nprobes-1] =1.
Hi[nprobes-2] =1.
Rf = Hr/(Hi+1e-15)


#show()
fig4 = figure(4)
fig4.add_subplot(1,1,1)
line1 = plot(probex[:],100.*Havg/H0/9810/Ap,'ko',lw=2)
line2 = plot(probex[:],100.*Hi/H0/Ap,'k-',lw=2)
line3 = plot(probex[:],100.*Hr/H0/Ap,'k--',lw=2)
line4 = plot(probex[:],100*Rf,'x',lw=2)
fig4.legend((line1,line2,line3,line4), ('H','Hi','Hr','Rf'),'right') 

xlim(5.0,30.0)
ylim(0,110)
grid()
title('Wave height error along the waveflume')
xlabel('$x$ (m)')
ylabel('$H/H_o$ (%)')        

savefig('Waveheight_stats.pdf')
savefig('Waveheight_stats.png')

H1 = Hi[0]




#Writing analysis
info = open('FSEProbe_stats.txt','w')

info.write('Analysis performed between '+ str(stime) +' and ' + str(etime) + '\n')
info.write('Probe No'+'\t'+'x'+'\t\t'+'Period'+'\t\t'+'Mean Height'+'\t\t' + 'Hi'+'\t'+'H1/3'+'\t'+'H1/10'+'\t'+'Hmax'+'\n')
#At FSEProbe_stats.txt it is printed: 1 Probe No / 2 Probe x/ 3 Probe y / 4 Avg Period / 5 Avg H / 6 Input H error / Input T error
for j in range(nprobes):
    info.write(str(j+1)+'\t\t'+str(probex[j])+'\t\t'+str(Tavg[j])+'\t'+str(Havg[j]/9810/Ap,)+'\t' + str(Hi[j])+'\t'+str(H3[j])+'\t'+str(H10[j])+'\t'+str(Hmax[j])+'\t'+'\n')

info.close()


Ref= mean(Rf[0:-41])
fir=open("refl.txt","w")
fir.write(str(Ref))
fir.close


print('Finished zero crossing analysis')


