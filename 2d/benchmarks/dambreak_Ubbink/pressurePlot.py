from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put relative path below
filename='pointGauge_pressure.csv'

# Reading file
with open (filename, 'rb') as csvfile:
    data=csv.reader(csvfile, delimiter=",")
    a=[]
    time=[]
    probes=[]
    nRows=0

    for row in data:
# Time steps
        if nRows!=0:
            time.append(float(row[0]))

# Probes location ######################
        if nRows==0:                   #
            for i in row:              #
                if i!= '      time':   #
                    i=float(i[14:24])  #
                    probes.append(i)   #
########################################

        row2=[]
        for j in row:
            if j!= '      time' and nRows>0.:
                j=float(j)
                row2.append(j)
        a.append(row2)
        nRows+=1

#####################################################################################

#Choose to plot in time or in space
    choose=(raw_input('Plotting pressure in time or in space? (ENTER t or s) : '))

    if choose=='s':  
# Choose which time step to plot  
        pressure=[]
        print('Number of rows : '+ str(nRows))
        line = int(raw_input('Enter which line (time step) to plot : '))
        pressure=a[line-1][1:]  
 # Plot pressure in space
        import matplotlib.pyplot as plt
        plt.plot(probes,pressure)
        plt.xlabel('probes location [m]')    
        plt.ylabel('pressure[Pa]')
        plt.suptitle('timestep = %d [sec]' %(float(time[line-2])))
	plt.grid(True)
        plt.show()
        savefig('pressure_in_space_%d.png' %(line))

    if choose=='t':    
# Choose which probes to plot  
        print('Number of probes : '+ str(len(probes)))
        x = int(raw_input('Enter which probes to plot (range from 1 to number of probes: ')) 
        pressure2=[]
        for k in range(1,nRows):
            pressure2.append(a[k][x])    
# Plot pressure in time
	import matplotlib.pyplot as plt
	plt.plot(time,pressure2)
	plt.xlabel('time step [sec]')    
	plt.ylabel('pressure [Pa]')
	plt.suptitle('probe = %d [m]' %(float(probes[x-1])))
        #plt.ylim((4000,6500))
	plt.grid(True)
        plt.show()
        savefig('pressure_in_time_%d.png' %(x))

#####################################################################################

# Print an output file
    info = open('probes.txt','w')
    string1=str(probes)
    string2=string1.replace('[',' ')
    string3=string2.replace(']',' ')   
    string4=string3.replace(',','\t') 
    info.write(str('x')+'\t'+string4+'\n')
    for j in range(1,nRows):
        string5=str(a[j])
        string6=string5.replace('[',' ')
        string7=string6.replace(']',' ')   
        string8=string7.replace(',','\t')   
        info.write(string8+'\n')
    info.close()




