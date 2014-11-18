#from scipy import *
from pylab import *
from numpy import *
from scipy import *
import collections as cll
import re

#Definitions
filename='combined_gauge_0_0.5_sample_all.txt'
headerline = 1

fprobes1 = "probe_0.5_0.5_0.csv"




# Open probe Data

fid = open(filename,'r')
prdata=loadtxt(fid,skiprows=headerline,delimiter=",")
fid.seek(0)
D = fid.readlines()
header = D[:headerline]

#Acquire header info
headerInfo = header[0].split(",")

headerInfo = headerInfo[1:]

#Load data from paraview
paraview1 = loadtxt(open(fprobes1,"r"),skiprows=1,delimiter=",")


# Definitions
ndt=len(prdata)
time = prdata[:,0]
prdata = prdata[:,1:]
nprobes = len(prdata[0,:])
variables = []
coord = []

# Writing header info to lists
for kk in range(nprobes):
    variables.append(headerInfo[2*kk])
    coord.append(headerInfo[2*kk+1])
    

#Plotting
for kkk in range(0,nprobes,2):
    
    var = variables[kkk]
    if(var == "p"):
        line2 = plot(paraview1[:,15],paraview1[:,5],"k--")
        theor = 1.2*9.81*1.2 + 0.1*9.81*1000
        line3 = plot(linspace(0,1,50),ones(50)*theor,"k-.")
    if(var == "u"):
        line2 = plot(paraview1[:,15],paraview1[:,1],"k--")
        theor = 0
        line3 = plot(linspace(0,1,50),ones(50)*theor,"k-.")
    if(var == "v"):
        line2 = plot(paraview1[:,15],paraview1[:,2],"k--")
        theor = 0
        line3 = plot(linspace(0,1,50),ones(50)*theor,"k-.")
    line1 = plot(time,prdata[:,kkk],"k-")
    title("%s (%s %s %s)" %(var,coord[kkk].split()[1],coord[kkk].split()[2],coord[kkk].split()[3][0]))
    xlabel("Time")
    unit = "(m/s)"
    if var =="p":
        unit = "(Pa)"
        lim1 = 950
        lim2 = 1050
    if var =="u" or var=="v":
        unit = "(m/s)"
        lim1 = -0.03
        lim2 = 0.03
    grid()
    figlegend((line1,line2,line3),("Probe function","Paraview","Theoretical"),"upper center",bbox_to_anchor=(0.5,0.9))
    ylabel("%s %s" %(var,unit))
    xlim(0.001,1)
    ylim(lim1,lim2)
    savefig(str(var)+".pdf")
    show()
        
        
    
