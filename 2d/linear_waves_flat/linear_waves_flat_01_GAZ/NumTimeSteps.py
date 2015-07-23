#from scipy import *
#from pylab import *
#import numpy 
#import collections as cll
#from tank import *
#from math import pi

filename_PR=[]
filename_PR.append('combined_gauge_0_0.5_sample_all.txt')
fid = open(filename_PR[0],'r')
fid.seek(0)
headerline = 1
D = fid.readlines()
header = D[:headerline]
n=len(D)-1

print n
 
