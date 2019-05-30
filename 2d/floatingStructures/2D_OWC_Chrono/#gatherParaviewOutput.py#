#!/usr/bin/env python

#Pulls Timestep data into an array for further manipulation.
import csv
import numpy as np
#from builtins import range

def gatherParaviewOutput(filename,timesteps,variable):
    """
    Collect Data from the numbered .csv files output by Paraview's save data for all time steps
    """
    Array = np.zeros(timesteps)
    for i in range(timesteps):
        if timesteps[-4:] is '.csv':
	    filename_str = filename[:-4] + '.' + str(i) + '.csv'
	else:
	    filename_str = filename + '.' + str(i) + '.csv'
        with open(filename_str, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Array[i] = float(row[variable]) 
    return(Array)

#gatherParaviewOutput('Paraview_Data',301,'water')
