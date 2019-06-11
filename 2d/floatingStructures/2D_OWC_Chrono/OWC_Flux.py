#!/usr/bin/env python
#Generate a Time series of water height level at the OWC midpoint
from paraview.simple import *
import csv
import os
#from gatherParaviewOutput import *
#import gatherParaviewOutput as AT
#import numpy as np ##NUMPY INCOMPATIBLE WITH PVPYTHON

def gatherParaviewOutput(filename,timesteps,variable):
    """
    Collect Data from the numbered .csv files output by Paraview's save data for all time steps
    """
    Array = [0]*timesteps
    #Array = np.zeros((timesteps))
    for i in range(timesteps):
        if filename[-4:] is '.csv':
	    filename_str = filename[:-4] + '.' + str(i) + '.csv'
	else:
	    filename_str = filename + '.' + str(i) + '.csv'
        with open(filename_str, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Array[i] = float(row[variable])
    return(Array)


#Generate Paraview Output
def GenerateWaterLevel(file_location):
    #Load State
    files = [(file_location + "floating2D.xmf")]
    reader = OpenDataFile(files)
    #Fluid
    clip = Clip(Input=reader)
    clip.ClipType = 'Scalar'
    clip.Scalars = ('POINTS', 'vof')
    clip.Value = 0.5
    clip.InsideOut = 1           #We want to visualise the water (vof > 0.5)
    
    # Interrogate record.csv for Barycentr (Not robust)
    csv_file = file_location + "custom1.csv"
    with open(csv_file, 'rb') as csvfile:
	reader = csv.DictReader(csvfile)
        rownum = 0
	for row in reader:
	    if rownum == 1:
                Barycentre_X = float(row['x'])
	        Barycentre_Y = float(row['y'])
            rownum += 1
    ##
    print("Read .csv file, Barycenter obtained")

    #Inside OWC clip
    clip2 = Clip()
    clip2.ClipType = 'Cylinder'
    clip2.ClipType.Center = [Barycentre_X,Barycentre_Y,0.]
    clip2.ClipType.Axis = [0,1,0]
    clip2.ClipType.Radius = 0.17
    clip2.InsideOut = 1

    #Calculate base slice of OWC
    #If Amplitude is 0.0025, water height is 1.03, and device is 0.18 below FS:
    slice_y = 1.03 - 0.18 + 0.02       # ensure fully inside the device 
    
    #Mid_OWC Column Slice
    slice = Slice()
    slice.SliceType = 'Plane'
    slice.SliceType.Origin = [Barycentre_X,slice_y,0.]
    slice.SliceType.Normal = [0,1,0]
    
    ## Integrate velocity over slice, and calculate the water flux. Then the Water Height from this.
    
    Integrator = IntegrateVariables(slice)
    #What is the tank diameter
    Tank_Diameter = 0.18    
    Calculator1=Calculator(Integrator)
    Fun_Str = "v/"+str(Tank_Diameter)
    Calculator1.Function = Fun_Str
    Calculator1.ResultArrayName = "Flux"
    
    #Calculator2=Calculator(Integrator)
    #Fun_Str_2 = "1.03+v"
    #Calculator2.Function = Fun_Str_2
    #Calculator2.ResultArrayName = "WaterHeight"
    ##Set up views for inspection
    #Set up views/plots
    SaveState(file_location + "Auto_Flux_Calc.pvsm")
    print("Saved Paraview State")

    #Export csv files
    #SaveData("test.csv", Calculator1, Precision = 2, WriteAllTimeSteps = True)    #OLD DO NOT USE
    SetActiveSource(Calculator1)
    print("Saving .csv files of data")

    #SAVE OUTPUT DATA IN A DIRECTORY TO ADD
    SaveData(file_location + "test.csv",proxy=Calculator1,WriteAllTimeSteps=1)
    
## RUN
file_location = "Trials/20_05_SpringTest/1_204/"
# gmsh, gmsh_ref
print("Running")
GenerateWaterLevel(file_location)

#Count TimeSteps
fileCounter = 0
for root, dirs, files in os.walk(file_location):
    for file in files:
        if file.startswith('test.'):
            fileCounter += 1

filename = file_location + 'test'

#print("expecting ", fileCounter, " .csv files to parse")
#Pull Data together
WaterLevel = gatherParaviewOutput(filename,fileCounter,'Flux')

#Remove the evidence
print("Deleting temporary paraview output files")
for i in range(fileCounter):
    if filename[-4:] is '.csv':
        filename_str = filename[:-4] + '.' + str(i) + '.csv'
    else:
        filename_str = filename + '.' + str(i) + '.csv'
    os.remove(filename_str)

#print(WaterLevel)

#TIME FOUND MANUALLY FOR THIS SET

#WRITE A .CSV OF THE DATA
with open(file_location + "Flux.csv", 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(WaterLevel)
