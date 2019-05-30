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
    #Mid_OWC Column Slice
    slice = Slice(Input=reader)
    slice.SliceType = 'Plane'
    ##
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
    print("Read .csv, Barycentre obtained")
    slice.SliceType.Origin = [Barycentre_X,Barycentre_Y,0.]
    slice.SliceType.Normal = [1,0,0]
    ## Integrate vof over slice, and calculate the water height.
    Integrator = IntegrateVariables(slice)
    #What is the tank height???
    Tank_Height = 2.06    
    Calculator1=Calculator(Integrator)
    Fun_Str = str(Tank_Height) + '-vof'
    Calculator1.Function = Fun_Str
    Calculator1.ResultArrayName = "WaterLevel"

    ##Set up views for inspection
    #Set up views/plots
    SaveState(file_location + "Auto_FS_Calc.pvsm")
    print("Saved Paraview State")

    #Export csv files
    #SaveData("test.csv", Calculator1, Precision = 2, WriteAllTimeSteps = True)    #OLD DO NOT USE
    SetActiveSource(Calculator1)
    print("Saving .csv files of data")

    #SAVE OUTPUT DATA IN A DIRECTORY TO ADD
    SaveData(file_location + "test.csv",proxy=Calculator1,WriteAllTimeSteps=1)
    
## RUN
file_location = "Trials/20_03_T_gmsh/1_185/"
# gmsh, gmsh_ref
print("Running")
GenerateWaterLevel(file_location)

#Count TimeSteps
fileCounter = 0
for root, dirs, files in os.walk(file_location):
    for file in files:
        if file.startswith('test.'):
            fileCounter += 1

#print("expecting ", fileCounter, " .csv files to parse")
#Pull Data together
WaterLevel = gatherParaviewOutput(file_location + 'test',fileCounter,'WaterLevel')
#print(WaterLevel)

#TIME FOUND MANUALLY FOR THIS SET

#WRITE A .CSV OF THE DATA
with open(file_location + "Results.csv", 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(WaterLevel)
