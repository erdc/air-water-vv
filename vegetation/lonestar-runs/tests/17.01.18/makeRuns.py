import numpy as np
import os
import shutil

params = np.loadtxt("tank_parameters.csv", delimiter=',', skiprows=1)
os.system("mkdir runs")
for i in range(params.shape[0]):
    filepath = "runs/"+ "run" + `i`
    os.system("mkdir " + filepath)
    os.system("cp *.py " + filepath)
    os.system("cp tank.stampede.slurm " + filepath)
    os.system("cp petsc.options.asm " + filepath)
    f = open(filepath + "/context.options",'w')
    f.write("parallel=True ")
    f.write("wave_type='single-peaked' ")
    f.write("gauges=True ")
    f.write("depth=" + `params[i][0]` + " ")
    f.write("wave_height=" + `params[i][1]` + " ")
    f.write("peak_period=" + `params[i][2]` + " ")
    f.write("peak_wavelength=" + `params[i][3]` + " ")
    f.write("use_veg=" + "True ")
    f.write("veg_density=400 ")
    f.write("input_filepath=/scratch/01082/smattis/air-water-vv/2d/waves_vegetation/lonestar-runs/runs/run" + `i` + "/ ")
    f.write("tank_height=" + `params[i][4]`)
    f.close()
    
