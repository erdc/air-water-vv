import os
import numpy as np
import time
fid = open("tests_lcs.csv","r")
ff = fid.readlines()
ff = ff[1:]
for test in ff:
    test = test.replace("[","")
    test = test.replace("]","")
    test = test.replace('"',"")
    test = test.replace('\n',"")
    test = test.split(",")
    test = map(str.strip,test)
    test = map(np.float,test)
    i = int(test[0])-1
    H = test[1]
    T = test[2]
    L = test[3]
    Yc = "["
    Bc = "["
    for ii in range(7):
        Yc += str(test[4+ii])+","
        Bc += str(test[12+ii])+","
    Yc+=str(test[11])+"]"
    Bc+=str(test[19])+"]"
    walltime = '100:00:00'
    so_file = 'tank_so.py'
    petsc_option = '$HOME/air-water-vv/inputTemplates/petsc.options.asm'
    context = 'wave_height=%s wave_period=%s wavelength=%s Ycoeff=%s Bcoeff=%s' %(H,T,L,Yc,Bc)
    folder = '../FOLDER_Test%s_structure3_d_0.315_H_%s_T_%s' %(i+1,H,T)
    pbs_name = 'PBS_LCS_Test%s_struttura3_d_0.315_H_%s_T_%s.pbs' %(i+1,H,T)
    run_name = 'RUN_Test_%s' %(i+1)
    master_folder = os.getcwd()

    try:
        os.system("mkdir %s" %folder)
    except:
        pass
    os.chdir(folder)
    pbs = open(pbs_name, 'w')
    nn=4
    os.system("cp %s/*.py ." %master_folder)
    pbs.write("""#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -q standard
#PBS -N {6}
#PBS -l select={0}:ncpus=32:mpiprocs=32
#PBS -l walltime={1}
#PBS -M asd@hrwallingford.com

cd $PBS_O_WORKDIR/
aprun -n {2} parun {3} -O {4} -l 2 -v -C "{5}" 
""".format(str(nn), str(walltime), str(nn*32), str(so_file), str(petsc_option), str(context),str(run_name)))
    pbs.close()
    os.system("qsub "+pbs_name)
    time.sleep(1)
    os.chdir(master_folder)
