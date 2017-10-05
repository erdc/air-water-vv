#!/bin/bash
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.4 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf0
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.2 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf1
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.1 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf2
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.05 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf3
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.025 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf4
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.0125 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1" -O petsc/petsc.options.schur -D pcd_rf5
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.4 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.2 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf0
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.2 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.2 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf1
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.1 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.2 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf2
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.05 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf3
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.025 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf4
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.0125 dt_fixed=0.01 fixed_step=True T=1.0 dt_out=0.1 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf5
echo "Iteration Statistics" > statsf.txt
for f in *_rf*/*.log; do echo $f >> statsf.txt; grep "converged= True" $f >> statsf.txt; done
cat statsf.txt
./avg.py statsf.txt
