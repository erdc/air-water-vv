#!/bin/bash
parun broad_crested_weir_so.py -l 5 -v -C "he=0.1 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10" -O petsc/petsc.options.schur -D pcd_rf1
parun broad_crested_weir_so.py -l 5 -v -C "he=0.05 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10" -O petsc/petsc.options.schur -D pcd_rf2
parun broad_crested_weir_so.py -l 5 -v -C "he=0.025 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10" -O petsc/petsc.options.schur -D pcd_rf3
parun broad_crested_weir_so.py -l 5 -v -C "he=0.0125 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10" -O petsc/petsc.options.schur -D pcd_rf4
parun broad_crested_weir_so.py -l 5 -v -C "he=0.1 dt_fixed=0.01 fixed_step=True T=5.0 dt_out=0.2 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf1
parun broad_crested_weir_so.py -l 5 -v -C "he=0.05 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf2
parun broad_crested_weir_so.py -l 5 -v -C "he=0.025 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf3
parun broad_crested_weir_so.py -l 5 -v -C "he=0.0125 dt_fixed=0.01 fixed_step=True T=5.0 nsave=10 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_rf4
echo "Iteration Statistics" > statsf.txt
for f in *_rf*/*.log; do echo $f >> statsf.txt; grep "converged= True" $f >> statsf.txt; done
cat statsf.txt
./avg.py statsf.txt
