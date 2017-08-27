#!/bin/bash
parun broad_crested_weir_so.py -l 5 -v -C "he=0.1 cfl=0.9 T=5.0" -O petsc/petsc.options.schur -D pcd_r2
parun broad_crested_weir_so.py -l 5 -v -C "he=0.05 cfl=0.9 T=5.0" -O petsc/petsc.options.schur -D pcd_r3
parun broad_crested_weir_so.py -l 5 -v -C "he=0.025 cfl=0.9 T=5.0" -O petsc/petsc.options.schur -D pcd_r4
parun broad_crested_weir_so.py -l 5 -v -C "he=0.0125 cfl=0.9 T=5.0" -O petsc/petsc.options.schur -D pcd_r5
parun broad_crested_weir_so.py -l 5 -v -C "he=0.1 cfl=0.9 T=5.0 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r2
parun broad_crested_weir_so.py -l 5 -v -C "he=0.05 cfl=0.9 T=5.0 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r3
parun broad_crested_weir_so.py -l 5 -v -C "he=0.025 cfl=0.9 T=5.0 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r4
parun broad_crested_weir_so.py -l 5 -v -C "he=0.0125 cfl=0.9 T=5.0 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r5
echo "Iteration Statistics" > stats.txt
for f in *_r*/*.log; do echo $f >> stats.txt; grep "converged= True" $f >> stats.txt; grep failed $f; done
cat stats.txt
./avg.py stats.txt
