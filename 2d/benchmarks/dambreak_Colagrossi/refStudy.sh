#!/bin/bash
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.4" -O petsc/petsc.options.schur -D pcd_r0
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.2" -O petsc/petsc.options.schur -D pcd_r1
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.1" -O petsc/petsc.options.schur -D pcd_r2
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.05" -O petsc/petsc.options.schur -D pcd_r3
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.025" -O petsc/petsc.options.schur -D pcd_r4
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.0125" -O petsc/petsc.options.schur -D pcd_r5
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.4 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r0
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.2 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r1
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.1 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r2
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.05 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r3
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.025 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r4
parun dambreak_Colagrossi_so.py -l 5 -v -C "he=0.0125 schur_solver='selfp_petsc'" -O petsc/petsc.options.schur.selfp_petsc -D selfp_r5
echo "Iteration Statistics" > stats.txt
for f in *_r*/*.log; do echo $f >> stats.txt; grep "converged= True" $f >> stats.txt; done
cat stats.txt
./avg.py stats.txt
