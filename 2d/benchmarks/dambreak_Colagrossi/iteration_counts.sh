#!/bin/bash
echo "Iteration Statistics" > stats.txt
for f in output_files/out_selfp_boomerAMG_*.all; do echo $f >> stats.txt; grep "Linear rans2p_ solve" $f >> stats.txt; grep "Linear rans2p_fieldsplit_velocity_ solve" $f >> stats.txt; grep "Time (sec):" $f >> stats.txt; grep "Flops:  " $f >> stats.txt;grep "PCSetUp " $f >> stats.txt; grep "PCApply" $f >> stats.txt;  grep "KSPSolve" $f >> stats.txt; done
echo "LaTeX formatted results:" > latex_results.txt
./avg.py stats.txt > results.txt
./latex_tables.py latex_results.txt > latex_tables.txt
cat results.txt
