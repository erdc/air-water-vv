#!/bin/bash
echo "Iteration Statistics" > stats.txt
for f in output_files/out_selfp_boomerAMG_*.all; do echo $f >> stats.txt; grep "converged= True" $f >> stats.txt; done
cat stats.txt
./avg.py stats.txt
