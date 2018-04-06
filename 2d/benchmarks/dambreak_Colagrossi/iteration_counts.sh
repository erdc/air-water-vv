#!/bin/bash
echo "Iteration Statistics" > stats.txt
for f in output_files/out_selfp_boomerAMG_*.all; do echo $f >> stats.txt; grep "Linear rans2p_ solve converged" $f >> stats.txt; done
cat stats.txt
./avg.py stats.txt
