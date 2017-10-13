#!/bin/bash
echo "Iteration Statistics" > stats.txt
for f in *_mpi*/*.log; do echo $f >> stats.txt; grep "converged= True" $f >> stats.txt; grep failed $f; done
cat stats.txt
./avg.py stats.txt
