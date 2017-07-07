#!/bin/bash 

export LD_LIBRARY_PATH=$PROTEUS/linux2/lib:$LD_LIBRARY_PATH
mpirun -np 4 parun floating2D_so.py -O petsc.options.superlu_dist -l 2 -v -D RUN01_TestAdapt

