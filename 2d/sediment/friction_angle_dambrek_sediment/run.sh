#!/bin/bash

mkdir $1
cp tank.py $1/tank.py
mpirun -np 4 parun -l 5 -v tank_so.py -O ../../../inputTemplates/petsc.options.asm -D $1