#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -l walltime=001:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -q standard
#PBS -N crump2d
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
cd $PBS_O_WORKDIR
mkdir $WORKDIR/crump_weir.$PBS_JOBID
aprun -n 32  parun crump_weir_so.py -l 5 -O ../../../inputTemplates/petsc.options.asm -D $WORKDIR/crump_weir.$PBS_JOBID
