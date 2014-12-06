#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -l walltime=008:00:00
#PBS -l select=2:ncpus=32:mpiprocs=32
#PBS -q standard
#PBS -N sluice2d
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
cd $PBS_O_WORKDIR
mkdir $WORKDIR/sluice_gate.$PBS_JOBID
aprun -n 64  parun sluice_gate_so.py -l 5 -O ../../inputTemplates/petsc.options.asm -D $WORKDIR/sluice_gate.$PBS_JOBID
