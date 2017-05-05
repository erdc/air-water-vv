#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -q standard
#PBS -N bDuck
#PBS -l walltime=05:00:00
#PBS -l select=10:ncpus=32:mpiprocs=32
#PBS -M asd@hrwallingford.com
#PBS -m bea
source /opt/modules/default/etc/modules.sh
module unload PrgEnv-pgi
module load  PrgEnv-gnu
module unload cray-libsci
module load acml

cd $PBS_O_WORKDIR
export OUTPUT_FOLDER= $PBS_O_WORKDIR/$PBS_JOBNAME.$PBS_JOBID
mkdir /work/adimako/$PBS_JOBID
aprun -n 320 parun tank3D_so.py -l 2 -O /u/adimako/air-water-vv/inputTemplates/petsc.options.asm -F -D /work/adimako/$PBS_JOBID
