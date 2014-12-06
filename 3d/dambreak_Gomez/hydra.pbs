#!/bin/sh

#PBS -S /bin/sh
#PBS -N cas1215_Dambreak_Gomez
#PBS -o dambreak.output.o   # stdout file
#PBS -e dambreak.output.e   # stderr file
#PBS -l nodes=1:ppn=10 # nodes required / processors per node
#PBS -q highp          # queue name

source /etc/profile.d/modules.sh


module load proteus/0.9.0/
cd $PBS_O_WORKDIR

export NCPUS="10"
mpirun -n $NCPUS  parun dambreak_Gomez_so.py -m -p -l 7 -v -O /apps/proteus/0.9.0/proteus-mprans/benchmarks/inputTemplates/petsc.options.superlu_dist -D dambreak_Gomez.$PBS_JOBID | tee log.$PBS_JOBID
