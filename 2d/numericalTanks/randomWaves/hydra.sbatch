#!/bin/bash
#SBATCH -e caisson.output.e    #imb.%j.e
#SBATCH -p highp #name of the queue
#SBATCH -J cas1234_movC        #job_name
#SBATCH -N 1                   #num of nodes
#SBATCH --ntasks-per-node 12  #cores per node

source /etc/profile.d/modules.sh
export PROTEUS=~/PROTEUS/proteus
export PROTEUS_ARCH=centos
export PROTEUS_PREFIX=$PROTEUS/$PROTEUS_ARCH
export PROTEUS_MPRANS=~/PROTEUS/proteus-mprans

export PATH=$PROTEUS/$PROTEUS_ARCH/bin:$PROTEUS_MPRANS/scripts:$PATH
export LD_LIBRARY_PATH=$PROTEUS/$PROTEUS_ARCH/lib:$LD_LIBRARY_PATH
export INCLUDE_PATH=$PROTEUS/$PROTEUS_ARCH/include:$INCLUDE_PATH
module load gcc/6.1.0 openmpi/1.10.2/gcc/6.1.0

export NCPUS="12"
./batchScript
