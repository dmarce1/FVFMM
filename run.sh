#!/bin/bash
set -x

cd /work/dmarce1/FVFMM/Release

cp $PBS_NODEFILE hosts
awk 'NR == 0 || NR % 20 == 0' hosts > hosts1
export PBS_NODEFILE=hosts1
NPROC=$(wc -l < hosts1)

export I_MPI_FABRICS=shm:ofa
export I_MPI_DAPL_PROVIDER="ofa-v2-mlx4_0-1u"
export I_MPI_OFA_ADAPTER_NAME=mlx4_0

mpirun -np $NPROC ./FVFMM -t20

