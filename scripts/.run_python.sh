#!/bin/bash
PYTHON_SCRIPT=$1
MPI_NUM_PROCESSES=$2
OMP_NUM_THREADS=$3

echo "######################################################################"
echo $PYTHON_SCRIPT
echo "######################################################################"
python ${PYTHON_SCRIPT} ${MPI_NUM_PROCESSES} ${OMP_NUM_THREADS}
