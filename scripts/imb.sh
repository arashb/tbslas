#!/bin/bash
WORK_DIR=${PWD}
cd ${TBSLAS_DIR}/examples
# make clean
make -j
cd ${WORK_DIR}

# NODES=$1
# MPI_NUM_PROCESSES=$2
OMP_NUM_THREADS=10
TOTAL_TIME=$1

# NODE_LIST=(4 8 16 32)
NODE_LIST=(1 2 4 8 16 32)

for node in "${NODE_LIST[@]}"
do
./.submit_python_job.sh imb.py $node $((2*$node)) $OMP_NUM_THREADS $TOTAL_TIME
sleep 10s
done
