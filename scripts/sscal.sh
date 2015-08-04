#!/bin/bash
WORK_DIR=${PWD}
cd ${TBSLAS_DIR}/examples
# make clean
make -j
cd ${WORK_DIR}

# NODES=$1
# MPI_NUM_PROCESSES=$2
OMP_NUM_THREADS=$1
TOTAL_TIME=$2

#JOB_LIST=(advection diffusion advdiff zalesak cubic)
NP_LIST+=(1   2   4   8   16   32   64   128   256   )
# NP_LIST+=(            1         8        16    )

for np in "${NP_LIST[@]}"
do
./.submit_python_job.sh sscal.py ${np} ${np} $OMP_NUM_THREADS $TOTAL_TIME
sleep 10 
done
