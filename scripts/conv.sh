#!/bin/bash
WORK_DIR=${PWD}
cd ${TBSLAS_DIR}/examples
# make clean
make -j
cd ${WORK_DIR}

NODES=$1
MPI_NUM_PROCESSES=$2
OMP_NUM_THREADS=$3
TOTAL_TIME=$4

JOB_LIST=(advection diffusion advdiff zalesak cubic)

for job in "${JOB_LIST[@]}"
do

case $HOSTNAME in
    *stampede*) #stampede.tacc.utexas.edu
    module load fftw3
    module load python
    sbatch -N${NODES}  -n${MPI_NUM_PROCESSES} \
	-p normal \
        --time=${TOTAL_TIME} -J ${job} \
	./.run_python.sh ./conv-${job}.py ${MPI_NUM_PROCESSES} ${OMP_NUM_THREADS}
    ;;
    *maverick*) #maverick.tacc.utexas.edu
    module load python
    sbatch -N${NODES}  -n${MPI_NUM_PROCESSES} \
	-p vis \
        --time=${TOTAL_TIME} -J ${job} \
	./.run_python.sh ./conv-${job}.py ${MPI_NUM_PROCESSES} ${OMP_NUM_THREADS}
    ;;
    *ronaldo*) #ronaldo.ices.utexas.edu
	# TODO
    ;;
    *zico*) #zico.ices.utexas.edu
    module load fftw3
    module load python
    $EXEC $MPI_NUM_PROCESSES $OMP_NUM_THREADS
    ;;
    *) # none of the known machines
esac

done