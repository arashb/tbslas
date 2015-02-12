#!/bin/bash

JOB=$1
NODES=$2
MPI_NUM_PROCESSES=$3
OMP_NUM_THREADS=$4
TOTAL_TIME=$5

case $HOSTNAME in
    *stampede*) #stampede.tacc.utexas.edu
    module load fftw3
    module load python
    sbatch -N${NODES}  -n${MPI_NUM_PROCESSES} \
	-p normal \
	-o ${JOB}-$(date +%Y%m%d-%H%M%S).out \
        --time=${TOTAL_TIME} -J ${JOB} \
	./.run_python.sh ${JOB} ${MPI_NUM_PROCESSES} ${OMP_NUM_THREADS}
    ;;
    *maverick*) #maverick.tacc.utexas.edu
    module load python
    sbatch -N${NODES}  -n${MPI_NUM_PROCESSES} \
	-p vis \
	-o ${JOB}-$(date +%Y%m%d-%H%M%S).out \
        --time=${TOTAL_TIME} -J ${JOB} \
	./.run_python.sh ${JOB} ${MPI_NUM_PROCESSES} ${OMP_NUM_THREADS}
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
