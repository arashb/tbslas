#!/bin/bash
#SBATCH -J TESTJOB           # job name
#SBATCH -o TESTJOB.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 256              # total number of mpi tasks requested
#SBATCH -N 256  
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00        # run time (hh:mm:ss) - 1.5 hours
ibrun tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 32768 -tol 1.0 -d 15 -dt 0.5 -tn 1.0 -vs 1 -omp 16 -cubic 1 -cuf 4
