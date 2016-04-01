#!/bin/bash
#SBATCH -J TESTJOB           # job name
#SBATCH -o TESTJOB.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 8              # total number of mpi tasks requested
#SBATCH -N 8  
#SBATCH -p vis     # queue (partition) -- normal, development, etc.
#SBATCH -t 04:00:00        # run time (hh:mm:ss) - 1.5 hours
ibrun -n 8 -o 0 tacc_affinity /home/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 512 -tol 1e-5 -q 14 -d 15 -dt 0.00625 -tn 160 -test 5 -omp 20 -merge 3 -diff 1e-3 -adap 0 -vsr 10

