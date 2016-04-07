#!/bin/bash
#SBATCH -J TESTJOB           # job name
#SBATCH -o TESTJOB.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1024              # total number of mpi tasks requested
#SBATCH -N 1024  
#SBATCH -p large     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00        # run time (hh:mm:ss) - 1.5 hours
ibrun -n 1024 -o 0 tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 2048 -tol 1e-05 -d 15 -q 14 -cuf 4 -tn 1 -dt 0.00625 -merge 3 -omp 16 -vsr 0 -test 7 -diff 0.001 -ea 160 -m 8
# ibrun -n 1024 -o 0 tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 2048 -tol 1e-05 -d 15 -q 14 -cuf 4 -tn 1 -dt 0.00625 -merge 3 -omp 16 -vsr 0 -test 7 -diff 0.001 -ea 160 -m 8
# ibrun -n 8 -o 0 tacc_affinity /home/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 512 -tol 1e-5 -q 14 -d 15 -dt 0.00625 -tn 160 -test 5 -omp 20 -merge 3 -diff 1e-3 -adap 0 -vsr 10

