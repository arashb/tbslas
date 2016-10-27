#!/bin/bash
#SBATCH -J TESTJOB           # job name
#SBATCH -o TESTJOB.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -N 32  
#SBATCH -p vis     # queue (partition) -- normal, development, etc.
#SBATCH -t 04:00:00        # run time (hh:mm:ss) - 1.5 hours
# ibrun -n 32 -o 0 tacc_affinity ../examples/bin/ns -N 64.0 -tol 1e-03 -d 4 -q 14 -tn 2 -dt 0.001 -merge 3 -omp 16 -vsr 1 -diff 1e-3 -test 2
# ibrun -n 32 -o 0 tacc_affinity ../examples/bin/ns -N 64.0 -tol 1e-03 -d 15 -q 14 -tn 0 -dt 0.0628 -merge 3 -omp 20 -vsr 1 -diff 1e-3 -test 2

ibrun -n 32 -o 0 tacc_affinity ../examples/bin/ns -N 64 -tol 1e-02 -d 3 -q 14 -tn 100 -dt 0.0628 -merge 3 -omp 20 -vsr 1 -diff 1e-3 -test 3
# ibrun -n 32 -o 0 tacc_affinity ../examples/bin/ns -N 64.0 -tol 1e-03 -d 6 -q 14 -tn 100 -dt 0.0628 -merge 3 -omp 20 -vsr 1 -diff 1e-3 -test 3

# ibrun -n 1024 -o 0 tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 2048 -tol 1e-05 -d 15 -q 14 -cuf 4 -tn 1 -dt 0.00625 -merge 3 -omp 16 -vsr 0 -test 7 -diff 0.001 -ea 160 -m 8
# ibrun -n 1024 -o 0 tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 2048 -tol 1e-05 -d 15 -q 14 -cuf 4 -tn 1 -dt 0.00625 -merge 3 -omp 16 -vsr 0 -test 7 -diff 0.001 -ea 160 -m 8

# ibrun -n 8 -o 0 tacc_affinity /home/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 512 -tol 1e-5 -q 14 -d 15 -dt 0.00625 -tn 160 -test 5 -omp 20 -merge 3 -diff 1e-3 -adap 0 -vsr 10

