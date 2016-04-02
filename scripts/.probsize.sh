#!/bin/bash
#SBATCH -J TESTJOB           # job name
#SBATCH -o TESTJOB.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -N 32  
#SBATCH -p vis     # queue (partition) -- normal, development, etc.
#SBATCH -t 04:00:00        # run time (hh:mm:ss) - 1.5 hours

TOL_LIST+=( 1e-05 )
# TOL_LIST+=( 1e-11 1e-12 1e-13 1e-14 1e-16 )
# Q_LIST+=( 04 06 14 )
Q_LIST+=( 14 )
# A_LIST+=( 10 20 40 80 160 320)
A_LIST+=( 15 30 60 120 240 480)
#A_LIST+=( 1280 )

for q in "${Q_LIST[@]}"
do

for tol in "${TOL_LIST[@]}"
do

for a in "${A_LIST[@]}"
do

ibrun -n 32 -o 0 tacc_affinity /home/03237/ga29hoj/code/tbslas/examples/bin/advdiff-ss -N 512 -tol ${tol} -q ${q} -d 15 -dt 0.00625 -tn 1 -test 6 -omp 20 -merge 3 -diff 1e-3 -adap 0 -vsr 0 -m 8 -ea ${a} > max_oct_${a}.out

done

done

done

