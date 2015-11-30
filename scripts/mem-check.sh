#!/bin/bash



# TOL_LIST+=( 1e-05 1e-06 1e-07 1e-08 1e-09 1e-10 )
TOL_LIST+=( 1e-11 1e-12 1e-13 1e-14 1e-16 )
# Q_LIST+=( 04 06 14 )
Q_LIST+=( 14 )

for q in "${Q_LIST[@]}"
do

for tol in "${TOL_LIST[@]}"
do
# ibrun -n 1 -o 0 tacc_affinity valgrind --tool=massif --massif-out-file=valgrind-output-${tol} /home1/03237/ga29hoj/code/tbslas/examples/bin/advection -N 64.0 -tol ${tol} -d 15 -q 6 -cuf 2 -tn 1 -dt 1e-9 -vs 1 -merge 3 -omp 16 -cubic 1 > max_oct_tol${tol}.out

ibrun -n 1 -o 0 tacc_affinity /home1/03237/ga29hoj/code/tbslas/examples/bin/advection -N 64.0 -tol ${tol} -d 15 -q ${q} -cuf 2 -tn 1 -dt 1e-9 -vs 1 -merge 3 -omp 16 -cubic 1 > max_oct_q${q}_tol${tol}.out
done

done