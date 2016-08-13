#!/bin/bash

# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 2  > mod_stokes_q8_m02.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 4  > mod_stokes_q8_m04.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 6  > mod_stokes_q8_m06.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 8  > mod_stokes_q8_m08.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 10 > mod_stokes_q8_m10.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 12 > mod_stokes_q8_m12.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 14 > mod_stokes_q8_m14.out
# ibrun -n 1 -o 0 ../examples//bin/mod_stokes -N 8  -q 8 -tol 1e-4 -omp 20  -m 16 > mod_stokes_q8_m16.out

ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 2  > mod_stokes_q14_m02.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 4  > mod_stokes_q14_m04.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 6  > mod_stokes_q14_m06.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 8  > mod_stokes_q14_m08.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 10 > mod_stokes_q14_m10.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 12 > mod_stokes_q14_m12.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 14 > mod_stokes_q14_m14.out
ibrun -n 1 -o 0 ../examples/bin/mod_stokes -N 8  -q 14 -tol 1e-4 -omp 20  -m 16 > mod_stokes_q14_m16.out

