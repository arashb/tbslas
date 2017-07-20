#!/bin/bash

ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 2  > k1m02.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 4  > k1m04.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 6  > k1m06.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 8  > k1m08.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 10 > k1m10.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 12 > k1m12.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 14 > k1m14.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 1 -omp 20 -dist 0 -m 16 > k1m16.out

ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 2  > k2m02.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 4  > k2m04.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 6  > k2m06.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 8  > k2m08.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 10 > k2m10.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 12 > k2m12.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 14 > k2m14.out
ibrun -n 8 -o 0 ../examples/bin/kernel_tst -N 1e+6 -ker 2 -omp 20 -dist 0 -m 16 > k2m16.out