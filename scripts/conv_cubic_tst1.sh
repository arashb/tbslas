#!/bin/bash
# define file array
output_file=conv-cubic-tst1-$(date +%Y%m%d-%H%M%S).out
cuf=(2 4 8 16)
omp=4
EXAMPLES_DIR=../examples/;
cd ${EXAMPLES_DIR};
make -j;
cd -;
EXEC=${EXAMPLES_DIR}/bin/cubic;
# EXEC=./bin/cubic
for f in "${cuf[@]}"
do
    echo "**********************************************************************"
    echo "UPSAMPLING FACTOR: $f "
    ${EXEC} -N 8 -omp $omp -tol 1.0e-08 -cubic 1 -ca 1 -cuf $f
done > "$output_file"
