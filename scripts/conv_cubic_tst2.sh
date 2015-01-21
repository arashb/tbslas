#!/bin/bash
# define file array
output_file=conv-cubic-tst2-$(date +%Y%m%d-%H%M%S).out
cuf=(2 4 8 16)
tol=(1e-0      1e-1      1e-2      1e-3      1e-4      1e-5      1e-6)

omp=4
EXAMPLES_DIR=../examples/;
cd ${EXAMPLES_DIR};
make -j;
cd -;
EXEC=${EXAMPLES_DIR}/bin/cubic;

for t in "${tol[@]}"
do
    echo "======================================================================"
    echo "TOL: $t"
    echo "======================================================================"
    echo "*****************************"
    echo "CHEBYSEHV: "
    echo "*****************************"
    ${EXEC} -N 8 -omp $omp -tol $t
for f in "${cuf[@]}"
do
    echo "*****************************"
    echo "CUBIC: UPSAMPLING FACTOR: $f "
    echo "*****************************"
    ${EXEC} -N 8 -omp $omp -tol $t -cubic 1 -cuf $f
done #> "$output_file"

done > "$output_file"
