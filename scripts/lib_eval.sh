#!/bin/sh

binFile=$1
matrixPath=$2
bigMat=$2
ch=1
l=1
thBlas=4
div=7
th=4

echo "Running $binFile for the set in $matrixPath"
    for f in ${matrixPath}*.mtx; do
        echo "Time spent in factorization step (numfct) $f";
        #for th in {4}; do
		 #$binFile  $f $th 5 ${bigMat}${f}.ord
		 $binFile  $f $th 5
		 echo ""
    	#done
    done
