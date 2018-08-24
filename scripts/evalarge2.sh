#!/bin/bash
binFile=$1
matrixPath=""
bigMat=$2

ch=2
dr=6
thBlas=1
th=6
l=0
echo "Running $binFile for the set in $matrixPath"
    for f in ${matrixPath}*.mtx; do
        echo $f;
       	for l in {2,1,0,-1,-2,-3,-4,-5}; do
            $binFile $f $th $ch $th $l $thBlas $dr ${bigMat}${f}.ord
            echo ""
            #thread
    	done
    done
