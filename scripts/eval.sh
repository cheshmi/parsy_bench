#!/bin/bash
binFile=$1
matrixPath=$2
bigMat=$3
ch=1
dr=4
thBlas=1
th=4
l=-1
echo "Running $binFile for the set in $matrixPath"
    for f in ${matrixPath}*.mtx; do

        echo "#  $f";
       	#for ip in {6,5}; do
            for l in {2,1,0,-1,-2}; do
                for dr in {2,4}; do
                    $binFile $f $th $ch $th $l $thBlas $dr
                    echo ""
                done
    	    done
    	#done
    done
