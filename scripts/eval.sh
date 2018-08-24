#!/bin/bash
binFile=$1
matrixPath=$2
bigMat=$3
echo "Running $binFile for the set in $matrixPath"
for i in {1..3}; do
    for f in ${matrixPath}*.mtx; do
    	for th in {1..7}; do
		for ch in {3,5,7,10,15}; do
    			for cost in {100,10000}; do
    				for l in {1,4}; do
    					$binFile $f $th $ch $cost $l
        				echo ""
    				done
    			done
    			#level
		done
		#chunck
    	done
    	#thread
    done
done
