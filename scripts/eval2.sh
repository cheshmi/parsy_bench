#!/bin/bash
binFile=$1
matrixPath=$2
bigMat=$3
ch=2
echo "Running $binFile for the set in $matrixPath"
    for f in ${matrixPath}*.mtx; do
        echo $f;
        for th in {6,5,4,3,2,1}; do
            for thBlas in {6,1,4,2}; do
                    for l in {1..4}; do
                        for cost in {2,3,4,5,6,12}; do
                                $binFile $f $th $ch $cost $l $thBlas
                                echo ""
                        done
                    done
                    #level
            done
            #thread
    	done
    done
