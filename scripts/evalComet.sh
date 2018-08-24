#!/bin/sh
#SBATCH -A TG-CCR180004
#SBATCH --job-name="ParSY_DAG_TRNS_Par"
#SBATCH --output="ParSY_DAG_TRNS_Par.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH -t 23:00:00
#SBATCH --mail-user=kazem.cheshmi@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

binFile=$1
matrixPath=$2
bigMat=$3
ch=2
l=1
thBlas=1
th=12
echo "Running $binFile for the set in $matrixPath"
    for f in ${matrixPath}*.mtx; do
        echo $f;
    #    for th in {12}; do
		for l in {2,1,0,-1,-2}; do
			for div in {2,6,7}; do
		            $binFile  $f $th $ch $th $l $thBlas $div
		            echo ""
			done
		done
    #	done
    done
