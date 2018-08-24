#!/bin/sh
#SBATCH -A TG-ASC160007
#SBATCH --job-name="Sympiler_TRNS_Par"
#SBATCH --output="Sympiler_TRNS_Par.%j.%N.out"
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
bigMat=$2
ch=1
l=1
thBlas=6
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
