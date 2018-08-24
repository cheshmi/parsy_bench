#!/bin/sh
#SBATCH -A TG-ASC160007
#SBATCH --job-name="edgeCutCholesky"
#SBATCH --output="edgeCutCholesky.%j.%N.out"
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
ch=1
dr=7
thBlas=1
th=12
ip=12
echo "Running $binFile for the edge cut analysis in $matrixPath"
    for f in ${matrixPath}*.mtx; do
        echo "#  $f";
        for l in {-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9}; do
            $binFile $f $th $ch $ip $l $thBlas $dr
            echo ""
        done
    done
