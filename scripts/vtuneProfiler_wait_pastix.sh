#!/bin/bash
binFile=$1
reportPath=$2
cur=""

ch=1
dr=6
thBlas=1
th=6
l=-2
source /home/kazem/programs/intel/mkl/bin/mklvars.sh intel64;
source /home/kazem/programs/intelStudio/vtune_amplifier/amplxe-vars.sh;
for f in ${cur}*.mtx; do
for i in {1..5}; do
	amplxe-cl -collect locksandwaits -r something  $binFile -mm $f -t $th > ${reportPath}$f$i.txt;

	amplxe-cl -report top-down -r something/something.amplxe -column=time:total -report-output ${reportPath}$f$i.csv;

	rm -r something;
done
done
