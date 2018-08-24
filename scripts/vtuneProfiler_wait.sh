#!/bin/bash
binFile=$1
matrixPath=$2
reportPath=$3
cur=""

ch=1
dr=6
thBlas=1
th=6
l=-1

source /home/kazem/programs/intelStudio/vtune_amplifier/amplxe-vars.sh;

for f in ${cur}*.mtx; do
	for i in {1..5}; do
	amplxe-cl -collect locksandwaits -r something  $binFile $f $th $ch $th $l $thBlas $dr > ${reportPath}$f$i.txt;

	amplxe-cl -report top-down -r something/something.amplxe -column=time:total -report-output ${reportPath}$f$i.csv;

	rm -r something;
done
done
