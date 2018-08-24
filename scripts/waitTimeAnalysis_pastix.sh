#!/bin/bash
reportPath=$1

# From /sbin/start_udev by Greg KH (GPL v2 only)
# Does $1 contain $2 ?
strstr() {
  [ "${1#*$2*}" = "$1" ] && return 1
  return 0
}

postf1=.txt
postf2=.csv
name=""
srchVar='sy_sopalin_smp '
repSrchVar='sy_sopalin_smp:'
echo "Running $binFile for the set in $reportPath"
for f in ${reportPath}*; do
	actualTime=""
	percentage=""
	replaced=""
	if strstr $f $postf1;  then
		name="${f%$postf1}"
		name="${f#*/}"
		actualTime=$(cat $f | grep 'CPU Time: \|Wait Time: ') 
		
	fi

	if strstr $f $postf2; then
		name="${f%$postf2}"
		name="${f#*/}"
		percentage=$(cat $f | grep 'sy_sopalin_smp ')  ;
		replaced="${percentage/'sy_sopalin_smp'/$repSrchVar}"
	fi
	echo "$name"
	echo "$actualTime"
	echo "$replaced";

done
