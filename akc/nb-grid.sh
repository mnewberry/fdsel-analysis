#!/bin/bash

export SEARCH=$1
export SIM_ARGS=$2

function run {
SS=$1
MU=$2
TAG=nb$SEARCH-s$SS-mu$MU

if [ `echo "$SS * $MU > 0.00005" | bc -l` == "1" ] 
then echo Skipping simulation s=$SS mu=$MU that would take forever
else 
cat > out/params_$TAG <<EOF
param	ind	val	var	lci	uci
mu	NA	$MU	0.0	NA	NA
ne	NA	1.	0.	1.	1.
srep	NA	0.	NA	0.	0.
s1	0	$SS	0.0	NA	NA
EOF
echo Rscript sim-nb.R \
  -o out/wplog_$TAG.tsv \
  -t out/timeseries_$TAG \
  -p out/params_$TAG \
  -s $SS -u $MU \
  $SIM_ARGS
Rscript sim-nb.R \
  -o out/wplog_$TAG.tsv \
  -t out/timeseries_$TAG \
  -p out/params_$TAG \
  -s $SS -u $MU \
  $SIM_ARGS
fi
}

export -f run

echo " ================== Coarse grid search ===================== "
SSS="0.0000001 0.0000002 0.0000005
  0.000001 0.000002 0.000005
  0.00001 0.00002 0.00005
  0.0001 0.0002 0.0005
  0.001 0.002 0.005
  0.01 0.02 0.05
  0.1 0.2 0.5"

MUS="0.0000001 0.0000002 0.0000005
  0.000001 0.000002 0.000005
  0.00001 0.00002 0.00005
  0.0001 0.0002 0.0005
  0.001 0.002 0.005"

parallel run {1} {2} ::: $SSS ::: $MUS || exit 1

echo " ================== Fine grid search ======================= "
SSS="0.00000011 0.00000016  0.00000025 0.0000004 0.0000006
  0.0000011 0.0000016 0.0000025 0.000004 0.000006
  0.000011 0.000016 0.000025 0.00004 0.00006
  0.00011 0.00016 0.00025 0.0004 0.0006
  0.0011 0.0016 0.0025 0.004 0.006"

MUS="0.0000011 0.0000016 0.0000025 0.000004 0.000006
  0.000011 0.000016 0.000025 0.00004 0.00006
  0.00011 0.00016 0.00025 0.0004 0.0006
  0.0011 0.0016 0.0025 0.004 0.006"

parallel run {1} {2} ::: $SSS ::: $MUS || exit 1

if [ $SEARCH == "cpld" ]
then echo " ================== Super Fine ======================= "
SSS="0.00008 0.00009 0.00012 0.00014 0.00017 0.00019 0.00021"
MUS="0.00003 0.00005 0.00007 0.00008 0.00009 0.00012 0.00014
  0.00017 0.00019 0.00021"

parallel run {1} {2} ::: $SSS ::: $MUS || exit 1
fi
