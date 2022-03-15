#!/bin/bash

OUT=out/timeseries-ssa
if [[ ! -e out/ssa ]] ; then mkdir out/ssa ; fi
unzip -o -q -d out/ssa inp/names.zip

echo -e "gen\ttype\tcount" > ${OUT}raw.tsv

for FI in out/ssa/yob*.txt ; do 
  yes `basename -s .txt $FI | sed -e s/yob//g` | \
    head -n `wc -l $FI | cut -d ' ' -f 1` | \
    paste - $FI | \
    sed -e 's///g' | \
    sed -e 's/,\([MF]\),/,\1\t/g' \
    >> ${OUT}raw.tsv
done

cat ${OUT}raw.tsv | grep -v ,M > ${OUT}M.tsv
cat ${OUT}raw.tsv | grep -v ,F > ${OUT}F.tsv

function trunc35 {
LN=`grep -n '^1936' $1.tsv | head -1 | cut -d : -f 1`
head -1 $1.tsv > $2.tsv
tail +$LN $1.tsv >> $2.tsv
}

trunc35 ${OUT}raw ${OUT}d35
trunc35 ${OUT}M ${OUT}Md35
trunc35 ${OUT}F ${OUT}Fd35

function truncCE {
cat $1.tsv | \
  grep -v '	Christop,M	' | \
  grep -v '	Alexandr,F	' | \
  grep -v '	Alexande,M	' | \
  grep -v '	Jacqueli,F	' | \
  grep -v '	Cassandr,F	' | \
  grep -v '	Christia,M	' | \
  grep -v '	Nathanie,M	' | \
  grep -v '	Jacquely,F	' | \
  grep -v '	Infant,F	' > $2.tsv
}

truncCE ${OUT}d35  ${OUT}d35dCE 
truncCE ${OUT}Md35 ${OUT}Md35dCE
truncCE ${OUT}Fd35 ${OUT}Fd35dCE

function write_censored_counts {
fdsel timeseries -i $1.tsv -p out/inf-$3.pop_sizes.tsv
cp $1.tsv $2.tsv
cat out/inf-$3.pop_sizes.tsv \
  | cut -f 2 | tail +2 | sed 's/^/0.05*/' \
  | bc -l \
  | cut -d . -f 1 | sed 's/^/CENSORED\t/' \
  | sed '1 i type\tcount' | paste out/inf-$3.pop_sizes.tsv - \
  | cut -f 1,3,4 | tail +2 >> $2.tsv
}

write_censored_counts ${OUT}d35dCE ${OUT}Cd35dCE ssad35dCE
write_censored_counts ${OUT}Md35dCE ${OUT}MCd35dCE ssaMd35dCE
write_censored_counts ${OUT}Fd35dCE ${OUT}FCd35dCE ssaFd35dCE
if [[ -e out/timeseries-ssaKd35dCE.tsv ]]
then write_censored_counts ${OUT}Kd35dCE ${OUT}KCd35dCE ssad35dCE
fi
