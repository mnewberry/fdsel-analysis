#!/bin/bash

if [[ ! -e out ]] ; then mkdir out ; fi

# === Download AKC source data from Ghirlanda et al. Figshare ===
wget https://ndownloader.figshare.com/files/1080844 -O inp/AKC.csv
echo Should match:
echo 4fe024dd46ea4629d014e214f7292dd9	inp/AKC.csv
md5sum inp/AKC.csv

# === Convert data to "gen\ttype\tcount" timeseries formats ===
# This script creates out/timeseries.tsv and out/timeseries-*.tsv
# No other script should write to files matching out/timeseries-*.tsv
Rscript timeseries.R
echo Paper timeseries:
echo 00b2db140fbcb15339ea188d5686e3a7 out/timeseries-diyo45.tsv
echo e3c25bc08f8c969e1f4a2e9936a0998a out/timeseries-diyosamples.tsv
echo 2f4fe359e426809f67fcba775288fc50 out/timeseries-diyo.tsv
echo d16122fe66d442acc4a646404be888e2 out/timeseries-diy.tsv
echo 3f38fc4dd6599f88a44bd87968e389dc out/timeseries.tsv
echo Should match generated timeseries:
md5sum out/timeseries-*.tsv out/timeseries.tsv

# === Conduct inferences ===
echo "Inferring raw timeseries..."
fdsel infer -i out/timeseries.tsv -q 10 \
  -p out/inf.pop_sizes.tsv \
  -S out/inf.bootstrap.tsv \
  -b out/inf.breaks.tsv -o out/inf.params.tsv \
  -R out/inf.residuals.tsv \
  -P out/inf.update_residuals.tsv || exit 1

# The "series" is named according to the substring following "-" in the name of
# the timeseries. 
# The "inference" is named according to the concatenation of the series, which
# bins are used (q for quantile and l for logarithmic) and the 0-padded number
# of bins.
# The raw timeseries and the resulting inference have no name.
# Only the inference diyoq10 is actually used in the results.

function infer {
FLAGS=$1
BINNING=$2
NBINS=$3
TAG=$FLAGS$BINNING$NBINS
echo fdsel infer -i out/timeseries-$FLAGS.tsv -$BINNING $NBINS \
  -p out/inf-$TAG.pop_sizes.tsv \
  -S out/inf-$TAG.bootstrap.tsv \
  -b out/inf-$TAG.breaks.tsv \
  -o out/inf-$TAG.params.tsv \
  -R out/inf-$TAG.residuals.tsv \
  -P out/inf-$TAG.update_residuals.tsv \
  -U
fdsel infer -i out/timeseries-$FLAGS.tsv -$BINNING $NBINS \
  -p out/inf-$TAG.pop_sizes.tsv \
  -S out/inf-$TAG.bootstrap.tsv \
  -b out/inf-$TAG.breaks.tsv \
  -o out/inf-$TAG.params.tsv \
  -R out/inf-$TAG.residuals.tsv \
  -P out/inf-$TAG.update_residuals.tsv \
  -U
}

export -f infer

infer diyo q 10

# # Computationally expensive, optional:
# # Infer all series for all reasonable binnings.
#
# # All possible results:
# echo "Inferring timeseries-* timeseries..."
# FLAGS=""
# for FI in out/timeseries-*.tsv ; do 
#   FLAGS="$FLAGS `basename -s .tsv $FI | sed -e s/.*-//`"
# done
# parallel infer {1} {2} {3} ::: $FLAGS ::: q l ::: 06 07 08 09 10 11 12 \
#   || exit 1
# infer diyosamples q 10


# # === Conduct novelty bias grid search ===
# # EXTREMELY COMPUTATIONALLY INTENSIVE (~10 CPU-days)
# # GENERATES LARGE OUTPUT FILES (~2 GB, ~100k files)
# # TYPICALLY UNNECESSARY
# 
# function collect_logs {
# SEARCH=$1
# FN=out/wplog_nb$SEARCH
# echo -n -e "file\t" > $FN.tsv
# cat `ls $FN-s[0-9]* | head -1` | head -1 >> $FN.tsv
# for FI in $FN-s[0-9]* 
# do cat $FI | sed -e s/^/`basename $FI`\\t/ | tail -n +2 >> $FN.tsv
# done
# }
# 
# # census population size with low drift
# ./nb-grid.sh cpld "-B out/inf-diyoq10.breaks.tsv -n 1435737 -D 30"
# collect_logs cpld
# 
# # higher population size with low drift
# ./nb-grid.sh hpld "-B out/inf-diyoq10.breaks.tsv -n 3056784 -D 60"
# collect_logs hpld
# 
# # census population size with high drift
# ./nb-grid.sh cphd "-B out/inf-diyoq10.breaks.tsv -n 1435737 -D 100"
# collect_logs cphd
# 
# # higher population size with high drift
# ./nb-grid.sh hphd "-B out/inf-diyoq10.breaks.tsv -n 3056784 -D 300"
# collect_logs hphd
# 
# # === Simulate best-fit novelty bias model ===
# cat > out/params_nbcpld-s0.00021-mu0.00004 <<EOF
# param	ind	val	var	lci	uci
# mu	NA	0.00004	0.0	NA	NA
# ne	NA	1.	0.	1.	1.
# srep	NA	0.	NA	0.	0.
# s1	0	0.00021	0.0	NA	NA
# EOF
# 
# # simulate preferred parameters for a very long time
# Rscript sim-nb.R -B out/inf-diyoq10.breaks.tsv -n 1435737 -D 30 -m 500 \
#   -p out/params_nbcpld-s0.00021-mu0.00004 \
#   -o out/wplog_nbcpld-long-s0.00021-mu0.00004.tsv \
#   -t out/timeseries_nbcpld-long-s0.00021-mu0.00004 \
#   -s 0.00021 -u 0.00004
# 
# # Write a high-resolution s_ave(p) by stitching the chunks into one long
# # timeseries
# Rscript stitchts.R "out/timeseries_nbcpld-long-s0.00021-mu0.00004-r*.tsv" \
#   "out/timeseries_nbcpld-long-s0.00021-mu0.00004-stitched.tsv" \
#   "out/wp_nbcpld-long-s0.00021-mu0.00004-stitched.tsv" \
#   0.00021 0.00004 50
# 
# # Make an indepenent replicate of the above from one long simulation,
# # thereby testing stitching codepath and independence of 80-gen chunks.
# # simulate the same in one giant run of 500x80 generations
# fdsel simulate -N -s 1 -n 1435737 -D 30 -g 40000 \
#   -o out/timeseries_nbcpld-long-s0.00021-mu0.00004.tsv \
#   -i out/params_nbcpld-s0.00021-mu0.00004
# 
# # Write a high-resolution s_ave(p)
# Rscript - <<END
# source("../common.R")
# ts = read.delim("out/timeseries_nbcpld-long-s0.00021-mu0.00004.tsv")
# minp = min(0.5,ts_gtc_to_gtf(ts)\$freq)
# wp = make_novelty_wp(ts,logit_breaks(1.5*minp,50),0.00021,0.00004)
# write.table(wp, "out/wp_nbcpld-long-s0.00021-mu0.00004.tsv",
#   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
# END
# 
# # Infer s(p) for each chunk and and for the long timeseries
# function runchunk {
# chunk=$1
# fdsel infer -i out/timeseries_nbcpld-long-s0.00021-mu0.00004-r0$chunk.tsv \
#   -B out/inf-diyoq10.breaks.tsv \
#   -o out/inf-nbcpld-long-s0.00021-mu0.00004-r0$chunk.params.tsv \
#   -R out/inf-nbcpld-long-s0.00021-mu0.00004-r0$chunk.residuals.tsv
# }
# 
# export -f runchunk
# seq -w 1 500 | parallel runchunk {}
# 
# # Collect all the s(p)s into a single file
# head -1 out/inf-nbcpld-long-s0.00021-mu0.00004-r0001.params.tsv \
#   | sed -e "s/^/file\\t/" \
#   > out/inf-nbcpld-long-s0.00021-mu0.00004.params.tsv
# for FI in out/inf-nbcpld-long-s0.00021-mu0.00004-r0*.params.tsv ;
# do
# FILE=`basename $FI`
# cat $FI | grep "^s" | sed -e "s/^/$FILE\\t/" \
#   >> out/inf-nbcpld-long-s0.00021-mu0.00004.params.tsv
# done
