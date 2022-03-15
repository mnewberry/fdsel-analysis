#!/bin/bash

OUT=out/timeseries-netherlands

echo -e "type\tgen\tcount" > ${OUT}.tsv

zcat inp/corpus_1stName_NL.csv.gz > inp/corpus_1stName_NL.csv
cat inp/corpus_1stName_NL.csv | \
  grep -v ';2015;' | \
  sed -e 's/"M";/M/g' | \
  sed -e 's/"V";/F/g' | \
  sed -e 's/;/\t/g' \
  >> ${OUT}.tsv

cat ${OUT}.tsv | grep -v ^M > ${OUT}F.tsv
cat ${OUT}.tsv | grep -v ^F > ${OUT}M.tsv

# === Simulate timeseries at inferred parameters ===
fdsel timeseries -i out/timeseries-netherlands.tsv -p - \
  | cut -f 2 | tail +2 | sed 's/^/5*/' \
  | bc -l > out/netherlands5x.pop_sizes.nlsv

Rscript - <<EOM
options(stringsAsFactors = FALSE)
source("../common.R")
ts = read.delim("out/timeseries-netherlands.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlands.typefreqs.tsv")
ts10y = subset(ts, gen >= 1950 & gen < 2010)
ts10y = transform(ts10y, y10 = (gen - 1950) %/% 10)
ts10y = transform(ts10y, y10type = paste(type, y10), y10gen = y10 * 10 + 1950)
ts10y = transform(ts10y, y10count = ave(count, y10type, FUN=sum))
ts10y = unique(ts10y[,c("y10gen","type","y10count")])
names(ts10y) = c("gen","type","count")
write.tsv(ts10y, "out/timeseries-netherlands10y.tsv")
ts5y = subset(ts, gen >= 1950)
ts5y = transform(ts5y, y5 = (gen - 1950) %/% 5)
ts5y = transform(ts5y, y5type = paste(type, y5), y5gen = y5 * 5 + 1950)
ts5y = transform(ts5y, y5count = ave(count, y5type, FUN=sum))
ts5y = unique(ts5y[,c("y5gen","type","y5count")])
names(ts5y) = c("gen","type","count")
write.tsv(ts5y, "out/timeseries-netherlands5y.tsv")
ts2y = subset(ts, gen >= 1950 & gen < 2014)
ts2y = transform(ts2y, y2 = (gen - 1950) %/% 2)
ts2y = transform(ts2y, y2type = paste(type, y2), y2gen = y2 * 2 + 1950)
ts2y = transform(ts2y, y2count = ave(count, y2type, FUN=sum))
ts2y = unique(ts2y[,c("y2gen","type","y2count")])
names(ts2y) = c("gen","type","count")
write.tsv(ts2y, "out/timeseries-netherlands2y.tsv")
EOM

# The inp params file below comes from manually adjusting the timescale
# relative to the 5yl25 inference. An ne line must be present but is ignored.
# > params = read.delim("inp/inf-netherlands5yl25.params.tsv")
# > params = transform(params,val = val/5)
# > write.tsv(params, "inp/inf-netherlands5yl25.params.tsv")
fdsel simulate -i inp/inf-netherlands5yl25.params.tsv \
  -B inp/inf-netherlands5yl25.breaks.tsv \
  -T out/timeseries-netherlands5y.tsv \
  -o out/timeseries-netherlands5xsim.tsv \
  -d out/netherlands5x.pop_sizes.nlsv

# cat out/inf-netherlands2B.pop_sizes.tsv | tail +2 | cut -f 2 > inp/inf-netherlands2B.pop_sizes.nlsv

fdsel simulate -b 100 -i inp/inf-netherlands5yl25.params.tsv \
  -B inp/inf-netherlands5yl25.breaks.tsv \
  -o out/timeseries-netherlandssim5yp1ydM.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

fdsel simulate -b 100 -i inp/inf-netherlands5yl25.params.tsv \
  -B inp/inf-netherlands5yl25.breaks.tsv \
  -T out/timeseries-netherlands.tsv \
  -o out/timeseries-netherlandssim5yp1yd.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

Rscript - <<EOM
options(stringsAsFactors = FALSE)
source("../common.R")
ts = read.delim("out/timeseries-netherlandssim5yp1ydM.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandssim5yp1ydM.typefreqs.tsv")
ts = read.delim("out/timeseries-netherlandssim5yp1yd.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandssim5yp1yd.typefreqs.tsv")
EOM

fdsel simulate -b 5000 -i inp/inf-netherlands5yl25mu.params.tsv \
  -B inp/inf-netherlands5yl25.breaks.tsv \
  -o out/timeseries-netherlandsmusim5yp1ydM.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

fdsel simulate -b 5000 -i inp/inf-netherlands5yl25mu.params.tsv \
  -B inp/inf-netherlands5yl25.breaks.tsv \
  -T out/timeseries-netherlands.tsv \
  -o out/timeseries-netherlandsmusim5yp1yd.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

Rscript - <<EOM
options(stringsAsFactors = FALSE)
source("../common.R")
ts = read.delim("out/timeseries-netherlandsmusim5yp1ydM.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandsmusim5yp1ydM.typefreqs.tsv")
ts = read.delim("out/timeseries-netherlandsmusim5yp1yd.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandsmusim5yp1yd.typefreqs.tsv")
EOM

fdsel simulate -b 10000 -K -i inp/inf-netherlands5yl25mu.params.tsv \
  -o out/timeseries-netherlandsneusim5yp1ydM.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

fdsel simulate -b 10000 -K -i inp/inf-netherlands5yl25mu.params.tsv \
  -T out/timeseries-netherlands.tsv \
  -o out/timeseries-netherlandsneusim5yp1yd.tsv \
  -d inp/inf-netherlands2B.pop_sizes.nlsv

Rscript - <<EOM
options(stringsAsFactors = FALSE)
source("../common.R")
ts = read.delim("out/timeseries-netherlandsneusim5yp1ydM.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandsneusim5yp1ydM.typefreqs.tsv")
ts = read.delim("out/timeseries-netherlandsneusim5yp1yd.tsv")
typedist = make_type_dist(ts)
write.tsv(typedist, "out/netherlandsneusim5yp1yd.typefreqs.tsv")
EOM
