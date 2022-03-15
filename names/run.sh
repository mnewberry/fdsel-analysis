#!/bin/bash

if [[ ! -e out ]] ; then mkdir out ; fi
if [[ ! -e inp ]] ; then mkdir inp ; fi

# For subpopulation simulation (Fig S6), follow instructions in
# ./timeseries-subpopsim.sh and uncomment subpopsim lines below.
# For neutral and fd simulations with SSA demography (Fig S4-5), compile and
# run ml/simssa.ml, placing the timeseries-*.tsv files in out/
# before running this analysis.

## === generate timeseries ===
echo "Generating baby names timeseries..."
cat <<EOM | time parallel
time ./timeseries-ssa.sh
time ./timeseries-france.sh
time ./timeseries-netherlands.sh
#time ./timeseries-subpopsim.sh
time Rscript timeseries-norway.R
EOM
# some of the scripts must be done serially...
cat <<EOM | time parallel
time Rscript timeseries-ssabibcontrol.R 0
time Rscript timeseries-netherlands.R
time Rscript timeseries-france.R
#time Rscript timeseries-subpopsim.R
time Rscript timeseries-samples.R
time Rscript timeseries-ssa.R
EOM

function outflags {
TAG=$1
echo -p out/inf-$TAG.pop_sizes.tsv \
  -b out/inf-$TAG.breaks.tsv \
  -o out/inf-$TAG.params.tsv
}

function resflag {
TAG=$1
echo -P out/inf-$TAG.update_residuals.tsv
}

function bsflag {
TAG=$1
echo -S out/inf-$TAG.bootstrap.tsv
}

function infer {
FLAGS=$1
INFFLAGS=$2
shift 2
TAG=$FLAGS$INFFLAGS
COMMAND="fdsel infer -i out/timeseries-$FLAGS.tsv `outflags $TAG` $@"
echo $COMMAND
$COMMAND
}

function inferB {
FLAGS=$1
INFFLAGS=$2
shift 2
TAG=$FLAGS$INFFLAGS
COMMAND="fdsel infer -i out/timeseries-$FLAGS.tsv \
  `outflags $TAG` `bsflag $TAG` $@"
echo $COMMAND
$COMMAND
}

function inferBR {
FLAGS=$1
INFFLAGS=$2
shift 2
TAG=$FLAGS$INFFLAGS
COMMAND="fdsel infer -i out/timeseries-$FLAGS.tsv \
  `outflags $TAG` `bsflag $TAG` `resflag $TAG` $@"
echo $COMMAND
$COMMAND
}

function inferR {
FLAGS=$1
INFFLAGS=$2
shift 2
TAG=$FLAGS$INFFLAGS
COMMAND="fdsel infer -i out/timeseries-$FLAGS.tsv \
  `outflags $TAG` `resflag $TAG` $@"
echo $COMMAND
$COMMAND
}

export -f infer inferB inferBR inferR outflags resflag bsflag

# === Inferences ===
echo "Inferring from timeseries..."
cat <<EOM | time parallel
inferB ssaCd35dCE 2Bl -B inp/log2.breaks.tsv -c 5 -w CENSORED
inferBR ssaCd35dCE 2Bu -B inp/log2.breaks.tsv -c 5 -C -w CENSORED
infer ssasamples l20 -l 20
inferB ssabib 2Bl -B inp/log2.breaks.tsv -c 5 -w NONBIBCEN
inferB ssabib 2Bu -B inp/log2.breaks.tsv -c 5 -C -w NONBIBCEN
inferB ssanob 2Bl -B inp/log2.breaks.tsv -c 5 -w BIBCEN
inferB ssanob 2Bu -B inp/log2.breaks.tsv -c 5 -C -w BIBCEN
inferB ssaM 2Bl -B inp/log2.breaks.tsv -c 5 -w FCEN
inferB ssaM 2Bu -B inp/log2.breaks.tsv -c 5 -C -w FCEN
inferB ssaF 2Bl -B inp/log2.breaks.tsv -c 5 -w MCEN
inferB ssaF 2Bu -B inp/log2.breaks.tsv -c 5 -C -w MCEN
inferR ssabibctrl000 2Bl -B inp/log2.breaks.tsv -c 5 -w NONBIBCEN
inferR ssabibctrl001 2Bu -B inp/log2.breaks.tsv -c 5 -C -w NONBIBCEN
inferR ssanobctrl000 2Bl -B inp/log2.breaks.tsv -c 5 -w BIBCEN
inferR ssanobctrl001 2Bu -B inp/log2.breaks.tsv -c 5 -C -w BIBCEN
inferB simparlin 2B -B inp/log2.breaks.tsv
inferB simparlinm 2B -B inp/log2.breaks.tsv -w WT
inferB simparlinM 2B -B inp/log2.breaks.tsv -w WT
inferB france 2Bl -B inp/log2.breaks.tsv -c 3 -w CENSORED
inferBR france 2Bu -B inp/log2.breaks.tsv -c 3 -C -w CENSORED
infer francesamples l20 -l 20
infer francesamplescen 2Bl -B inp/log2.breaks.tsv -c 3 -w CENSORED
infer francesamplescen 2Bu -B inp/log2.breaks.tsv -c 3 -C -w CENSORED
inferBR netherlands 2B -B inp/log2.breaks.tsv
infer netherlandssamples 2B -B inp/log2.breaks.tsv
inferBR netherlands10y 2B -B inp/log2.breaks.tsv
inferBR netherlands5y 2B -B inp/log2.breaks.tsv
inferBR netherlands2y 2B -B inp/log2.breaks.tsv
infer netherlands5ysamples 2B -B inp/log2.breaks.tsv
infer netherlands5y l25 -l 25
inferB netherlandsmusim5yp1yd 2B -B inp/log2.breaks.tsv
inferB netherlandsmusim5yp1yd l25B -B inp/inf-netherlands5yl25.breaks.tsv
infer norwayd45 2Bl -B inp/log2.breaks.tsv -c 4 -w CENSORED
inferBR norwayd45 2Bu -B inp/log2.breaks.tsv -c 4 -C -w CENSORED
infer norwayd45samplescen 2Bl -B inp/log2.breaks.tsv -c 4 -w CENSORED
infer norwayd45samplescen 2Bu -B inp/log2.breaks.tsv -c 4 -C -w CENSORED
#infer ssaKCd35dCE 2Bl -B inp/log2.breaks.tsv -c 5 -w CENSORED
#infer ssaKCd35dCE 2Bu -B inp/log2.breaks.tsv -c 5 -C -w CENSORED
#inferB ssaKCd35dCE q20u -q 20 -f 0.00001 -c 5 -C -w CENSORED
#inferB ssaKCd35dCE l25u -l 25 -f 0.00001 -c 5 -C -w CENSORED
#inferB ssaKCd35dCE l2u -l 2 -f 0.00001 -c 5 -C -w CENSORED
#inferB ssaKCd35dCE 2Bfu -B inp/log2.breaks.tsv -f 0.00001 -c 5 -C -w CENSORED
EOM

# Report which bins are bad... bad if the sampling bias is present (>0.008) or
# if there are less than 5 types.
echo "Writing bad bins file..."
R << EOM
options(stringsAsFactors = FALSE)
datasets = data.frame(
  src = c("ssaCd35dCE2Bu", "france2Bu", "netherlands2B", "netherlands5y2B",
    "norwayd452Bu"),
  sampcen = c("ssasamplescen2Bu", "francesamplescen2Bu", 
   "netherlandssamples2B", "netherlands5ysamples2B", "norwayd45samplescen2Bu"),
  #name = c("Alice,F","1,FRANÇOIS","M90154","M90154","Bj\xf8rn,M"),
  name = c("Sarah,F","1,FRANÇOIS","M90154","M90154","Bj\xf8rn,M"),
  wtoffset = c(1, 1, 0, 0, 1))

source("../common.R")

bininfo = data.frame()
nameranges = data.frame()
for (fi in datasets$src) {
  dsinfo = subset(datasets, src == fi)
  sampcenparams = make_s_of_ind(read.delim(sprintf("out/inf-%s.params.tsv",
    dsinfo$sampcen)))
  ress = read.delim(sprintf("out/inf-%s.update_residuals.tsv", fi))
  namesubress = subset(ress, type==dsinfo$name)
  nameminf = min(namesubress$obsf,namesubress$initf)
  namemaxf = max(namesubress$obsf,namesubress$initf)
  nameranges = rbind(nameranges, data.frame(
    src=fi, name=dsinfo$name, min=nameminf, max=namemaxf))
  tyinds = unique(ress[,c("type","ind")])
  tyinds = transform(tyinds,
    ind = ind - dsinfo$wtoffset,
    ntypes = as.numeric(ave(type, ind, FUN=length)))
  tyinds = unique(tyinds[,c("ind","ntypes")])
  sampcenparams$maxind = max(subset(sampcenparams, !is.nan(val))$ind)
  sampcenparams = merge(sampcenparams, tyinds)
  bininfo = rbind(bininfo, transform(sampcenparams, src=fi)) 

  # output fraction of biblical names for ssa
  if(fi == "ssaCd35dCE2Bu") {
    bnames <- read.delim("inp/bible-names.txt",header=FALSE)$V1
    # put names in SSA format indiscriminately w.r.t gender:
    gbnames = c(paste(bnames, "M", sep=","), paste(bnames, "F", sep=","))
    count_biblical = function(df) {
      names = df$type
      biblical = names %in% gbnames
      return(data.frame(
        namesfrac = sum(unique(names) %in% gbnames)/length(unique(names)),
        nameyearsfrac = sum(names %in% gbnames)/length(names))) }
    fracs = rbindlist(by(ress, ress$ind, count_biblical),idcol="ind")
    write.tsv(fracs, "out/ssaCd35dCE2Bu.bibfracs.tsv")

    topnames = subset(ress, 
      type %in% c("Robert,M", "James,M", "Linda,F", "Michael,M", "John,M"))
    write.tsv(topnames, "out/ssaCd35dCE2Bu.topress.tsv")
  } }

bad_bins = subset(bininfo, 
  (abs(val) > 0.008 & ind < maxind) | ntypes < 5)[,c("src","ind","ntypes")]
  
write.tsv(bad_bins, "out/bad_bins.tsv")
write.tsv(nameranges, "out/nameranges.tsv")
EOM

# Write permuted biblical namesfrac
R << EOM
options(stringsAsFactors = FALSE)
source("../common.R")
ressbib = read.delim("out/inf-ssabibctrl0012Bu.update_residuals.tsv")
ressnob = read.delim("out/inf-ssanobctrl0012Bu.update_residuals.tsv")
ressbib\$ntrans = with(ressbib, ave(type,ind,FUN=length))
ressnob\$ntrans = with(ressnob, ave(type,ind,FUN=length))
nobntrans = unique(ressnob[,c("ind","ntrans")])
bibntrans = unique(ressbib[,c("ind","ntrans")])
ntrans = merge(bibntrans,nobntrans,by="ind")
ntrans = transform(ntrans, nameyearsfrac = 
  as.numeric(ntrans.x)/(as.numeric(ntrans.x)+as.numeric(ntrans.y)))
write.tsv(ntrans[,c("ind","nameyearsfrac")],
  "out/ssabibctrl0012bu.bibfracs.tsv")
EOM

# Write unique names and top name frequencies
R << EOM
options(stringsAsFactors = FALSE)
source("../common.R")
stats = function(src, centype) {
  ts = read.delim(sprintf("../names/out/timeseries-%s.tsv",src))
  ts$total = ave(ts$count, ts$type, FUN=sum)
  all_total = sum(ts$count)
  typects = unique(ts[,c("type","total")])
  typects = subset(typects, type != centype)
  typects = typects[order(-typects$total),]
  top100count = sum(typects[seq(1,100),"total"])
  top100pct = top100count / all_total
  nnames = length(unique(typects$type))
  return(data.frame(src=src,top100pct=top100pct,nnames=nnames,total=all_total)) }
write.tsv(rbind(
    stats("ssaCd35dCE","CENSORED"),
    stats("ssad35dCE",""),
    stats("france","CENSORED"),
    stats("netherlands",""),
    stats("netherlands5y",""),
    stats("norwayd45","CENSORED"),
    stats("norwayrawd45","")),
  "out/ts-stats.tsv")
EOM

# # Generate a distribution of ssabibcontrol timeseries
# # VERY SLOW AND CONSUMES HUGE AMOUNTS OF DISK SPACE
# function infer_nr {
# FLAGS=$1
# INFFLAGS=$2
# shift 2
# TAG=$FLAGS$INFFLAGS
# COMMAND="fdsel infer -i out/timeseries-$FLAGS.tsv \
#   -p out/inf-$TAG.pop_sizes.tsv \
#   -b out/inf-$TAG.breaks.tsv \
#   -o out/inf-$TAG.params.tsv $@"
# echo $COMMAND
# $COMMAND
# }
# function resample {
# IND=$1
# Rscript timeseries-ssabibcontrol.R $IND
# infer_nr ssabibctrl$IND 2Bl -B inp/log2.breaks.tsv -c 5 -w NONBIBCEN
# infer_nr ssanobctrl$IND 2Bl -B inp/log2.breaks.tsv -c 5 -w BIBCEN
# }
# export -f infer_nr resample
# time parallel resample {1} :::: <(seq -w 1 500)

