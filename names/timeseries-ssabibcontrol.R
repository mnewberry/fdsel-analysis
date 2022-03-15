#!/usr/bin/Rscript
source("../common.R")

ind = as.integer(commandArgs(trailingOnly=TRUE)[1])
#set.seed(ind)

ts <- read.delim("out/timeseries-ssaCd35dCE.tsv")
l2breaks = read_breaks("inp/log2")

# === Compute biblical and non-biblical name treatments ===
bnames <- read.delim("inp/bible-names.txt",header=FALSE)$V1
gbnames = c(paste(bnames, "M", sep=","), paste(bnames, "F", sep=","))

ts$biblical = ts$type %in% gbnames

# Scramble biblical and nonbiblical names roughly preserving the frequency
# distribution of scrambled "biblical" names.
bib_control = function (ts) {
  ts = transform(ts, typesum = ave(count, type, FUN=sum))
  dist = unique(ts[,c("type","typesum","biblical")])
  dist = transform(dist, freq = typesum / sum(typesum))
  distfreqs = dist$freq[dist$typesum > 10]
#  dist = transform(dist,
#    bin = cut(freq, breaks=quantile(distfreqs, probs=seq(0,1,0.05))))
  dist = transform(dist, bin = cut(freq, breaks=c(0,l2breaks,1)))
  dist = rbindlist(by(dist, dist$bin, FUN=function (df) {
    df$biblical = sample(df$biblical)
    return(df) }))
  dist = transform(dist,
    bibcat = ifelse(biblical, type, "NONBIBLICAL"),
    nobcat = ifelse(biblical, "BIBLICAL", type))
  dist = unique(dist[,c("type","bibcat","nobcat")])
  ts = merge(ts[,c("gen","type","count")], dist)
  ts = transform(ts,
    bibgc = paste(gen, bibcat),
    nobgc = paste(gen, nobcat))
  ts = transform(ts,
    bibcount = ave(count, bibgc, FUN=sum),
    nobcount = ave(count, nobgc, FUN=sum),
    ps = ave(count, gen, FUN=sum))
  return(ts) }

bib_slice = function (bibts) {
  bibts = unique(bibts[,c("gen","bibcat","bibcount")])
  names(bibts) <- c("gen","type","count")
  bibts$type[bibts$type == "NONBIBLICAL"] = "NONBIBCEN"
  return(bibts) }

nob_slice = function (bibts) {
  bibts = unique(bibts[,c("gen","nobcat","nobcount")])
  names(bibts) <- c("gen","type","count")
  bibcen = bibts$type %in% c("BIBLICAL", "CENSORED")
  bibts$type[bibcen] = "BIBCEN"
  bibts$count[bibcen] = ave(bibts$count[bibcen],bibts$gen[bibcen],FUN=sum)
  bibts = unique(bibts)
  return(bibts) }

bibtscontrol = bib_control(ts)
bibc = bib_slice(bibtscontrol)
nobc = nob_slice(bibtscontrol)
write.tsv(bibc, sprintf("out/timeseries-ssabibctrl%03d.tsv", ind))
write.tsv(nobc, sprintf("out/timeseries-ssanobctrl%03d.tsv", ind))

scrambledmale = function(ts) {
  ts = transform(ts, typesum = ave(count, type, FUN=sum),
    male=substr(type,nchar(type),nchar(type)) == "M")
  dist = unique(ts[,c("type","typesum","male")])
  dist = transform(dist, freq = typesum / sum(typesum))
  distfreqs = dist$freq[dist$typesum > 10]
  dist = transform(dist, bin = cut(freq, breaks=c(0,l2breaks,1)))
  dist = rbindlist(by(dist, dist$bin, FUN=function (df) {
      df$male = sample(df$male)
      return(df) 
    }))
  return(unique(subset(dist, male)$type)) }

ssanames = unique(ts$type)
namesM = scrambledmale(ts)
tsM = merge_wildtype(ts, "FCEN", setdiff(ssanames, namesM))
tsF = merge_wildtype(ts, "MCEN", c("CENDORED",namesM))

write.tsv(tsM, sprintf("out/timeseries-ssaMctrl%03d.tsv", ind))
write.tsv(tsF, sprintf("out/timeseries-ssaFctrl%03d.tsv", ind))
