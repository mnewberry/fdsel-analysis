options(stringsAsFactors = FALSE)
args = commandArgs(TRUE)
inglob = args[1]
tsout = args[2]
wpout = args[3]
ss= as.numeric(args[4])
mu= as.numeric(args[5])
nbins = as.numeric(args[6])

source("../common.R")

files = system(sprintf("ls %s",inglob), intern=TRUE)

last.ts = read.delim(files[1])
minp = 0.5
write.tsv(last.ts,tsout)
for(runfn in files[seq(2,length(files))]) {
  ts = read.delim(runfn)
  minp = min(0.5,ts_gtc_to_gtf(ts)$freq)

  lastgen = subset(last.ts, gen == max(gen))
  firstgen = subset(ts, gen == min(gen))
  offsetdf = cbind(lastgen, firstgen)
  names(offsetdf) <- c("g1","t1","c1","g2","t2","c2")
  offset = unique(offsetdf$t1 - offsetdf$t2)
  stopifnot(length(offset) == 1) # Make sure tss are stitchable

  ts$type = ts$type + offset
  ts = subset(ts, gen > 0)
  ts$gen = ts$gen + max(last.ts$gen)

  write.tsv(ts, tsout, append=TRUE, col.names=FALSE)

  last.ts = ts
}
fullts = read.delim(tsout)

wp = make_novelty_wp(fullts,logit_breaks(1.5*minp,nbins),ss,mu)
write.tsv(wp, wpout)
