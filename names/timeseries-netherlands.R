options(stringsAsFactors = FALSE)
source("../common.R")

# Censorship
tscen = read.delim("out/timeseries-netherlands.tsv")
tscen[tscen$count < 3,]$type = "CENSORED"
tscen = transform(tscen, gentype = paste(gen, type))
tscen = transform(tscen, count = ave(count, gentype, FUN=sum))
tscen = unique(tscen[,c("gen","type","count")])
write.tsv(tscen, "out/timeseries-netherlands3cen.tsv")

# Resampling
ts5xsim = read.delim("out/timeseries-netherlands5xsim.tsv")
pops = by(ts5xsim, ts5xsim$gen, function (df) {
  ps = sum(df$count)
  freqs = df$count/ps
  df$count = rmultinom(1, ps/5, freqs)
  return(subset(df, count > 0)) })
tssimsamp = rbindlist(pops)
write.tsv(subset(tssimsamp, as.integer(gen) > 0),
  "out/timeseries-netherlandssimsamp.tsv")
