source("../common.R")

# === no immigration timeseries ===

for (fn in c("simpar", "simlin", "simparlin", "simlinoff")) {
  ts = read.delim(sprintf("out/timeseries-%s.tsv", fn))
  write.tsv(make_type_dist(ts), sprintf("out/%s.typefreqs.tsv", fn))
  
  mintypes = unique(subset(ts, type %% 3 == 0)$type)
  majtypes = unique(subset(ts, type %% 3 != 0)$type)
  
  tsmin = merge_wildtype(ts, "WT", majtypes)
  tsmaj = merge_wildtype(ts, "WT", mintypes)
  write.tsv(tsmin, sprintf("out/timeseries-%sm.tsv", fn))
  write.tsv(tsmaj, sprintf("out/timeseries-%sM.tsv", fn))
}
