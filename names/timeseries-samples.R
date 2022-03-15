source("../common.R")

ssa = make_samples_control(read.delim("out/timeseries-ssaCd35dCE.tsv"))
write.tsv(ssa, "out/timeseries-ssasamples.tsv")
write.tsv(subset(ssa, count > 4), "out/timeseries-ssasamplescen.tsv")

fra = make_samples_control(read.delim("out/timeseries-france.tsv"))
write.tsv(fra, "out/timeseries-francesamples.tsv")
write.tsv(subset(fra, count > 2), "out/timeseries-francesamplescen.tsv")

net = make_samples_control(read.delim("out/timeseries-netherlands.tsv"))
write.tsv(net, "out/timeseries-netherlandssamples.tsv")
net5 = make_samples_control(read.delim("out/timeseries-netherlands5y.tsv"))
write.tsv(net5, "out/timeseries-netherlands5ysamples.tsv")

nor = make_samples_control(read.delim("out/timeseries-norwayd45.tsv"))
write.tsv(nor, "out/timeseries-norwayd45samples.tsv")
write.tsv(subset(nor, count > 3), "out/timeseries-norwayd45samplescen.tsv")
