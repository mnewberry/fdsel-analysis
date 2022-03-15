# Common ggplot extensions for plotting timeseries
source("../common.R")
source("../common_plot.R")

mmperpt = 1/72 * 25.4
pt = mmperpt

blank_theme = theme(
  text = element_text(size=7),
  legend.title = element_text(size=7),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.ticks.length = unit(0,"pt"),
  axis.text = element_text(color="#000000"),
  legend.key=element_blank(),
  plot.margin=unit(c(0,0,0,0),c("cm","cm","cm","cm")),
  panel.margin=unit(0,"cm"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.spacing.x=unit(1.25,"mm"),
  panel.spacing.y=unit(1.25,"mm"),
  strip.background=element_rect(color=NA,fill=NA),
  strip.text.x=element_text(margin=margin(t=unit(3,"mm"),r=0,b=0,l=0)),
  panel.background = element_rect(
    fill="transparent",
    colour = NA),
  plot.background = element_rect(
    fill="transparent",
    colour = NA))

compact_legend = theme(
    legend.margin=margin(t=0,r=0,b=0.5,l=0,unit="mm"),
    legend.spacing=unit(1,"pt"),
    legend.key.size=unit(1,"mm"),
    legend.key.width=unit(3,"mm"),
    legend.box.margin=margin(t=0,r=0,b=0,l=0,unit="pt"),
    legend.background=element_blank(),
    legend.box.background=element_blank())

xaxis = theme(
  text = element_text(size=8),
  axis.ticks.length = unit(2,"pt"),
  axis.ticks = element_line(size=0.75*pt,color="#000000"),
  axis.line = element_line(size=0.75*pt,color="#000000"),
  axis.line.x = element_line(size=0.75*pt,color="#000000"))


# Timeseries plot
ts = read.delim("../names/out/timeseries-ssad35dCE.tsv")
ts = annotate_frequencies(ts)
ts = rbindlist(by(ts, ts$gen, function(df) {
  freqs = sort(df$freq,decreasing=TRUE)
  df$avg50 = mean(freqs[seq(1,50)])
  return(df) }))
topnth = transform(unique(ts[,c("gen","avg50")]), src="ssad35dCE")
ggplot(ts) + minilog10y + 
  geom_line(aes(x=gen,y=freq,group=type),alpha=0.05,size=0.75*pt) + 
  geom_line(data=topnth,aes(x=gen,y=avg50,color="red"),size=0.75*pt) +
  blank_theme + xaxis + labs(y="Frequency", x="Year")

# Frequency distribution of types
ssatypedist = read.delim("out/ssad35dCE.typefreqs.tsv")
minfreq = 2.292212e-06
ssatypedist = subset(ssatypedist, freq > minfreq)
ssabibfrac = nrow(subset(ssatypedist, biblical))/nrow(ssatypedist)
ssaecdf = ecdf(ssatypedist$freq)
bibecdf = ecdf(subset(ssatypedist, biblical)$freq)
bibfraccdf = with(environment(bibecdf), 
  make_discrete_cdf(x, (1 - y)*ssabibfrac/(1 - ssaecdf(x))))
bibfraccdf = bibfraccdf[2:nrow(bibfraccdf),]
ssaecdfs = rbind(
  transform(make_ecdf_df(ssatypedist$freq), src="All names"),
  transform(make_ecdf_df(subset(ssatypedist, biblical)$freq),
    y = y * ssabibfrac, src="Biblical names"),
  transform(bibfraccdf, src="Fraction biblical of names >f"))
ggplot(ssaecdfs, aes(x=x,y=y,color=src)) + 
  geom_line() + coord_fixed() + minilog10x + minilog10y +
  labs(color="", x="Name frequency, f",
    y="Fraction greater than f among all names") +
  ggtitle("Empirical tail distribution of name frequencies, SSA")

ssatypedist = read.delim("out/ssad35dCE.typefreqs.tsv")
minfreq = 2.292212e-06
ssatypedist = subset(ssatypedist, freq > minfreq)
ssatypedist$male = with(ssatypedist,
  substr(type,nchar(type),nchar(type)) == "M")
ssaMfrac = nrow(subset(ssatypedist, male))/nrow(ssatypedist)
ssaecdf = ecdf(ssatypedist$freq)
ssaMecdf = ecdf(subset(ssatypedist, male)$freq)
ssaMfraccdf = with(environment(ssaMecdf), 
  make_discrete_cdf(x, (1 - y)*ssaMfrac/(1 - ssaecdf(x))))
ssaMfraccdf = ssaMfraccdf[2:nrow(ssaMfraccdf),]
ssaecdfs = rbind(
  transform(make_ecdf_df(ssatypedist$freq), src="All names"),
  transform(make_ecdf_df(subset(ssatypedist, male)$freq),
    y = y * ssaMfrac, src="Male names"),
  transform(ssaMfraccdf, src="Fraction Male of names >f"))
ggplot(ssaecdfs, aes(x=x,y=y,color=src)) + 
  geom_line() + coord_fixed() + minilog10x + minilog10y +
  labs(color="", x="Name frequency, f",
    y="Fraction greater than f among all names") +
  ggtitle("Empirical tail distribution of name frequencies, SSA")

# Combined Fig 2c and Fig S11b
ecdfplot = function () {
  src_colors = data.frame(
    src = c("ssaCd35dCE2Bu","france2Bu",
      "netherlands2B","netherlands2y2B","netherlands5y2B",
      "netherlandsneusim5yp1ydM","netherlandsneusim5yp1yd",
      "netherlandsmusim5yp1ydM","netherlandsmusim5yp1yd",
      "norwayd452Bu"),
    color=c("#3d219b","#4981c3ff",
      "#e69f00","#e68900","#e17400",
      "#e14f00","#e14f80",
      "#e14f40","#e14fb0",
      "#19ab68"))
  zerobins = data.frame(bound = c(8e-4,3.2e-3))

  mkecdf = function (tsname,ds) {
    typedist = read.delim(sprintf("../names/out/%s.typefreqs.tsv",tsname))
    return(cbind(make_ecdf_df(typedist$freq),src=ds)) }

  ecdfs = rbind(
    mkecdf("ssad35dCE","ssaCd35dCE2Bu"),
    mkecdf("france","france2Bu"),
    mkecdf("netherlandsmusim5yp1ydM","netherlandsmusim5yp1ydM"),
    mkecdf("netherlandsmusim5yp1yd","netherlandsmusim5yp1yd"),
    mkecdf("netherlandsneusim5yp1ydM","netherlandsneusim5yp1ydM"),
    mkecdf("netherlandsneusim5yp1yd","netherlandsneusim5yp1yd"),
    mkecdf("netherlands","netherlands5y2B"),
    mkecdf("norwayd45","norwayd452Bu"))

  ecdfs = merge(ecdfs, src_colors)

  pl = ggplot(ecdfs, aes(x=x,y=y,color=color)) + 
    coord_fixed(expand=FALSE) +
    geom_vline(data=zerobins,aes(xintercept=bound),color="grey") +
    geom_line() +
    scale_discrete_identity(aesthetics = c("colour", "fill")) +
    labs(x="Frequency, f", y="Fraction >f",color="Timeseries") +
    minilog10x + minilog10y
  return(pl) }
print(ecdfplot())

write_corrected_params = function (paramstag, correctiontag, outputtag) {
  system(sprintf(
    "for FI in out/inf-%s.* ; do cp $FI `echo $FI | sed 's/%s/%s/'` ; done",
    paramstag, paramstag, outputtag))
  p1 = read_params(sprintf("out/inf-%s",paramstag))
  p2 = read_params(sprintf("out/inf-%s",correctiontag))
  inds = p1$ind %in% seq(0,22)
  p1$val[inds] = p1$val[inds] - p2$val[inds] + get_param(p2, "srep")
  p1$var[inds] = p1$var[inds] + p2$var[inds]
  p1$lci[inds] = p1$val[inds] - 1.96*sqrt(p1$var[inds])
  p1$uci[inds] = p1$val[inds] + 1.96*sqrt(p1$var[inds])
  write.tsv(p1, sprintf("out/inf-%s.params.tsv",outputtag)) }

# === SSA ===
# Fig 2a censorship controls, extended plot range
plot_inf_comparison(c("out/inf-ssaCd35dCE2Bl","out/inf-ssaCd35dCE2Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("SSA, bounds, dpre1935 dCodingErrors, Hahn bins")

plot_inf_comparison(c(
    "out/inf-ssaCd35dCE2Bl", "out/inf-ssaCd35dCE2Bu",
    "out/inf-ssabib2Bl", "out/inf-ssabib2Bu",
    "out/inf-ssanob2Bl", "out/inf-ssanob2Bu",
    "out/inf-ssabibctrl0002Bl", "out/inf-ssabibctrl0002Bu",
    "out/inf-ssanobctrl0002Bl", "out/inf-ssanobctrl0002Bu")) + 
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("SSA by biblical, non-biblical, all, and randomly permuted controls")

# === Gendered
plot_inf_comparison(c(
    "out/inf-ssaF2Bl", "out/inf-ssaF2Bu",
    "out/inf-ssaM2Bl", "out/inf-ssaM2Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Male versus female names")

# === Gendered data control
plot_inf_comparison(
  c(sprintf("out/inf-ssaFctrl%03d2Bl", seq(1,500)),
    sprintf("out/inf-ssaMctrl%03d2Bl", seq(1,500)))) +
  theme(legend.position="none") +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Resampled ssa randomly permuted gender controls")

# === Bible data control
plot_inf_comparison(
  c(sprintf("out/inf-ssabibctrl%03d2Bl", seq(1,500)),
    sprintf("out/inf-ssanobctrl%03d2Bl", seq(1,500)))) +
  theme(legend.position="none") +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Resampled ssa randomly permuted bible controls")

# === Bible model control ===
ps = 20000
mu = 0.0003
expunit = function (x) { return((exp(x) + mu - 1)*100) }
truess = rbind(
  transform(data.frame(x=seq(1:(ps/5))/ps),
    y=0.01 * (-log10(x) - 4) * (log10(x) + 2),
    src="minority"),
  transform(data.frame(x=seq(1:(ps/5))/ps),
    y=0.012 - 0.006 * (log10(x) + 4),
    src="majority"))

# # For use with ml/parlin.ml
# 
# plot_inf_comparison(
#     c("out/inf-simparlin2B", "out/inf-simparlinm2B", "out/inf-simparlinM2B")) +
#   geom_line(data=truess,aes(x=x,y=expunit(y),linetype=src)) +
#   scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
#   labs(linetype="Absolute fitness", color="Inference", fill="Inference") + 
#   ggtitle("Recovering subpopulation fitness, minority type parabolic")

# === France ===
# Fig 2a censorship controls, extended plot range
plot_inf_comparison(c("out/inf-france2Bl","out/inf-france2Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("France, log2 bins")

# France sampling error controls
write_corrected_params("france2Bl","francesamplescen2Bl","francecorrected2Bl")
write_corrected_params("france2Bu","francesamplescen2Bu","francecorrected2Bu")
plot_inf_comparison(c(
    "out/inf-france2Bl", "out/inf-france2Bu",
    "out/inf-francesamplescen2Bl", "out/inf-francesamplescen2Bu",
    "out/inf-francesamplesl20",
    "out/inf-francecorrected2Bl", "out/inf-francecorrected2Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("France controls, true timeseries vs resampled constant freqs")

# === Netherlands treatments ===
plot_inf_comparison_dx(c("out/inf-netherlands2B","out/inf-netherlands5y2B")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Netherlands, log2 bins")
write_corrected_params(
  "netherlands2B","netherlandssamples2B","netherlandscorrected2B")
write_corrected_params(
  "netherlands5y2B","netherlands5ysamples2B","netherlands5ycorrected2B")
plot_inf_comparison_dx(c(
  "out/inf-netherlands2B",
  "out/inf-netherlandssamples2B",
  "out/inf-netherlandscorrected2B",
  "out/inf-netherlands5yl25",
  "out/inf-netherlandsmusim5yp1yd2B",
  "out/inf-netherlands5y2B",
  "out/inf-netherlands5ysamples2B",
  "out/inf-netherlands5ycorrected2B"))+
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Netherlands controls, true timeseries vs resampled constant freqs")

# === Norway treatments ===
plot_inf_comparison(c("out/inf-norwayd452Bl","out/inf-norwayd452Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Norway, log2 bins")

write_corrected_params("norwayd452Bl","norwayd45samplescen2Bl","norwayd45corrected2Bl")
write_corrected_params("norwayd452Bu","norwayd45samplescen2Bu","norwayd45corrected2Bu")
plot_inf_comparison(c(
    "out/inf-norwayd452Bl", "out/inf-norwayd452Bu",
    "out/inf-norwayd45samplescen2Bl", "out/inf-norwayd45samplescen2Bu",
    #"out/inf-norwayd45samplesl20",
    "out/inf-norwayd45corrected2Bl", "out/inf-norwayd45corrected2Bu")) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Norway controls, true timeseries vs resampled constant freqs")

# === cross-country inferences ===
plot_inf_comparison_dx(c(
    "out/inf-france2Bl","out/inf-france2Bu",
    "out/inf-ssaCd35dCE2Bl","out/inf-ssaCd35dCE2Bu",
    "out/inf-netherlands5y2B",
    "out/inf-norwayd452Bl","out/inf-norwayd452Bu"
    )) +
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Country comparisons")
