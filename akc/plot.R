# Common ggplot extensions for plotting timeseries
source("../common.R")
source("../common_plot.R")

# === Produce main figure results ===
mainfig = function () {
  series = "out/inf-diyoq10"
  params = read_params(series)
  boots = read_bootstrap(series)

  # In order for coarsely-binned s(p) to track high resolution s_ave(p), we
  # need to use midbins that represent the data rather than just the midpoint
  # of the frequency interval
  ressbin = annotate_residual_bins(read_residuals(series),read_breaks(series))
  midbins = make_mean_freq_midbins(ressbin)

  sp = merge(make_s_of_ind(params), midbins)
  bssp = merge(make_bootstrap_cis(boots), midbins)
  ypt = growth_units(1)

  avesps = read.delim("out/wplog_nbcpld-long-s0.00021-mu0.00004.tsv")
  avesp = merge(make_avg_s_cis(avesps),midbins)

  highresavespstitched = transform(read.delim(
      "out/wp_nbcpld-long-s0.00021-mu0.00004-stitched.tsv"),
    val = avg_s, midbin=meanfreq)
  highresavesp = transform(read.delim(
      "out/wp_nbcpld-long-s0.00021-mu0.00004.tsv"),
    val = avg_s, midbin=meanfreq)

  return(ggplot() +
    # s = 0 line
    geom_hline(yintercept=0,linetype="11",color="#000000",alpha=0.5) +

    # Green ribbon from empirical middle 95% quantiles of model s(p)
    gg_params_ribbon(avesp,aes_pr(ypt),fill="#19ab68",alpha=0.4) +

    # Blue and dotted ribbons for data inference
    gg_params_ribbon(sp,aes_pr(ypt),fill="#004492",alpha=0.4) +
    gg_params_ribbon(bssp,aes_pr(ypt),fill=NA,color="black",linetype="11") +

    # mean model s(p) over 80-gen chunks
    gg_params_line(avesp,aes_pl(ypt),color="#19ab68") +

    # High-res s_ave(p), stitched version
    gg_params_line(highresavespstitched,aes_pl(ypt),color="black",alpha=0.5) +

    # High-res s_ave(p), continuous version, independent replicate of above
    gg_params_line(highresavesp,aes_pl(ypt),color="black",alpha=0.5) +

    # data inference line
    gg_params_line(sp,aes_pl(ypt),color="#004492") +

    # labels
    scale_x_log10("Frequency"))
}

# main figure 
mainfig() + coord_cartesian(xlim=c(0.00001,1),ylim=c(-2,20)) +
  ylab("Growth rate (%/year)") +
  ggtitle("Figure with extra replicate of NB long run mean") 

# Full scope
mainfig() + ylab("Growth rate (%/year)") +
  ggtitle("Figure with no plot limit")

# Representative timeseries results
plot_timeseries_file(
    "out/timeseries_nbcpld-long-s0.00021-mu0.00004-r0003.tsv") +
  scale_y_log10("Frequency") + xlab("Generation") + 
  ggtitle("Representative example of novelty bias simulation.")

# === Plot the residuals and inference before and after data cleanup ===
# the parameter 5,000 is an N_e value that is a coefficient of the P-value,
# hence the P-values before and after cleaning are in comparable units,
# though the color scale differs between the plots
plot_residual_p_files("out/inf",5000) + 
  scale_y_log10("Breed Registration Frequency") +
  xlab("Year") + ggtitle("Raw data timeseries pseudo-P-value") 

plot_residual_p_files("out/inf-diyoq10",5000) +
  scale_y_log10("Breed Registration Frequency") +
  xlab("Year") + ggtitle("Timeseries pseudo-P-value after outlier removal") 

# === Plot inference before and after cleanup ===
# Plot interpreation is exactly as paper figure (inference +- model (blue) and
# bootstrap (dotted) CIs)
# "out/inf" indicates the raw inference files out/inf.*
plot_inference_files("out/inf",alpha=0.8,color="#004492",fill="#004492") +
  scale_x_log10("Frequency") + ylab("Intrinsic growth rate (% per year)") +
  ggtitle("Inference on raw timeseries")

# diyoq10 indicates removal of initial years, outliers and 10 quantile bins
plot_inference_files("out/inf-diyoq10",
    alpha=0.8,color="#004492",fill="#004492") +
  scale_x_log10("Frequency") + ylab("Intrinsic growth rate (% per year)") +
  ggtitle("Inference after data cleaning")

# === Sensitivity of diyoq10 series to binning ===
binnings = paste("out/inf-diyoq",sprintf("%02d",seq(6,12)),sep="")
plot_inf_comparison(binnings) + 
  scale_x_log10("Frequency") + ylab("Growth rate (%/year)") +
  ggtitle("Inferences w/ 6 to 12 quantile bins, bin bounds and bootstrap CIs (rectangles).")

# === Novelty bias model parameter fits ===
wps = subset(read.delim("out/wplog_nbcpld.tsv"), run > 3)
avgwps = filewise_make_avg_s_cis(wps)

sp = make_s_of_ind(read_params("out/inf-diyoq10"))
# hack to compute avg_ntypes from timeseries within bins
tsnt = make_novelty_wp(
  transform(read.delim("out/timeseries-diyo.tsv"),
    type=as.numeric(factor(type))),
  read_breaks("out/inf-diyoq10"),0,0)[,c("ind","avg_ntypes")]
sp = merge(sp,tsnt)

fileavgwps = unique(annotate_errs_avgwp_sp(avgwps,sp)[,
  c("file","smeanerr2","ntypesmeanerr2")])
fileavgwps = annotate_s_mu_of_file(fileavgwps)

sim_params = geom_point(data=data.frame(x=log10(0.00004),y=log10(0.00021)),
  aes(x=x,y=y),color="green")

plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(smeanerr2))) + sim_params +
  coord_fixed() +
  ggtitle(expression(paste(
    "Census N, low drift: Chosen ", s, ", ", mu, " in heatmap of error in ", s(p)))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(err^2)/n)))

 
plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(ntypesmeanerr2)))+sim_params+
  coord_fixed() +
  ggtitle(expression(paste(
    "Census N, low drift: Chosen ", s, ", ", mu, " in error in type frequency distribution"))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(log(n[sim]/n[AKC])^2)/n)))

# === Novelty bias model parameter fits ===
wps = subset(read.delim("out/wplog_nbhpld.tsv"), run > 3)
avgwps = filewise_make_avg_s_cis(wps)

sp = make_s_of_ind(read_params("out/inf-diyoq10"))
# hack to compute avg_ntypes from timeseries within bins
tsnt = make_novelty_wp(
  transform(read.delim("out/timeseries-diyo.tsv"),
    type=as.numeric(factor(type))),
  read_breaks("out/inf-diyoq10"),0,0)[,c("ind","avg_ntypes")]
sp = merge(sp,tsnt)

fileavgwps = unique(annotate_errs_avgwp_sp(avgwps,sp)[,
  c("file","smeanerr2","ntypesmeanerr2")])
fileavgwps = annotate_s_mu_of_file(fileavgwps)

sim_params = geom_point(data=data.frame(x=log10(0.00004),y=log10(0.00021)),
  aes(x=x,y=y),color="green")

plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(smeanerr2))) + sim_params +
  coord_fixed() +
  ggtitle(expression(paste(
    "High N, low drift: Chosen ", s, ", ", mu, " in heatmap of error in ", s(p)))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(err^2)/n)))

 
plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(ntypesmeanerr2)))+sim_params+
  coord_fixed() +
  ggtitle(expression(paste(
    "High N, low drift: Chosen ", s, ", ", mu, " in error in type frequency distribution"))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(log(n[sim]/n[AKC])^2)/n)))

# === Novelty bias model parameter fits ===
wps = subset(read.delim("out/wplog_nbcphd.tsv"), run > 3)
avgwps = filewise_make_avg_s_cis(wps)

sp = make_s_of_ind(read_params("out/inf-diyoq10"))
# hack to compute avg_ntypes from timeseries within bins
tsnt = make_novelty_wp(
  transform(read.delim("out/timeseries-diyo.tsv"),
    type=as.numeric(factor(type))),
  read_breaks("out/inf-diyoq10"),0,0)[,c("ind","avg_ntypes")]
sp = merge(sp,tsnt)

fileavgwps = unique(annotate_errs_avgwp_sp(avgwps,sp)[,
  c("file","smeanerr2","ntypesmeanerr2")])
fileavgwps = annotate_s_mu_of_file(fileavgwps)

sim_params = geom_point(data=data.frame(x=log10(0.00004),y=log10(0.00021)),
  aes(x=x,y=y),color="green")

plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(smeanerr2))) + sim_params +
  coord_fixed() +
  ggtitle(expression(paste(
    "Census N, high drift: Chosen ", s, ", ", mu, " in heatmap of error in ", s(p)))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(err^2)/n)))

 
plot_nb_grid_heatmap(fileavgwps,aes_v(fill=-log10(ntypesmeanerr2)))+sim_params+
  coord_fixed() +
  ggtitle(expression(paste(
    "Census N, high drift: Chosen ", s, ", ", mu, " in error in type frequency distribution"))) +
  xlab(expression(mu)) + ylab("s") + 
  labs(fill=expression(-log[10](sum(log(n[sim]/n[AKC])^2)/n)))
