library(ggplot2)
# Requires: common.R

minilog10x = scale_x_log10(
  breaks=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001),
  labels=c(expression(paste(phantom(1^1), 1, phantom(1^1))),
  expression(10^-1),
  expression(10^-2),
  expression(10^-3),
  expression(10^-4),
  expression(10^-5),
  expression(10^-6),
  expression(10^-7),
  expression(10^-8)))

minilog10y = scale_y_log10(
  breaks=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001),
  labels=c(expression(paste(phantom(1^1), 1, phantom(1^1))),
  expression(10^-1),
  expression(10^-2),
  expression(10^-3),
  expression(10^-4),
  expression(10^-5),
  expression(10^-6),
  expression(10^-7),
  expression(10^-8)))

segs_to_line = function(df) {
  start = df[,c("x","y")]
  end = df[,c("xend","yend")]
  names(end) <- c("x","y")
  return(unique(rbind(start, end))) }

# geometries for plotting residuals as a timeseries
aes_r = function(dx,...) aes(y=initf,x=gen,xend=(gen + dx),yend=obsf,...)
gg_ress = function(ress,aes,...) geom_segment(data=ress,aes,...)

# geometry for plotting frequency bin boundaries as horizontal lines
gg_hbreaks = function (breaks,...) 
  geom_hline(data=data.frame(freq=breaks),aes(yintercept=freq),alpha=0.2,...)

# Plot raw data residuals p-value
plot_residual_p_files = function(basefn,ne=NA,wt=FALSE) {
  params = read_params(basefn)
  ress = read_udresiduals(basefn)
  breaks = read_breaks(basefn)
  if (wt) {
    ress = subset(ress, ind != 0) }
  if (is.na(ne)) {
    ne = get_param(params, "ne") }

  # A rough and ready p-value of each residual based on the neutral model
  # ne is arbitrary and sets the scale of p
  twosided = Vectorize(function (p) ifelse(p < 0.5, 2*p, 2*(1-p)))
  nmean = exp(mean(log(unique(ress[,c("gen","fin_ps")])$fin_ps)))
  noverne = nmean/ne
  ress$p = twosided(pbinom(floor(ress$fin_ps/noverne)*ress$obsf,
    floor(ress$fin_ps/noverne), ress$expf))

  ress = ress[order(-ress$p),]
  dx = gen_interval(ress)
  minfreq = get_param(params, "minf")
  
  return(ggplot() +
    gg_hbreaks(breaks) +
    geom_hline(yintercept=minfreq) +
    gg_ress(subset(ress, type!="MUT"),aes_r(dx,color=log10(p)),alpha=0.7) +
    scale_color_viridis_c("log10 pseudo-P",option="inferno")) }

aes_tf = function(...) aes(x=gen,y=freq,group=type,...)
gg_timeseries_gtf = function(gtf,aest,...) geom_line(data=gtf,aest,...)

plot_timeseries_file = function(fn,...) {
  gtf = ts_gtc_to_gtf(read.delim(fn))
  return(ggplot() + gg_timeseries_gtf(gtf,aes_tf(),...)) }

# These expect "sp" columns: ind, val, lci, uci, and midbin
aes_pl = function (yunits, ...) aes(x=midbin, y=yunits(val), ...)
aes_pl_raw = function (...) aes(x=midbin, y=val, ...)
gg_params_line = function(sp,aesp, ...) geom_line(data=sp,aesp, ...)
gg_params_point = function(sp,aesp, ...) geom_point(data=sp,aesp, ...)

aes_pr = function (yunits, ...) aes(x=midbin,
  ymin=yunits(lci), ymax=yunits(uci), ...)
aes_pr_raw = function (...) aes(x=midbin, ymin=lci, ymax=uci, ...)
gg_params_ribbon = function(sp,aesp, ...) geom_ribbon(data=sp,aesp, ...)

# Requires columns lbb, ubb, uci, lci
aes_prect = function (yunits, ...) aes(
  xmin=lbb, xmax=ubb,
  ymin=yunits(lci), ymax=yunits(uci), ...)
aes_prect_raw = function (...) aes(xmin=lbb, xmax=ubb, ymin=lci, ymax=uci, ...)
gg_params_rect = function(sp,aesp, ...) geom_rect(data=sp,aesp,...)
  
# Plot inference for an arbitrary number of series
plot_inference = function(basefn,...) {
  params = read_params(basefn)
  pops = read_pops(basefn) # smallest output file containing all gens.
  breaks = read_breaks(basefn)
  midbins = make_params_logmidbins(params, breaks)
  sp = merge(make_s_of_ind(params),midbins)
  dx = gen_interval(pops)
  ypt = growth_units(dx)

  pl = ggplot() +
    geom_hline(yintercept=0,linetype="11",color="#000000",alpha=0.5) +
    gg_params_ribbon(sp,aes_pr(ypt),...) +
    gg_params_line(sp,aes_pl(ypt),...) +
    gg_params_point(sp,aes_pl(ypt),...) +
    ylab("% per timeseries unit") +
    xlab("Frequency range midpoint")
    
  if(exists_bootstrap(basefn)) {
    bs = read_bootstrap(basefn)
    bssp = merge(make_bootstrap_cis(bs),midbins)
    
    pl = pl +
      gg_params_ribbon(bssp,aes_pr(ypt),linetype="11",fill=NA,
        color="black") }

  return(pl) }
plot_inference_files = plot_inference # weird backwards compatibility

plot_inf_comparison_dx = function(basefns,...) {
  nplots=length(basefns)

  # Assume comparison is happening over comparable units
  get_dx = Vectorize(function (bfn) { return(gen_interval(read_pops(bfn))) })
  dxs = get_dx(basefns)

  # plot accumulator, initially mostly blank
  pl = ggplot() +
    geom_hline(yintercept=0,linetype="11",color="#000000",alpha=0.5)

  sps = data.frame()
  for (ix in 1:nplots) {
    basefn = basefns[ix]
    dx = dxs[ix]
    units = growth_units(dx)
    params = read_params(basefn)
    breaks = read_breaks(basefn)
    midbins = make_params_logmidbins(params, breaks)
    this_sps = cbind(merge(make_s_of_ind(params),midbins),fn=basefn)
    sps = rbind(sps, 
      transform(this_sps, val=units(val), lci=units(lci), uci=units(uci)))
 
    # Add bootstrap rectangles to the plot if the bootstrap exists
    if(exists_bootstrap(basefn)) {
      bs = read_bootstrap(basefn)
      bssp = merge(make_bootstrap_cis(bs),midbins)
      bssp = transform(bssp, lcu = units(lci), uci = units(uci))
      pl = pl + 
        gg_params_rect(bssp,aes_prect_raw(group=fn),
          fill="black",alpha=0.5/nplots, color=NA) }}

  # Add the s(p)s to the plot
  pl = pl +
    gg_params_ribbon(sps,aes_pr_raw(fill=fn),alpha=1/nplots,...) +
    gg_params_line(sps,aes_pl_raw(color=fn),alpha=3/nplots,...) +
    ylab("% per timeseries unit") +
    xlab("Frequency range midpoint")

  return(pl) }

plot_inf_comparison = function(basefns,...) {
  nplots=length(basefns)

  # Assume comparison is happening over comparable units
  dx = gen_interval(read_pops(basefns[1]))
  ypt = growth_units(dx)

  # plot accumulator, initially mostly blank
  pl = ggplot() +
    geom_hline(yintercept=0,linetype="11",color="#000000",alpha=0.5)

  sps = data.frame()
  for (ix in 1:nplots) {
    basefn = basefns[ix]
    params = read_params(basefn)
    breaks = read_breaks(basefn)
    midbins = make_params_logmidbins(params, breaks)
    sps = rbind(sps, cbind(merge(make_s_of_ind(params),midbins),fn=basefn))
 
    # Add bootstrap rectangles to the plot if the bootstrap exists
    if(exists_bootstrap(basefn)) {
      bs = read_bootstrap(basefn)
      if(length(get_param(params, "swt")) > 0) {
        bs$ind = bs$ind - 1
        bs$ind[bs$ind == -1] = NA }
      maxind = max(subset(params, !is.nan(val) & !is.na(ind))$ind)
      bs$ind[bs$ind > maxind] = NA 
      bssp = merge(make_bootstrap_cis(bs),midbins)
      pl = pl + 
        gg_params_rect(bssp,aes_prect(ypt,group=fn),
          fill="black",alpha=0.5/nplots, color=NA) }}

  # Add the s(p)s to the plot
  pl = pl +
    gg_params_ribbon(sps,aes_pr(ypt,fill=fn),alpha=0.8/nplots,...) +
    gg_params_line(sps,aes_pl(ypt,color=fn),alpha=5/nplots,...) +
    gg_params_point(sps,aes_pl(ypt,color=fn),alpha=1/nplots,...) +
    ylab("% per timeseries unit") +
    xlab("Frequency range midpoint")

  return(pl) }

library(deldir)
library(sqldf)

aes_v = function (...) aes(x=x,y=y,group=vind,...)
gg_voronoi = function (hulls,aesv,...) geom_polygon(data=hulls,aesv,...)

plot_nb_grid_heatmap = function(mps,aesv,...) {
  mps$vind = 1:nrow(mps)
  vnoi = deldir(log10(mps$mu), log10(mps$s))
  deldir=vnoi$dirsgs
  allcoords = unique(rbind(
    sqldf("select mps.*, x1 as x, y1 as y from mps,deldir where vind = ind1"),
    sqldf("select mps.*, x2 as x, y2 as y from mps,deldir where vind = ind1"),
    sqldf("select mps.*, x1 as x, y1 as y from mps,deldir where vind = ind2"),
    sqldf("select mps.*, x2 as x, y2 as y from mps,deldir where vind = ind2")))
  
  allhulls = data.frame()
  for (thisvind in unique(allcoords$vind)) {
    thiscoords = subset(allcoords, vind == thisvind)
    allhulls = rbind(allhulls, thiscoords[with(thiscoords,chull(x,y)),]) }

  return(ggplot() +
    gg_voronoi(allhulls,aesv,...) +
    scale_fill_viridis_c(option="inferno")) }
