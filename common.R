library("data.table") # for rbindlist
options(stringsAsFactors = FALSE)

# Convert dataframe from years as columns to year x type as rows
normalize_columns = function (df, datatitle, coltitle, cols) {
  header = names(df)
  keep = setdiff(header,cols)
  new = data.frame()
  for (i in cols) {
    temp = df[c(keep,i)]
    names(temp) <- c(keep,datatitle)
    temp[coltitle] <- i
    new <- rbind(new, temp) }
  return(new) }

write.tsv = function (df, fn, ...) write.table(df, fn,
  sep="\t", row.names=FALSE, quote=FALSE, ...)

# logit transform and its inverse
logit = Vectorize(function (p) { return(log(p) - log(1 - p)) })
expit = Vectorize(function (x) { return(exp(x)/(exp(x) + 1)) })

# geometric mean
geomean = function (x) exp(mean(log(x)))

# convert selection to units of percent per period
growth_units = function (time_units_per_timestep)
  Vectorize(function(x) { return((exp(x/time_units_per_timestep) - 1)*100) })

# no zeros
nozero = function (x) { return(x[x!=0]) }

countunique = function (x) length(unique(x))
fac2int = function (x) as.numeric(as.character(x))

# logit-evenly-spaced breaks
logit_breaks = function (minp, nbins) 
  expit(seq(logit(minp), logit(1-minp), length.out=nbins-1))

merge_wildtype = function(ts, wtt, types) {
  ts$wt = ts$type %in% types
  tswt = subset(ts, wt)
  ts = subset(ts, !wt)
  tswt = transform(tswt, count = ave(count, gen, FUN=sum), type = wtt)
  gtc = c("gen","type","count")
  ts = rbind(ts[,gtc], unique(tswt[, gtc]))
  return (ts) }

annotate_frequencies = function(ts) {
  ts = transform(ts, ps = ave(count, gen, FUN=sum))
  ts = transform(ts, freq = count/ps)
  return (ts) }


# Impute missing data into a timeseries, replacing absent values with an
# unobservable frequency. The argument missingdf is a function that takes a
# data frame with columns gen, type, prevcount, nextcount that describes
# missing data that is adjacent to data that is present. It should return a
# vector (another column of the same dataframe) describing the imputed counts.
impute_missing = function (ts, missingdf) {
  library(sqldf) # slow to load, load only as necessary
  # merge previous, next and current generation counts 
  # into a complete gen x type timeseries
  df = sqldf(paste(
    "SELECT gts.gen as gen, gts.type as type,",
      "ts.count as count,",
      "tsp.count as prevcount,",
      "tsn.count as nextcount",
    "FROM",
      "(SELECT * FROM",
        "(SELECT distinct gen FROM ts) CROSS JOIN",
        "(SELECT distinct type FROM ts)) AS gts",
    "LEFT OUTER JOIN ts AS tsp ON",
      "tsp.gen + 1 = gts.gen AND tsp.type = gts.type",
    "LEFT OUTER JOIN ts AS tsn ON",
      "tsn.gen - 1 = gts.gen AND tsn.type = gts.type",
    "LEFT OUTER JOIN ts ON",
      "ts.gen = gts.gen AND ts.type = gts.type")) #2m43s
  imputable = subset(df, is.na(count) & !(is.na(nextcount) & is.na(prevcount)))
  imputable$count = missingdf(imputable)
  res = sqldf(paste(
    "SELECT df.gen as gen, df.type as type,",
      "CASE WHEN df.count IS NULL THEN ic.count ELSE df.count END AS count",
    "FROM df",
    "LEFT OUTER JOIN imputable AS ic",
    "ON ic.gen = df.gen AND ic.type = df.type",
    "ORDER BY gen ASC, count DESC, type ASC")) #46s
  return(subset(res, !is.na(count))) }

# add gen, type, count rows to a timeseries so that sum(count) by generation in
# the timeseries equals pss (header: gen, pop_size), adding the minimum number
# of rows (imputed types) given that count <= maxcount.
pad_counts = function (ts, pss, maxcount) {
  ts = merge(ts, pss)
  ts = transform(ts, ps = ave(count, gen, FUN=sum))
  ts = transform(ts, diff = pop_size - ps)
  ts = transform(ts, nrows = diff %/% maxcount, rem = diff %% maxcount,
    nrrow = ifelse(diff %% maxcount > 0, 1, 0))
  ts = merge(ts, transform(unique(ts[,c("gen","nrows","nrrow")]),
    typen = cumsum(nrows + nrrow) - nrows - nrrow)[,c("gen","typen")])

  gen_pad_counts = function (gents) {
    gendata = unique(gents[,c("gen","nrows","rem","typen")])
    stopifnot(nrow(gendata) == 1)
    gg = gendata$gen ; nrows = gendata$nrows
    rem = gendata$rem ; typen = gendata$typen
    rows = data.frame()
    remdf = data.frame()
    if (nrows > 0) {
      rows = data.frame(gen=gg, count=maxcount,
        type=paste("IT",seq(typen, typen+nrows-1),sep="")) }
    if (rem > 0) {
      remdf = data.frame(gen=gg, count=rem,
        type=paste("IT",typen+nrows,sep="")) }
    imputed_deficit = rbind(rows,remdf)
    typen = typen + nrow(imputed_deficit)
    cols = c("gen","type","count")
    return(rbind(gents[,cols],imputed_deficit[,cols])) }
  return(as.data.frame(rbindlist(
    by(ts,ts$gen,gen_pad_counts)))) }

# compute average selection coefficient over timeseries `ts` in bins delimited
# by `breaks`, assuming delta s `ss`. Specifying mutation rate `mu` modulates
# the replacement fitness; mu is otherwise *not* empirically estimated.
# `breaks` should not include 0 or 1.
make_novelty_wp = function (ts,breaks,ss,mu=0) {
  ts$mintype = ave(ts$type, ts$gen, FUN=min)
  ts$s = with(ts,ss*(type-mintype))
  ts$ps= with(ts,ave(count,gen,FUN=sum))
  ts$p = with(ts,count/ps)
  ts$Z = with(ts,ave(exp(s)*p,gen,FUN=sum))
  ts$bin = cut(ts$p,c(0,breaks,1),labels=FALSE,right=TRUE,include.lowest=TRUE)
  wp = data.frame()
  for (byn in unique(ts$bin)) {
    subts = subset(ts, bin == byn)
    ss = log(
      sum(with(subts,
        exp(s) * p * ps / Z)) /
      sum(with(subts,
        p * ps / Z)))
    avg_nt = with(subts, length(unique(paste(gen, type))))
    wp = rbind(wp, data.frame(meanfreq = mean(subts$p),
      geomeanfreq = exp(mean(log(subts$p))), 
      avg_s = ss, bin=byn, ind = as.numeric(byn) - 1,
      avg_ntypes = avg_nt)) }
  s_rep = with(unique(ts[,c("gen","Z")]), mean(log(mu + Z)))
  wp$avg_s = wp$avg_s - s_rep
  wp$avg_ntypes = wp$avg_ntypes/length(unique(ts$gen))
  wp = wp[order(wp$ind),]
  return(wp) }

make_type_dist = function (ts) {
  ts = transform(ts, 
   typesum = ave(count, type, FUN=sum),
   ps = ave(count, gen, FUN=sum))
  dist = unique(ts[,c("type","typesum")])
  dist = transform(dist, freq = typesum / sum(typesum))
  return(dist[,c("type","freq")]) }

# make_discrete_cdf - create data frame with points for the top and bottom of
# discontinuities in a discrete cumulative distribution function
make_discrete_cdf = function(xs, ys) {
  return(data.frame(
    x=rep(xs,each=2)[1:(2*length(xs)-1)],
    y=c(1,rep(ys[1:(length(ys)-1)],each=2)))) }

# make_ecdf_df - create a dataframe in the same form as make_discrete_cdf for
# the cdf of an empirical sample
make_ecdf_df = function(vs) {
  ecd = ecdf(vs)
  xs = environment(ecd)$x
  ys = environment(ecd)$y
  return(make_discrete_cdf(xs, 1-ys)) }

make_samples_control = function (ts) {
  ts = transform(ts, 
   typesum = ave(count, type, FUN=sum),
   ps = ave(count, gen, FUN=sum))
  dist = unique(ts[,c("type","typesum")])
  dist = transform(dist, freq = typesum / sum(typesum))
  
  acc = expand.grid(type=dist$type, gen=unique(ts$gen))
  acc = merge(acc, unique(ts[,c("gen","ps")]))
  
  alltcs = by(acc, acc$gen, function (df) {
    df$count = rmultinom(1, unique(df$ps), prob=dist$freq)
    return(df) })
  acc = subset(rbindlist(alltcs), count > 0)
  return(acc[,c("gen","type","count")]) }

# outputfile naming conventions and formats
read_params = function (basefn) read.delim(sprintf("%s.params.tsv", basefn))
read_nes = function (basefn) read.delim(sprintf("%s.ne_ts.tsv", basefn))
read_bininfo = function(basefn) read.delim(sprintf("%s.bin_info.tsv", basefn))
read_pops = function(basefn) read.delim(sprintf("%s.pop_sizes.tsv", basefn))
read_residuals = function(bn) read.delim(sprintf("%s.residuals.tsv", bn))
read_udresiduals = function(bn) read.delim(sprintf("%s.update_residuals.tsv", bn))
read_bootstrap = function(bn) read.delim(sprintf("%s.bootstrap.tsv", bn))
read_breaks = function (basefn) { # convert one-column table to vector
  return(read.delim(sprintf("%s.breaks.tsv", basefn),
    header=FALSE)$V1) }

exists_bootstrap = function(bn) file.exists(sprintf("%s.bootstrap.tsv", bn))

get_param = function (params,name) { return(subset(params,param==name)$val) }

get_chi2_wt = function (params,g=NA) {
  ll = get_param(params, "ll")
  swt0ll = get_param(params, "swt0ll")
  if(is.na(g)) { 
    ne = get_param(params, "ne")
    nhmean = get_param(params, "nhmean")
    g = nhmean/ne }
  return((ll - swt0ll) / g) }

get_chi2 = function (params,g=NA) {
  ll = get_param(params, "ll")
  s0ll = get_param(params, "s0ll")
  if(is.na(g)) { 
    ne = get_param(params, "ne")
    nhmean = get_param(params, "nhmean")
    g = nhmean/ne }
  return((ll - s0ll) / g) }

# argument is a data.frame w/ "gen" col, e.g. pop_sizes, timeseries, residuals
gen_interval = function (pops) {
  gens = sort(unique(pops$gen))
  return(gens[2] - gens[1]) }

speaker_counts_df = function (ress) { return(rbind(
    transform(
      unique(ress[,c("gen","init_ps")]), type = "10", freq=10/init_ps),
    transform(
      unique(ress[,c("gen","init_ps")]), type = "100", freq=100/init_ps),
    transform(
      unique(ress[,c("gen","init_ps")]), type = "1000", freq=1000/init_ps)
  )[,c("gen","type","freq")]) }

ts_gtc_to_gtf = function (gtc) {
  ts = transform(gtc,
    ps = ave(count,gen,FUN=sum))
  ts = transform(ts,
    freq = count/ps)
  return(ts) }

# Decorate residuals with their bin information
annotate_residual_bins = function(ress,breaks) {
  cols = names(ress)
  ress = transform(ress,
    ibin = cut(ress$initf, c(0,breaks,1)),
    fbin = cut(ress$obsf, c(0,breaks,1)))
  ress = transform(ress,
    iind = as.integer(ibin) - 1,
    find = as.integer(fbin) - 1)
  ress = transform(ress,
    binspan = abs(iind-find))
  return(ress[,c(cols,"iind","find","binspan")]) }

# returns an "ss" dataframe with columns ind, val, lci, and uci that can
# be merged with a midbins df to create an sp for some sense of p
make_s_of_ind = function (params) {
  srep = get_param(params, "srep")
  return(transform(subset(params, !is.na(ind)),
         val = val - srep,
         uci = uci - srep,
         lci = lci - srep)) }

# input dataframe should contain bin boundaries
make_stepwise_sp = function(ss) {
  sstep <- data.frame()
  for (i in ss$ind) {
    row = subset(ss, ind == i)
    sstep <- rbind(sstep, data.frame(x=row$lbb,y=row$val))
    sstep <- rbind(sstep, data.frame(x=row$ubb,y=row$val)) }

  return(sstep) }

# Compute midbins and average number of types in the bin
make_mean_freq_midbins = function(ressbin,FUN=geomean) {
  # copy last-gen obsf to initf
  allgens = rbind(ressbin,
    transform(subset(ressbin, gen == max(gen)),
      initf=obsf, iind=find, gen=max(gen) + gen_interval(ressbin)))
  midbins = unique(transform(allgens,
      ind = iind,
      midbin = ave(initf,iind,FUN=FUN))[,c("ind","midbin")])
  ntypes = with(allgens, tapply(paste(gen, type), iind, length)) /
    length(unique(allgens$gen))
  midbins = merge(subset(midbins, !is.na(ind)),
    data.frame(ind=names(ntypes),avg_ntypes=ntypes))
  return(midbins[order(midbins$ind),]) }

make_logmidbins = function(minf,maxf,breaks) {
  breaks = c(minf, breaks, maxf)
  return(data.frame(
    ind=seq(0,length(breaks)-2),
    midbin=sqrt(breaks[2:length(breaks)]*breaks[1:length(breaks) - 1]),
    lbb=breaks[1:(length(breaks)-1)],
    ubb=breaks[2:(length(breaks))])) }

make_params_logmidbins = function(params,breaks) {
  minf = ifelse(length(get_param(params, "minf")) == 0,
    get_param(params, "dminf"),
    get_param(params, "minf"))
  minf = ifelse(minf == min(breaks), 0.5*minf, minf)
  return(make_logmidbins( minf, get_param(params, "maxf"), breaks)) }

vquantile = function (q) function (x) quantile(x, q, na.rm=TRUE)

make_bootstrap_cis = function (bootstrap,FUN=mean) {
  # Subtract each run's replacement fitness.
  bootstrap = transform(merge(
      subset(bootstrap, !is.na(ind)),
      subset(bootstrap, param=="srep")[,c("run","val")],
      by = "run"),
    val = val.x - val.y)
  # set uci and lci based on quantiles
  bootstrap = transform(bootstrap,
    lci = ave(val, param, FUN=vquantile(0.025)),
    uci = ave(val, param, FUN=vquantile(0.975)),
    val = ave(val, param, FUN=FUN))
  bootstrap = unique(bootstrap[,c("ind","val","lci","uci")])
  return(bootstrap[order(bootstrap$ind),]) }

make_avg_s_cis = function (avg_s_of_ind) {
  wp = transform(avg_s_of_ind,
    lci = ave(avg_s, ind, FUN=vquantile(0.025)),
    uci = ave(avg_s, ind, FUN=vquantile(0.975)),
    meanavg_s = ave(avg_s, ind, FUN=mean),
    meanavg_ntypes = ave(avg_ntypes, ind, FUN=mean),
    geomeanavg_ntypes = ave(avg_ntypes, ind, FUN=geomean),
    meanfreq = ave(meanfreq, ind, FUN=mean),
    geomeanfreq = ave(geomeanfreq, ind, FUN=geomean))
  wp$val = wp$meanavg_s
  return(unique(wp[,c("val","lci","uci","ind","meanavg_s",
    "meanavg_ntypes","geomeanavg_ntypes","meanfreq","geomeanfreq")])) }

filewise_make_avg_s_cis = function (wplog) as.data.frame(rbindlist(
    by(wplog,wplog$file,make_avg_s_cis),idcol="file"))

# Note that sp needs an extra avg_ntypes column than it usually has,
# hackishly obtainable from e.g. make_novelty_wp or inf.bin_types.tsv
annotate_errs_avgwp_sp = function (avgwps,sp) {
  avgwps = merge(avgwps[,!(names(avgwps) %in% c("val","lci","uci"))],
    sp[,c("ind","val","lci","uci","avg_ntypes")])
  avgwps = transform(avgwps,
    serr2 = (meanavg_s - val)^2,
    ntypeserr2 = log(meanavg_ntypes/avg_ntypes)^2)
  avgwps = transform(avgwps,
    smeanerr2 = ave(serr2, file, FUN=mean),
    ntypesmeanerr2 = ave(ntypeserr2, file, FUN=mean) )
  return(avgwps) }

annotate_s_mu_of_file = function (wps) {
  matches = regmatches(wps$file,
    regexec("s([0-9.-]*)-mu([0-9.-]*)\\.", wps$file))
  wps$s = as.numeric(sapply(matches, "[", 2))
  wps$mu = as.numeric(sapply(matches, "[", 3))
  return(wps) }

