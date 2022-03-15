# sim-nb.R - Simulate the novelty bias model in 80-character chunks and record
# average selection coefficient in frequency bins until the rate of approach to
# equilibrium is indistinguishable from statistical error.

options(stringsAsFactors = FALSE)
source("../common.R")

# Parse command-line arguments
suppressPackageStartupMessages(library("optparse"))
cli_opts = list(
  make_option(c("-o", "--logfile"), metavar="filename",
    help="output file to log average selection coefficients"),
  make_option(c("-t", "--timeseries"), metavar="filename",
    help="timeseries basename output (-rNNNN.tsv will be appended)"),
  make_option(c("-p", "--params"), metavar="filename",
    help="simulation parameters"),
  make_option(c("-B", "--breaks"), metavar="filename",
    help="bin boundaries (headerless nlsv not including endpoints 0 and 1)"),
  make_option(c("-m", "--min-chunks"), type="integer", default=5,
    help="minimum number of chunks (at least 3)", metavar="NN"),
  make_option(c("-M", "--max-chunks"), type="integer", default=12,
    help="maximum number of chunks (at least 5)", metavar="NN"),
  make_option(c("-c", "--chunk-length"), type="integer", default=80,
    help="length of simulation chunks", metavar="NN"),
  make_option(c("-n", "--pop-size"), type="integer", metavar="NN",
    help="constant population size"),
  make_option(c("-D", "--dilation"), type="integer", default=0,
    help="number of neutral timesteps per generation", metavar="NN"),
  make_option(c("-s", "--delta-s"), metavar="val",
    help="selection benefit to new mutants"),
  make_option(c("-u", "--mu"), metavar="val",
    help="mutation rate per generation per capita"),
  make_option(c("-R", "--rngseed"), default=1, 
    help="random number generator seed", metavar="NN"))


#opt <- parse_args(args=c("-o", "out/wplog_nb-s0.0000001-mu0.0000001.tsv", "-t", "out/timeseries_nb-s0.0000001-mu0.0000001", "-p", "out/params_nb-s0.0000001-mu0.0000001", "-s", "0.0000001", "-u", "0.0000001"), OptionParser(option_list=cli_opts, description="Simulate novelty bias model in chunks"), convert_hyphens_to_underscores=TRUE)
opt <- parse_args(OptionParser(option_list=cli_opts,
    description="Simulate novelty bias model in chunks"),
  convert_hyphens_to_underscores=TRUE)

stopifnot(!is.null(opt$breaks))
infbreaks = read.delim(opt$breaks, header=FALSE)$V1

novelty_wp = function(fn) make_novelty_wp(
  read.delim(fn,as.is=TRUE), infbreaks, 
  as.numeric(opt$delta_s),
  as.numeric(opt$mu))

tsfn = function(run) sprintf("%s-r%04d.tsv", opt$timeseries, run)

# Run the simulation chunks
run = 0
# initial chunk command, starting from monomorphic initial population
simparams = paste(sep=" ",
  "-s", opt$rngseed,
  "-i", opt$params,
  "-n", opt$pop_size,
  "-g", opt$chunk_length,
  "-D", opt$dilation)
system(paste(sep=" ", "fdsel simulate -N", "-o", tsfn(0), simparams))
wplog = transform(novelty_wp(tsfn(0)),
          run = run)
# initiate the log file
write.tsv(wplog, opt$logfile)
while (TRUE) {
  run = run + 1
  Sys.sleep(0.1) # impossible to interrupt if absent
  system(paste(sep=" ",
    "fdsel simulate -N",
      "-F -T", tsfn(run - 1),
      "-o", tsfn(run),
      simparams))
  wp = transform(novelty_wp(tsfn(run)), run = run)
  wplog = rbind(wplog, wp) # this gets slow when long
  # *append* wp to output
  write.tsv(wp, opt$logfile, append=TRUE, col.names=FALSE)

  # Decide whether to continue simulating
  if(run >= max(opt$max_chunks, opt$min_chunks)) { break }
  if(run >= max(opt$min_chunks, 3)) {
    thisrun = run
    recent = subset(wplog, run >= thisrun - 5)
    awp = subset(recent, run <= thisrun - 3)
    bwp = subset(recent, run > thisrun - 3)
    bins = unique(recent$bin)
    pco = 1 - 0.5^(1/(length(bins) + 1)) # p cutoff
    if(length(unique(awp$avg_ntypes)) >= 2 & 
         length(unique(bwp$avg_ntypes)) >= 2) {
      ntypes.p = t.test(awp$avg_ntypes,bwp$avg_ntypes)$p.value 
      if (ntypes.p < pco) { next } }
    pass = TRUE 
    for (thisbin in bins) {
      thisawp = subset(awp, bin == thisbin)
      thisbwp = subset(bwp, bin == thisbin)
      if(nrow(thisawp) < 3 | nrow(thisbwp) < 3) { pass = FALSE ; break }
      if(length(unique(thisawp$avg_s)) >= 2 &
           length(unique(thisbwp$avg_s)) >= 2) {
        avg_s.p = t.test(thisawp$avg_s,thisbwp$avg_s)$p.value 
        if(avg_s.p < pco) { pass = FALSE ; break } } }
    if(pass) { break }
  }
}
