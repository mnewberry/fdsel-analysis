options(stringsAsFactors = FALSE)
source("../common.R")

tsraw <- read.delim("out/timeseries-ssaraw.tsv")
write.tsv(make_type_dist(tsraw),
  "out/ssaraw.typefreqs.tsv")
tsrawM = subset(tsraw, substr(type,nchar(type),nchar(type)) == "M")
tsrawF = subset(tsraw, substr(type,nchar(type),nchar(type)) == "F")
ntypes = length(unique(tsraw$type))
ntypesM = length(unique(tsrawM$type))
ntypesF = length(unique(tsrawF$type))
pop = sum(tsraw$count)
popM = sum(tsrawM$count)
popF = sum(tsrawF$count)
tsrawsums = data.frame(
  param=c("ntypes","pop","ntypesM","ntypesF","popM","popF"),
  val = c(ntypes,  pop,  ntypesM,  ntypesF,  popM,  popF)) 
write.tsv(tsrawsums, "out/ssaraw.sums.tsv")

tsraw = annotate_frequencies(tsraw)
write.tsv(subset(tsraw[,c("gen","freq")], runif(nrow(tsraw)) < 0.01),
  "out/ssaraw.genfreqs1%.tsv")
write.tsv(subset(tsraw[,c("gen","freq")], runif(nrow(tsraw)) < 0.05),
  "out/ssaraw.genfreqs5%.tsv")
ts <- read.delim("out/timeseries-ssaCd35dCE.tsv")

# === Compute biblical and non-biblical name treatments ===
bnames <- read.delim("inp/bible-names.txt",header=FALSE)$V1
# put names in SSA format indiscriminately w.r.t gender:
gbnames = c(paste(bnames, "M", sep=","), paste(bnames, "F", sep=","))

typedist = make_type_dist(subset(ts, type != "CENSORED"))
typedist$biblical = typedist$type %in% gbnames
write.tsv(typedist, "out/ssad35dCE.typefreqs.tsv")

bib_annotate = function (ts) {

  ts = transform(ts, biblical = type %in% gbnames)
  ts = transform(ts,
    bibcat = ifelse(biblical, type, "NONBIBLICAL"),
    nobcat = ifelse(biblical, "BIBLICAL", type))
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

bibts = bib_annotate(ts)

bib = bib_slice(bibts)
nob = nob_slice(bibts)

write.tsv(bib, "out/timeseries-ssabib.tsv")
write.tsv(nob, "out/timeseries-ssanob.tsv")

ssanames = unique(ts$type)
namesMCEN = ssanames[substr(ssanames,nchar(ssanames),nchar(ssanames)) != "F"]
namesFCEN = ssanames[substr(ssanames,nchar(ssanames),nchar(ssanames)) != "M"]

tsM = merge_wildtype(ts, "FCEN", namesFCEN)
tsF = merge_wildtype(ts, "MCEN", namesMCEN)

write.tsv(tsM, "out/timeseries-ssaM.tsv")
write.tsv(tsF, "out/timeseries-ssaF.tsv")
