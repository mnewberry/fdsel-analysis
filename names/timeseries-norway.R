source("../common.R")

read_norway = function (filename) read.delim(filename, 
  na.string=c(".",".."), sep=";",strip.white=TRUE)

tsF <- read_norway("inp/ssb.no-Personer-1880-2019F.csv")
tsM <- read_norway("inp/ssb.no-Personer-1880-2019M.csv")

tsFnorm <- normalize_columns(tsF, "count", "Xyear",
  names(tsF)[2:length(names(tsF))])
tsMnorm <- normalize_columns(tsM, "count", "Xyear",
  names(tsM)[2:length(names(tsM))])

# filter NA values
tsFnorm <- tsFnorm[complete.cases(tsFnorm),]
tsMnorm <- tsMnorm[complete.cases(tsMnorm),]

xyeartoint = Vectorize(function (x) { str = as.character(x)
  return(as.integer(substr(str,nchar(str)-3, nchar(str)))) })

tsFnorm <- data.frame(type=tsFnorm$first.name,
  count=tsFnorm$count,gen= xyeartoint(tsFnorm$Xyear))
tsMnorm <- data.frame(type=tsMnorm$first.name,
  count=tsMnorm$count,gen= xyeartoint(tsMnorm$Xyear))
tsFnorm$type = paste(tsFnorm$type, "F", sep=",")
tsMnorm$type = paste(tsMnorm$type, "M", sep=",")

ts = rbind(tsFnorm, tsMnorm)[,c("gen","type","count")]
ts = transform(ts, ps = ave(count, gen, FUN=sum))
pss = unique(ts[,c("gen","ps")])
cencounts = transform(pss,
  type="CENSORED", count=round(0.05*ps))
ts = rbind(ts, cencounts)[,c("gen","type","count")]

# Throw away 1880-1945 as it is mostly NA due to sparse data.
tsd45 = subset(ts, gen >= 1946)

write.tsv(ts, "out/timeseries-norway.tsv")
write.tsv(tsd45, "out/timeseries-norwayd45.tsv")
write.tsv(subset(tsd45, type != "CENSORED"), "out/timeseries-norwayrawd45.tsv")
typedist = make_type_dist(subset(tsd45, type != "CENSORED"))
write.tsv(typedist, "out/norwayd45.typefreqs.tsv")
