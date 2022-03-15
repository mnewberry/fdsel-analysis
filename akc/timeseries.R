## Convert source data to timeseries readable by fdsel
#
# timeseries.tsv is the most raw plausible interpretation of AKC.csv; other
# timeseries are generated with other treatments or interpretations of the
# data. Timeseries with -d are named according to what data is removed, eg,
# timeseries-diy.tsv (delta initial years) omits data for the year of a breed's
# introduction which is contaminated by other factors (partial years or initial
# popularity effects)
# timeseries-diyo.tsv omits iy as well as some outliers.

options(stringsAsFactors = FALSE)
counts <- read.csv("inp/AKC.csv", header=TRUE, sep =";")

counts2 <- counts[,-1]
rownames(counts2) <- counts[,1]

source("../common.R")

# Exclude these types even in "raw" timeseries, since they are incomplete data
# subcategorizations of other types that exist in the dataset. A lot of special
# effort would be required to include these timeseries. Excluding these retains
# the timeseries that sum them, eg. Spaniel (all Cocker), which is retained,
# sums all the incomplete Cocker Spaniels.
exclude = c(
  "Dachshund-XX",
  "Dachshund-XX  (long-haired)",
  "Daschund-XX (smooth)",
  "Daschund-XX (wire-haired)",
  "English Toy Spaniel (Blenheim and Prince Charles)",
  "English Toy Spaniel (King Charles and Ruby)",
  "Fox Terrier-XX",
  "Fox Terrier (Smooth)",
  "Fox Terrier (Wire)",
  "Manchester Terrier-XX",
  "Manchester Terrier (Toy)",
  "Poodle-XX",
  "Poodle (miniature)",
  "Poodle (standard)",
  "Poodle (toy)",
  "Spaniel (Cocker)",
  "Spaniel (Cocker), American type, any solid color other than black, including black and tan",
  "Spaniel (Cocker), American Type, parti-color",
  "Spaniel (Cocker), American Type, solid black color",
  "Spaniel (English Cocker)")

norm <- normalize_columns(counts, "count", "year", names(counts2))

make_gtc = function(counts){
  gen <- as.numeric(gsub("X","",counts$year))
  type <- as.character( counts$Breed )
  count<- as.numeric(counts$count)
  df <- data.frame(gen=gen,type=type,count=count)
  df = subset(df, !is.na(count) & !(type %in% exclude) & count > 0)
  return(df)}

df = make_gtc(norm)
write.tsv(df, "out/timeseries.tsv")

## It is possible to impute ones or zeros for missing data, but this creates
## many more problems than it solves, so we leave it commented here.
# 
#complete_ts = expand.grid(gen=unique(df$gen),type=unique(df$type),
#  stringsAsFactors=FALSE)
#complete_ts = transform(complete_ts, gentype = paste(gen, type))
#chooseany = function (x) { return(x[!is.na(x)][1]) }
#inc_ts = transform(df, 
#  gentype = paste(gen, type),
#  mingen = ave(gen, type, FUN=min),maxgen = ave(gen, type, FUN=max))
#complete_ts = transform(
#  merge(complete_ts, inc_ts[,c("gentype","count","mingen","maxgen")],
#    by="gentype", all.x=TRUE),
#  oldcount = count,
#  mingen=ave(mingen,type,FUN=chooseany),maxgen=ave(maxgen,type,FUN=chooseany))
#inc_ts = transform(complete_ts, gentype = paste(gen-1, type), nextcount=count)
#complete_ts = merge(complete_ts, inc_ts[,c("gentype","nextcount")], all.x=TRUE)
#inc_ts = transform(complete_ts, gentype = paste(gen+1, type), prevcount=count)
#complete_ts = merge(complete_ts, inc_ts[,c("gentype","prevcount")], all.x=TRUE)
#complete_ts = transform(complete_ts,
#  count = ifelse(
#    gen > mingen & is.na(count) & !(is.na(nextcount) & is.na(prevcount)),
#    ifelse(!is.na(nextcount), 1, 0),
#    count))
#
#imp = subset(complete_ts[,c("gen","type","count")], !is.na(count))
#
#write.table(df, "out/timeseries-imp.tsv", sep="\t", col.names=TRUE,
#  row.names=FALSE, quote=FALSE)

iy = subset(transform(df, iyr = ave(gen,type,FUN=min)),
          iyr == gen & iyr > 1926)
iyears = with(iy, paste(gen, type))
dfdiy = subset(df, !(paste(gen, type) %in% iyears))

# Write a timeseries that omits the count for each year of introduction for
# breeds not present at the beginning of the dataset. This removes a lot of
# noisy, seemingly spurious data.
write.tsv(dfdiy, "out/timeseries-diy.tsv")

# Write a timeserise delta seven jumpy/discontinuous outlier transitions:
dfdiyo = dfdiy
# Remove certain bad data years
badcells = c("1939 Poodle", "1958 Manchester Terrier")
dfdiyo = subset(dfdiyo, !(paste(gen, type) %in% badcells))
# Cut certain transitions
dfdiyo$type[dfdiyo$type == "Manchester Terrier" & dfdiyo$gen <= 1945] = 
  "Manchester Terrier PRE"
dfdiyo$type[dfdiyo$type == "Greyhound" & dfdiyo$gen <= 1934] = 
  "Greyhound PRE"
dfdiyo$type[dfdiyo$type == "Curly-Coated Retriever" & dfdiyo$gen <= 1932] = 
  "Curly-Coated Retriever PRE"
write.tsv(dfdiyo, "out/timeseries-diyo.tsv")

diyosamples = make_samples_control(read.delim("out/timeseries-diyo.tsv"))
write.tsv(diyosamples, "out/timeseries-diyosamples.tsv")

# Write the same timeseries but drop years less than 1945
write.tsv(subset(dfdiyo, gen > 1945), "out/timeseries-diyo45.tsv")
