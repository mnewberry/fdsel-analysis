options(stringsAsFactors = FALSE)
source("../common.R")

ts = read.delim("out/france/nat2018.csv",sep=";",na.strings=c())

# Immediately throw away data prior to 1946: "l’exhaustivité n’est pas garantie
# sur toute la période, notamment pour les années antérieures à 1946"
ts = subset(ts, annais >= "1946")
nentries = nrow(ts)
totpop = sum(ts$nombre)
pops = unique(transform(ts, ps=ave(nombre,annais,FUN=sum))[,c("annais","ps")])

# Years XXXX count types that were censored due to low counts in some year
# "Les effectifs [counts] des prénoms remplissant [passing] la condition 2
# [occuring at least 20 times after 1946] mais pas la condition 3 [at least
# three times in any given year] sont regroupés (pour chaque sexe et chaque
# prénom) dans un enregistrement dont le champ année de naissance (ANNAIS)
# prend la valeur «XXXX»."
tsXXXX = subset(ts, annais == "XXXX")
ts = subset(ts, annais != "XXXX")

tsRARE = subset(ts, preusuel == "_PRENOMS_RARES")
ts = subset(ts, preusuel != "_PRENOMS_RARES")

stopifnot(nentries == nrow(tsRARE) + nrow(tsXXXX) + nrow(ts)) # no sneaks

# # Useful for establishing percent censored to use with other datasets.
# totXXXX = sum(tsXXXX$nombre)
# totcen = totXXXX + sum(tsRARE$nombre)
# totl5 = totcen + sum(subset(ts, nombre < 5)$nombre)
# totl4 = totcen + sum(subset(ts, nombre < 4)$nombre)
# cat(sprintf("Fraction of the data less then 5: %f\n", totl5 / totpop))
# cat(sprintf("Fraction of the data less then 4: %f\n", totl4 / totpop))
# cat(sprintf("Fraction of the data less then 3: %f\n", totcen / totpop))

# Put the timeseries in gen, type, count format
ts = transform(ts,
  gen = as.integer(annais),
  type = paste(sexe, preusuel, sep=","), # le ,  n'existe pas dans les prenoms
  count = nombre)[,c("gen","type","count")]

psXXXX = subset(pops, annais == "XXXX")$ps
cenperuncen = psXXXX / (totpop - psXXXX)
tscen = transform(subset(pops, annais != "XXXX"),
  gen = as.integer(annais),
  countXXXX = round(ps * cenperuncen))
tsRAREnosexe = unique(transform(tsRARE, 
  countPR = ave(nombre, annais, FUN=sum))[,c("annais","countPR")])
tscen = transform(merge(tscen, tsRAREnosexe, by="annais"),
  count = countXXXX + countPR,
  type = "CENSORED")[,c("gen","type","count")]

ts = rbind(ts, tscen)
ts = ts[order(ts$gen),]
# total population size including censored counts differ from original by a
# rounding error of at most one.
stopifnot(abs(totpop - sum(ts$count)) <= 1)

write.tsv(ts, "out/timeseries-france.tsv")
typedist = make_type_dist(subset(ts, type != "CENSORED"))
write.tsv(typedist, "out/france.typefreqs.tsv")
write.tsv(subset(ts,gen>1993), "out/timeseries-franced93.tsv")
write.tsv(subset(ts,gen<1993), "out/timeseries-francedp93.tsv")


ts5y = subset(ts, gen >= 1947 & gen <= 2016)
ts5y = transform(ts5y, y5 = (gen - 1947) %/% 5)
ts5y = transform(ts5y, y5type = paste(type, y5),
                       y5gen = y5 * 5 + 1947)
ts5y = transform(ts5y, y5count = ave(count, y5type, FUN=sum))
ts5y = unique(ts5y[,c("y5gen","type","y5count")])
names(ts5y) = c("gen","type","count")

write.tsv(ts5y, "out/timeseries-france5y.tsv")
