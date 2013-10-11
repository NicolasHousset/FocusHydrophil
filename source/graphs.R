# Retention time is predicted, now we can play around with the notion of hydrophobicity and number of peptides satisfying the condition

load(file = paste0(projectPath, proteinsPath, imageFile))

load(file = test3)


fooDT <- peptides[freqPep==1 & Predicted_RT <= 0.5, list(nbHyd050=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
fooDT[is.na(nbHyd050), nbHyd050 := 0L]
fooDT[, nbProt := NULL]
proteins <- proteins[fooDT]


prefix <- "rtlt_"
for(i in -20:100){
  j <- i/100 + 0.2
  condition <- paste0("freqPep==1 & Predicted_RT <= ",j)
  varCondition <- paste0(prefix,j)
  eval(parse(text=paste0("proteins[,",prefix,j," := NULL]")))
  eval(parse(text=paste0("fooDT <- peptides[",condition,", list(",varCondition,"=.N),keyby=parentIndex][proteins[,list(nbProt=.N),keyby=parentIndex]]")))
  eval(parse(text=paste0("fooDT[is.na(",prefix,j,"), ",prefix,j," := 0L]")))
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
}

a <- data.table(1:120)

for(i in 1:5){
  for(j in -20:99){
    k <- j/100 + 0.2
    eval(parse(text=paste0("a[", j+21,", atLeast_", i, " := NROW(proteins[rtlt_", k," >= ",i,"])]")))
  }
}
a[40, atLeast_3] / 16625
ggplot(a, aes(V1, atLeast_3)) + geom_line() + xlim(0,100) + ylim(0,20265)
a[1, nbPep_0.01 := NROW(proteins[rtlt_0.01 >= 1])]

