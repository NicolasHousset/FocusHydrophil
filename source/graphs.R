# Retention time is predicted, now we can play around with the notion of hydrophobicity and number of peptides satisfying the condition

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
proteinsPath <- "/data/Protein/MouseGold"
imageFile <- "/R_Image/test.RData"

load(file = paste0(projectPath, proteinsPath, imageFile))

load(file = test3)


fooDT <- peptides[freqPep==1 & Predicted_RT <= 0.5, list(nbHyd050=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
fooDT[is.na(nbHyd050), nbHyd050 := 0L]
fooDT[, nbProt := NULL]
proteins <- proteins[fooDT]


# For each level of retention time threshold, we count the number of peptides satisfying the criteria per protein
prefix <- "rtlt_"
for(i in 1:100){
  j <- i/100
  condition <- paste0("freqPep==1 & Predicted_RT <= ",j)
  varCondition <- paste0(prefix,j)
  eval(parse(text=paste0("proteins[,",prefix,j," := NULL]")))
  eval(parse(text=paste0("fooDT <- peptides[",condition,", list(",varCondition,"=.N),keyby=parentIndex][proteins[,list(nbProt=.N),keyby=parentIndex]]")))
  eval(parse(text=paste0("fooDT[is.na(",prefix,j,"), ",prefix,j," := 0L]")))
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
}

a <- data.table(1:100)

for(i in 1:5){
  for(j in 1:100){
    k <- j/100 
    eval(parse(text=paste0("a[", j,", atLeast_", i, " := NROW(proteins[rtlt_", k," >= ",i,"])]")))
  }
}

a[40, atLeast_3] / 16625
ggplot(a, aes(V1, atLeast_3/NROW(proteins))) + geom_line() + xlim(0,100) + ylim(0,1)
a[100, atLeast_3] / 6621

ggplot(peptides, aes(Predicted_RT)) + geom_histogram(binwidth=0.01)
a[1, nbPep_0.01 := NROW(proteins[rtlt_0.01 >= 1])]

