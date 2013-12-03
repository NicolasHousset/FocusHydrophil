# Retention time is predicted, now we can play around with the notion of hydrophobicity and number of peptides satisfying the condition

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
proteinsPath <- "/data/Protein/Yeast"
imageFile <- "/R_Image/test.RData"

load(file = paste0(projectPath, proteinsPath, imageFile))

load(file = test)


fooDT <- peptides[freqPep==1 & Predicted_RT <= 0.5, list(nbHyd050=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
fooDT[is.na(nbHyd050), nbHyd050 := 0L]
fooDT[, nbProt := NULL]
proteins <- proteins[fooDT]

ggplot(peptides, aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01)
ggplot(peptides[Probability < 0.15], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.15 & Probability < 0.25], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.25 & Probability < 0.35], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.35 & Probability < 0.45], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.45 & Probability < 0.55], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.55 & Probability < 0.65], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)
ggplot(peptides[Probability >= 0.65 & Probability < 0.75], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)

ggplot(peptides[Probability > 0.75], aes(Predicted_RT)) + geom_histogram(aes(y= ..density..), binwidth=0.01) + xlim(-0.25,1.1)

ggplot(peptides, aes(Probability)) + geom_histogram(aes(y= ..density..), binwidth=0.01)
# For each level of retention time threshold, we count the number of peptides satisfying the criteria per protein
prefix <- "rtlt_"
for(i in 1:100){
  j <- (i/100) * 2 -1
  condition <- paste0("freqPep==1 & Predicted_RT <= ",j, " & Predicted_RT > -1")
  varCondition <- paste0(prefix,i)
  eval(parse(text=paste0("proteins[,",prefix,i," := NULL]")))
  eval(parse(text=paste0("fooDT <- peptides[",condition,", list(",varCondition,"=.N),keyby=parentIndex][proteins[,list(nbProt=.N),keyby=parentIndex]]")))
  eval(parse(text=paste0("fooDT[is.na(",prefix,i,"), ",prefix,i," := 0L]")))
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
}

a <- data.table(1:100)

for(i in 1:5){
  for(j in 1:100){
    eval(parse(text=paste0("a[", j,", atLeast_", i, " := NROW(proteins[rtlt_", j," >= ",i,"])]")))
  }
}

a[40, atLeast_3] / 16625
png(filename = paste0(projectPath,plotPath,"/B_protsTheo.png"),
    width = 800, height = 800, units = "px")
ggplot(a, aes(V1, atLeast_2)) + geom_line() + xlim(0,100) + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()

peptides <- peptides[freqPep==1]

setkey(peptides, Predicted_RT)
peptides[, nbPepIdentified := 1:NROW(peptides)]
ggplot(peptides, aes(Predicted_RT, nbPepIdentified)) + geom_line() + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
a[100, atLeast_3] / 6621

ggplot(proteins, aes(rtlt_0.25)) + geom_histogram(binwidth=3)
a[1, nbPep_0.01 := NROW(proteins[rtlt_0.01 >= 1])]

