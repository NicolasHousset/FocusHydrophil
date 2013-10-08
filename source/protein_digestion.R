

# projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
yeastPath <- "/data/Protein/Yeast"

sequences <- data.table(read.table(file = paste0(projectPath, yeastPath, "/Yeast_01TS.txt"),
                       sep = "", col.names=c("Type","Sequence","Probability"),
                       fill = TRUE
                       ))
# Fill allows to set empty values for lines with a different number of columns

# A bit of reprocessing to make the data readable
sequences[Type!="PEPTIDE", Sequence := Type]
sequences[Type!="PEPTIDE", Type := "PROTEIN"]

sequences[, index := 1:NROW(sequences)]
# We consider the starting points of each protein
list_index <- sequences[Type=="PROTEIN", index]
list_index[NROW(list_index)+1] <- NROW(sequences)+1
# This for will link each peptide to its "parent" protein
for (i in 1:(NROW(list_index)-1)){
  print(i)
  sequences[list_index[i]:(list_index[i+1]-1), parentIndex := as.character(i)]
}

setkey(sequences, parentIndex)
# Then we can split the data into protein and peptide databases
proteins <- sequences[Type=="PROTEIN"]
peptides <- sequences[Type=="PEPTIDE"]

setkey(proteins, parentIndex)
setkey(peptides, parentIndex)

# Studying the length of the protein is not interesting
proteins[, length := nchar(as.character(Sequence))]

ggplot(proteins, aes(length)) + geom_histogram()
proteins[length < 100, length_class := 1]
proteins[length >= 100 & length < 400, length_class := 2]
proteins[length >= 400 & length < 1600, length_class := 3]
proteins[length >= 1600 & length < 6400, length_class := 4]
table(proteins[, length_class])



nbPepPerProt <- summary(peptides[, factor(parentIndex)], maxsum = 7000)
test <- as.data.frame(nbPepPerProt)
ggplot(test, aes(nbPepPerProt))+geom_histogram()

# At the global level, we count the occurences of each peptide
# Purpose : identifying unique peptides
# On Yeast they are most of the time unique though. No idea what will happen on human yet.
setkey(peptides, Sequence)
freqPep <- summary(peptides[, factor(Sequence)], maxsum = 250000)
setkey(peptides, parentIndex)
setkey(proteins, parentIndex)

# Function does not work, grrr
myFunction <- function(parentIndex,...){
  return(summary(factor(freqPep[peptides[parentIndex][,Sequence]])))
}
proteins[, nbUniqPep := summary(factor(freqPep[peptides[parentIndex][,Sequence]])), allow.cartesian = TRUE]
list_nb <- summary(factor(freqPep[peptides[,Sequence, by = parentIndex]]))
list_nb

freqPep[peptides["1234"][,Sequence]]==1
test <- data.table(data.frame(freqPep), keep.rownames=TRUE)
setkey(test, rn)
test[as.character(peptides[as.character(5)][,Sequence])]
truc <- NROW(test[as.character(peptides[as.character(5)][,Sequence])][freqPep==1])
for(i in 1:NROW(proteins)){
  print (i)
  proteins[as.character(i), nbUniqPep := summary(factor(freqPep[peptides[as.character(i)][,Sequence]]))]
}
for(i in 1:NROW(proteins)){
  print (i)
  proteins[as.character(i), nbUniqPep := NROW(test[as.character(peptides[as.character(5)][,Sequence])][freqPep==1])]
}
for(i in 1:NROW(proteins)){
  print (i)
  proteins[as.character(i), nbPep := NROW(peptides[as.character(i)])]
}

NROW(proteins[nbUniqPep>=1])
table(proteins[, nbUniqPep])
table(proteins[, nbPep])

test <- peptides[proteins[nbPep==1, parentIndex]]
sequences[test[is.na(Probability),parentIndex]]

