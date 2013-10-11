

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
yeastPath <- "/data/Protein/Yeast/Yeast_01TS.txt"

test <- myFunction("/mnt/compomics/Nicolas/R_Projects/FocusHydrophil","/data/Protein/Yeast/Yeast_01TS.txt")
test2 <- myFunction("/mnt/compomics/Nicolas/R_Projects/FocusHydrophil","/data/Protein/Yeast/Yeast_OneMiss.txt")

myFunction <- function(projectPath, proteinsPath,...){
  sequences <- data.table(read.table(file = paste0(projectPath, proteinsPath),
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
  # It takes time but it is necessary to speed up the computing afterwards
  for (i in 1:(NROW(list_index)-1)){
    sequences[list_index[i]:(list_index[i+1]-1), parentIndex := as.character(i)]
  }
  
  setkey(sequences, parentIndex)
  # Then we can split the data into protein and peptide databases
  proteins <- sequences[Type=="PROTEIN"]
  peptides <- sequences[Type=="PEPTIDE"]
  
  # To remove the last character (a ':'), useless and interfering with ELUDE
  # as.character is needed because strings are automatically converted to factors, and nchar does not do the conversion itself
  peptides[, Sequence := substr(Sequence, 1, nchar(as.character(Sequence))-1)]
  
  setkey(proteins, parentIndex)
  setkey(peptides, parentIndex)
  
  # We build 2 inside databases
  # One counts the number of peptides per protein
  # The other counts the number of proteins per protein
  # The point ? Not losing proteins with no peptide
  # At the end we add the information in the protein DB
  fooDT <- peptides[, list(nbPep=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbPep), nbPep := 0L]
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
  
  # At the global level, we count the occurences of each peptide
  # Purpose : identifying unique peptides
  # On Yeast they are most of the time unique though. No idea what will happen on human yet.
  freqPep <- summary(peptides[, factor(Sequence)], maxsum = 250000)
  
  # This version is better, it really counts the number of unique peptides
  freqPepDT <- data.table(data.frame(freqPep), keep.rownames=TRUE)
  setkey(freqPepDT, rn)
  setkey(peptides, Sequence)
  # We add the information of peptide frequency in the peptides DB
  peptides <- peptides[freqPepDT]
  
  fooDT <- peptides[freqPep==1, list(nbUniqPep=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbUniqPep), nbUniqPep := 0L]
  fooDT[, nbProt := NULL]
  table(fooDT[,nbUniqPep])
  proteins <- proteins[fooDT]
  return(proteins)
}






