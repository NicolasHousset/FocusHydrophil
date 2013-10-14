

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"

proteinsPath <- "/data/Protein/Yeast"

modelFile <- "/data/Elude/modelHydrophil.model"
digestionFile <- "/Peptides/Digestion/Yeast_01TS.txt"
eludeInputFile <- "/Peptides/RT_Prediction/pep_01TS.txt"
predictionFile <- "/Peptides/RT_Prediction/pep_01TS_pred.txt"

test <- myFunction(projectPath = "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil",
                   proteinsPath = "/data/Protein/Yeast",
                   modelFile = "/data/Elude/modelHydrophil.model",
                   digestionFile = "/Peptides/Digestion/YeastGold_01TS.txt",
                   eludeInputFile = "/Peptides/RT_Prediction/pep_01TS.txt",
                   predictionFile = "/Peptides/RT_Prediction/pep_01TS_pred.txt",
                   imageFile = "/R_Image/test.RData")

test2 <- myFunction(projectPath = "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil",
                   proteinsPath = "/data/Protein/MouseGold",
                   modelFile = "/data/Elude/modelHydrophil.model",
                   digestionFile = "/Peptides/Digestion/MouseGold_01TS.txt",
                   eludeInputFile = "/Peptides/RT_Prediction/pep_01TS.txt",
                   predictionFile = "/Peptides/RT_Prediction/pep_01TS_pred.txt",
                    imageFile = "/R_Image/test.RData")

test3 <- myFunction(projectPath = "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil",
                    proteinsPath = "/data/Protein/HumanGold",
                    modelFile = "/data/Elude/modelHydrophil.model",
                    digestionFile = "/Peptides/Digestion/HumanGold_01TS.txt",
                    eludeInputFile = "/Peptides/RT_Prediction/pep_01TS.txt",
                    predictionFile = "/Peptides/RT_Prediction/pep_01TS_pred.txt",
                    imageFile = "/R_Image/test.RData")

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
proteinsPath <- "/data/Protein/HumanGold"
modelFile <- "/data/Elude/modelHydrophil.model"
digestionFile <- "/Peptides/Digestion/HumanGold_01TS.txt"
eludeInputFile <- "/Peptides/RT_Prediction/pep_01TS.txt"
predictionFile <- "/Peptides/RT_Prediction/pep_01TS_pred.txt"
imageFile <- "/R_Image/test.RData"

myFunction <- function(projectPath, proteinsPath, modelFile, digestionFile, eludeInputFile, predictionFile,imageFile,...){
  sequences <- data.table(read.table(file = paste0(projectPath, proteinsPath, digestionFile),
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
  # Removing peptides containing ambiguous amino acid
  peptides <- peptides[!grepl("[BJOUXZ]",Sequence)]
  
  setkey(proteins, parentIndex)
  setkey(peptides, parentIndex)
  
  # At the global level, we count the occurences of each peptide
  freqPep <- summary(peptides[, factor(Sequence)], maxsum = 1000000)
  freqPepDT <- data.table(data.frame(freqPep), keep.rownames=TRUE)
  rm(freqPep)
  setkey(freqPepDT, rn)
  setkey(peptides, Sequence)
  # We add the information of peptide frequency in the peptides DB
  peptides <- peptides[freqPepDT]
  
  # We build 2 inside databases
  # One counts the number of peptides per protein
  # The other counts the number of proteins per protein
  # The point ? Not losing proteins with no peptide
  # At the end we add the information in the protein DB
  fooDT <- peptides[, list(nbPep=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbPep), nbPep := 0L]
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
  
  
  fooDT <- peptides[freqPep==1, list(nbUniqPep=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbUniqPep), nbUniqPep := 0L]
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
  
  eludeDT <- freqPepDT[,rn]
  write.table(eludeDT, file=paste0(projectPath, proteinsPath, eludeInputFile), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  rm(eludeDT)
  
  testData <- shQuote(paste0(projectPath, proteinsPath, eludeInputFile))
  savePredict <- shQuote(paste0(projectPath, proteinsPath, predictionFile))
  loadModel <- shQuote(paste0(projectPath, modelFile))
  
  verbFlag <- " -v "
  testFlag <- " -e "
  savePredictFlag <- " -o "
  ignoreNewTestPTMFlag <- " -p "
  verbLevel <- " 5"
  loadModelFlag <- " -l "
  
  strCommand <- "elude"
  
  # Probably does not work on Windows
  system2(strCommand, args = c(verbFlag, verbLevel, testFlag, testData,
                               loadModelFlag, loadModel, savePredictFlag, savePredict,
                               ignoreNewTestPTMFlag))
  
  results <- data.table(read.table(file=paste0(projectPath, proteinsPath, predictionFile), header = TRUE, sep = "\t"))
  setkey(results, Peptide)
  setkey(peptides, Sequence)
  peptides <- peptides[results]
  
  fooDT <- peptides[freqPep==1 & Predicted_RT <= 0.5, list(nbHyd050=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbHyd050), nbHyd050 := 0L]
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
  
  fooDT <- peptides[freqPep==1 & Predicted_RT <= 0.25, list(nbHyd025=.N),keyby=parentIndex][proteins[, list(nbProt=.N), keyby=parentIndex]]
  fooDT[is.na(nbHyd025), nbHyd025 := 0L]
  fooDT[, nbProt := NULL]
  proteins <- proteins[fooDT]
  
  save(list = c("proteins", "peptides"), file=paste0(projectPath, proteinsPath, imageFile),
       compress = "gzip", compression_level = 1)
  
  return(paste0(projectPath, proteinsPath, imageFile))
}


# Compute retention times for all potential peptides, otherwise a pain to integrate into peptides DB

