projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
testData <- "/data/DeNovoDigest.txt"
savePredict <- "/data/predictions.out"
loadModel <- "/data/modelHydrophil.model"
testData <- shQuote(paste0(projectPath, testData))
savePredict <- shQuote(paste0(projectPath, savePredict))
loadModel <- shQuote(paste0(projectPath, loadModel))

verbFlag <- " -v "
testFlag <- " -e "
savePredictFlag <- " -o "
ignoreNewTestPTMFlag <- " -p "
verbLevel <- " 5"
loadModelFlag <- " -l "

# Predicting retention time on the denovo digest
eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, testFlag,
                     testData, loadModelFlag, loadModel,
                     savePredictFlag, savePredict, ignoreNewTestPTMFlag)
shell(strCommand, translate = TRUE, wait = TRUE)

results <- data.table(read.table(file=paste0(projectPath, "/data/predictions.out"), header = TRUE, sep = "\t"))
setkey(results, Peptide)
setkey(peptides, Sequence)
peptides <- peptides[results]

