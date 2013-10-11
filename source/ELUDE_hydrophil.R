projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"

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
eludePath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, testFlag,
                     testData, loadModelFlag, loadModel,
                     savePredictFlag, savePredict, ignoreNewTestPTMFlag)
strCommand <- "elude"
# Currently not very portable between Windows and Unix, shell is for Windows, system for Unix
# Try system2 ? Grammar is very different though
shell(strCommand, translate = TRUE, wait = TRUE)
system(strCommand, wait = TRUE)

# Command works on UNIX, needs to be tested on Windows though
# Problem is that elude is not recognized as a system command
system2(strCommand, args = c(verbFlag, verbLevel, testFlag, testData,
                             loadModelFlag, loadModel, savePredictFlag, savePredict,
                             ignoreNewTestPTMFlag))

results <- data.table(read.table(file=paste0(projectPath, "/data/predictions.out"), header = TRUE, sep = "\t"))
setkey(results, Peptide)
setkey(peptides, Sequence)
peptides <- peptides[results]

