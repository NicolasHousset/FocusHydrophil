

# projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
yeastPath <- "/data/Protein/Yeast"

peptides <- data.table(read.table(file = paste0(projectPath, yeastPath, "/Yeast_OneMiss.txt"),
                       sep = "", col.names=c("Type","Sequence","Probability"),
                       fill = TRUE
                       ))