

# projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"

projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"
yeastPath <- "/data/Protein/Yeast"

peptides <- data.table(read.table(file = paste0(projectPath, yeastPath, "/Yeast_01TS.txt"),
                       sep = "", col.names=c("Type","Sequence","Probability"),
                       fill = TRUE
                       ))
# Fill allows to set empty values for lines with a different number of columns

peptides[Type!="PEPTIDE", Sequence := Type]
peptides[Type!="PEPTIDE", Type := "PROTEIN"]

length_protein <- data.table(as.data.frame(table(peptides[Type=="PROTEIN",nchar(as.character(Sequence))]), rownames="length"))

ggplot(length_protein, aes(Var1, Freq)) + geom_point()

proteins <- peptides[Type=="PROTEIN"]
sequences <- peptides[Type=="PEPTIDE"]

proteins[, length := nchar(as.character(Sequence))]

ggplot(proteins, aes(length)) + geom_histogram()
proteins[length < 100, length_class := 1]
proteins[length >= 100 & length < 400, length_class := 2]
proteins[length >= 400 & length < 1600, length_class := 3]
proteins[length >= 1600 & length < 6400, length_class := 4]
table(proteins[, length_class])

peptides[, index := 1:NROW(peptides)]
list_index <- peptides[Type=="PROTEIN", index]
list_index[NROW(list_index)+1] <- NROW(peptides)-1
for (i in 1:(NROW(list_index)-1)){
  print(i)
  peptides[list_index[i]:(list_index[i+1]-1), parent := peptides[list_index[i],Sequence]]
}

nbPepPerProt <- summary(peptides[, parent], maxsum = 7000)
