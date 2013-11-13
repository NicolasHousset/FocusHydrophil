dataPath <- "//uv2522.ugent.be/compomics/Nicolas/CPTAC/"

load(file= paste0(dataPath, "CPTAC_processed.RData"))

result <- result[!is.na(rtsec)]

# CPTAC data is divided into three labs
result[grepl("200803", fileOrigin), lab := "lab1"]
result[grepl("mam", fileOrigin), lab := "lab2"]
result[grepl("Orbi2", fileOrigin), lab := "lab3"]


result[, fileOrigin.f := factor(fileOrigin)]
list_experiments <- summary(result[, factor(fileOrigin)])
# Each lab has three replicates
setkey(result, fileOrigin.f)
result[labels(list_experiments[c(1:5, 16:20, 31:34)]), rep := "rep1"]
result[labels(list_experiments[c(6:10, 21, 23, 25, 27, 29 , 35, 37, 40, 42)]), rep := "rep2"]
result[labels(list_experiments[c(11:15, 22, 24, 26, 28, 30, 36, 38, 39, 41, 43)]), rep := "rep3"]

# Each experiment measures one of those 5 samples : A, B, C, D, E
result[labels(list_experiments[c(1,6,11,16,21,22,35,36)]), sample := "A"]
result[labels(list_experiments[c(2,7,12,17,23,24,31,37,38)]), sample := "B"]
result[labels(list_experiments[c(3,8,13,18,25,26,32,39)]), sample := "C"]
result[labels(list_experiments[c(4,9,14,19,27,28,33,40,41)]), sample := "D"]
result[labels(list_experiments[c(5,10,15,20,29,30,34,42,43)]), sample := "E"]

result[, exp := factor(paste0(lab, rep))]


result[, rtsec.f := cut(rtsec, 100)]
result[, q50_Intensity := quantile(MS1_Intensity, probs=0.50), by = list(exp, sample, rtsec.f)]

setkey(result, exp)
# Histogram of the intensity
graph <- ggplot(result[c("lab1rep1","lab1rep2","lab1rep3")], aes(rtsec, log(q50_Intensity))) + geom_point() + xlim(0,10000) +  ylim(6, 15) + facet_grid(sample ~ rep)
graph 
graph %+% result[c("lab1rep1","lab1rep2","lab1rep3")]
graph %+% result[c("lab2rep1","lab2rep2","lab2rep3")]
graph %+% result[c("lab3rep1","lab3rep2","lab3rep3")]



graph <- ggplot(result[c("lab1rep1","lab1rep2","lab1rep3")], aes(rtsec)) + geom_histogram(aes(y = ..count..), binwidth = 60) + xlim(0,10000) + ylim(0,150) + facet_grid(sample ~ rep)
graph

result[list_experiments[1]]
setkey(result, lab)
ggplot(result["lab1"], aes(rtsec)) + geom_histogram(aes(y = ..density..), binwidth = 60) + xlim(0,10000)
ggplot(result["lab2"], aes(rtsec)) + geom_histogram(aes(y = ..density..), binwidth = 60) + xlim(0,10000)
ggplot(result["lab3"], aes(rtsec)) + geom_histogram(aes(y = ..density..), binwidth = 60) + xlim(0,10000)



result_filtered <- result[as.double(as.character(q_value_percolator)) < 0.01][multi_pro == FALSE]

table(result_filtered[, multi_pro])

setkey(result_filtered, fileOrigin, modseq)
countsPerProject <- unique(result_filtered)[, list(fileOrigin,modseq)]
countsPerProject[, modseq.f := factor(modseq)]

nbProjPerPeptide <- summary(countsPerProject[, modseq.f], maxsum = 1000000)
id_peptide <- 1:NROW(nbProjPerPeptide)
dt <- data.table(id_peptide)
dt[, modified_sequence := labels(nbProjPerPeptide)]
dt[, nbProjPep := -nbProjPerPeptide]
# This ordering makes most common peptides appear first
setkey(dt, nbProjPep)
# We save the rank
dt[, rank_peptide := 1:NROW(nbProjPerPeptide)]
# Number of projects is brought back to a positive number
dt[, nbProjPep := -nbProjPep]


result_filtered_2 <- result_filtered[rtsec < 4000]
setkey(result_filtered_2, fileOrigin, protein)
countsPerProtein <- unique(result_filtered_2)[, list(fileOrigin,protein)]
countsPerProtein[, protein.f := factor(protein)]

nbProjPerProtein <- summary(countsPerProtein[, protein.f], maxsum = 1000000)
id_protein <- 1:NROW(nbProjPerProtein)
dt <- data.table(id_protein)
dt[, protein := labels(nbProjPerProtein)]
dt[, nbProjProt := -nbProjPerProtein]
# This ordering makes most common peptides appear first
setkey(dt, nbProjProt)
# We save the rank
dt[, rank_protein := 1:NROW(nbProjPerProtein)]
# Number of projects is brought back to a positive number
dt[, nbProjProt := -nbProjProt]

