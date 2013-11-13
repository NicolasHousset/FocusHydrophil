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
setkey(result, fileOrigin)

result[c("20080311_CPTAC6_04_6QC2",
         "20080311_CPTAC6_07_6A005",
         "20080311_CPTAC6_10_6B019",
         "20080311_CPTAC6_13_6C012",
         "20080311_CPTAC6_16_6D014",
         "20080311_CPTAC6_19_6E010",
         "20080311_CPTAC6_22_6QC1",
         "mam_042408o_CPTAC_study6_6QC2",
         "mam_042408o_CPTAC_study6_6A018",
         "mam_042408o_CPTAC_study6_6B011",
         "mam_042408o_CPTAC_study6_6C008",
         "mam_042408o_CPTAC_study6_6D004",
         "mam_042408o_CPTAC_study6_6E004",
         "mam_042408o_CPTAC_study6_6QC1",
         "Orbi2_study6a_W080314_6QC2_yeast_ft8_pc",
         "Orbi2_study6a_W080314_6B007_yeast_S48_ft8_pc",
         "Orbi2_study6a_W080314_6C001_yeast_S48_ft8_pc",
         "Orbi2_study6a_W080314_6D007_yeast_S48_ft8_pc",
         "Orbi2_study6a_W080314_6E008_yeast_S48_ft8_pc",
         "Orbi2_study6a_W080314_6QC1_sigma48_ft8_pc"), rep := "rep1"]

result[c("20080313_CPTAC6_04_6QC2",
         "20080313_CPTAC6_07_6A005",
         "20080313_CPTAC6_10_6B019",
         "20080313_CPTAC6_13_6C012",
         "20080313_CPTAC6_16_6D014",
         "20080313_CPTAC6_19_6E010",
         "20080313_CPTAC6_22_6QC1",
         "mam_050108o_CPTAC_study6_6QC2",
         "mam_050108o_CPTAC_study6_6A018",
         "mam_050108o_CPTAC_study6_6B011",
         "mam_050108o_CPTAC_study6_6C008",
         "mam_050108o_CPTAC_study6_6D004",
         "mam_050108o_CPTAC_study6_6E004",
         "mam_050108o_CPTAC_study6_6QC1",
         "Orbi2_study6b_W080321_6QC2_yeast_ft8_pc_01",
         "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6QC1_sigma48_ft8_pc_01"), rep := "rep2"]

result[c("20080315_CPTAC6_04_6QC2",
         "20080315_CPTAC6_07_6A005",
         "20080315_CPTAC6_10_6B019",
         "20080315_CPTAC6_13_6C012",
         "20080315_CPTAC6_16_6D014",
         "20080315_CPTAC6_19_6E010",
         "20080315_CPTAC6_22_6QC1",
         "mam_050108o_CPTAC_study6_6QC2_080504134857",
         "mam_050108o_CPTAC_study6_6A018_080504183404",
         "mam_050108o_CPTAC_study6_6B011_080504231912",
         "mam_050108o_CPTAC_study6_6C008_080505040419",
         "mam_050108o_CPTAC_study6_6D004_080505084927",
         "mam_050108o_CPTAC_study6_6E004_080505133441",
         "mam_050108o_CPTAC_study6_6QC1_080505181949",
         "Orbi2_study6b_W080321_6QC2_yeast_ft8_pc_02",
         "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_02",
         "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_02",
         "Orbi2_study6b_W080321_6C001_yeast_S48_ft8_pc_02",
         "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_02",
         "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_02",
         "Orbi2_study6b_W080321_6QC1_sigma48_ft8_pc_02"), rep := "rep3"]


# Each experiment measures one of those 7 samples : QC2, A, B, C, D, E, QC1
result[c("20080311_CPTAC6_04_6QC2",
         "20080313_CPTAC6_04_6QC2",
         "20080315_CPTAC6_04_6QC2",
         "mam_042408o_CPTAC_study6_6QC2",
         "mam_050108o_CPTAC_study6_6QC2",
         "mam_050108o_CPTAC_study6_6QC2_080504134857",
         "Orbi2_study6a_W080314_6QC2_yeast_ft8_pc",
         "Orbi2_study6b_W080321_6QC2_yeast_ft8_pc_01",
         "Orbi2_study6b_W080321_6QC2_yeast_ft8_pc_02"), sample := "1_QC2"]

result[c("20080311_CPTAC6_07_6A005",
         "20080313_CPTAC6_07_6A005",
         "20080315_CPTAC6_07_6A005",
         "mam_042408o_CPTAC_study6_6A018",
         "mam_050108o_CPTAC_study6_6A018",
         "mam_050108o_CPTAC_study6_6A018_080504183404",
         "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_02"), sample := "2_A"]

result[c("20080311_CPTAC6_10_6B019",
         "20080313_CPTAC6_10_6B019",
         "20080315_CPTAC6_10_6B019",
         "mam_042408o_CPTAC_study6_6B011",
         "mam_050108o_CPTAC_study6_6B011",
         "mam_050108o_CPTAC_study6_6B011_080504231912",
         "Orbi2_study6a_W080314_6B007_yeast_S48_ft8_pc",
         "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_02"), sample := "3_B"]

result[c("20080311_CPTAC6_13_6C012",
         "20080313_CPTAC6_13_6C012",
         "20080315_CPTAC6_13_6C012",
         "mam_042408o_CPTAC_study6_6C008",
         "mam_050108o_CPTAC_study6_6C008",
         "mam_050108o_CPTAC_study6_6C008_080505040419",
         "Orbi2_study6a_W080314_6C001_yeast_S48_ft8_pc",
         "Orbi2_study6b_W080321_6C001_yeast_S48_ft8_pc_02"), sample := "4_C"]

result[c("20080311_CPTAC6_16_6D014",
         "20080313_CPTAC6_16_6D014",
         "20080315_CPTAC6_16_6D014",
         "mam_042408o_CPTAC_study6_6D004",
         "mam_050108o_CPTAC_study6_6D004",
         "mam_050108o_CPTAC_study6_6D004_080505084927",
         "Orbi2_study6a_W080314_6D007_yeast_S48_ft8_pc",
         "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_02"), sample := "5_D"]

result[c("20080311_CPTAC6_19_6E010",
         "20080313_CPTAC6_19_6E010",
         "20080315_CPTAC6_19_6E010",
         "mam_042408o_CPTAC_study6_6E004",
         "mam_050108o_CPTAC_study6_6E004",
         "mam_050108o_CPTAC_study6_6E004_080505133441",
         "Orbi2_study6a_W080314_6E008_yeast_S48_ft8_pc",
         "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_01",
         "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_02"), sample := "6_E"]

result[c("20080311_CPTAC6_22_6QC1",
         "20080313_CPTAC6_22_6QC1",
         "20080315_CPTAC6_22_6QC1",
         "mam_042408o_CPTAC_study6_6QC1",
         "mam_050108o_CPTAC_study6_6QC1",
         "mam_050108o_CPTAC_study6_6QC1_080505181949",
         "Orbi2_study6a_W080314_6QC1_sigma48_ft8_pc",
         "Orbi2_study6b_W080321_6QC1_sigma48_ft8_pc_01",
         "Orbi2_study6b_W080321_6QC1_sigma48_ft8_pc_02"), sample := "7_QC1"]


result[, exp := factor(paste0(lab, rep))]


result[, rtsec.f := cut(rtsec, 100)]
result[, q50_Intensity := quantile(MS1_Intensity, probs=0.50), by = list(exp, sample, rtsec.f)]
result[, q_value_percolator := as.double(as.character(q_value_percolator))]

setkey(result, exp)
# Histogram of the intensity
graph <- ggplot(result[c("lab1rep1","lab1rep2","lab1rep3")], aes(rtsec, log(q50_Intensity))) + geom_point() + xlim(0,10000) +  ylim(6.5, 13) + facet_grid(sample ~ rep)
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


setkey(result, sample)


result_filtered <- result[q_value_percolator <= 0.01][multi_pro == FALSE]
result_filtered[, ups := grepl("ups",protein)]

setkey(result_filtered, sample)
table(result_filtered[, paste0(lab,sample)])


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

