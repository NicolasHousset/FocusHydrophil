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



convenient_vector <- 1:4000
setkey(result_filtered, fileOrigin, modseq, rtsec)
# Add an index : 1 for the first time a peptide is encountered in a LC-run, 2 the second time, etc...
# convenient_vector is automatically shrinked to the appropriate size : that is very convenient :)
result_filtered[, index_rt1 := convenient_vector, by = c("fileOrigin","modseq")]
# Slightly different index : number of times the peptide is identified in the LC-run.
result_filtered[, size_rt := .N, by = c("fileOrigin", "modseq")]

result_filtered[,MS1_Intensity := -MS1_Intensity]
setkey(result_filtered, fileOrigin, modseq, MS1_Intensity)
result_filtered[, index_rt2 := convenient_vector, by = c("fileOrigin","modseq")]
result_filtered[,MS1_Intensity := -MS1_Intensity]

setkey(result_filtered, modseq)
pepUnique <- unique(result_filtered)

mod_pep <- pepUnique[, modseq]
modif <- gregexpr("[ACDEFGHIKLMNPQRSTVWY]{1}<[-]*[ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789\\.]+>", mod_pep, ignore.case = TRUE)
modifN <- gregexpr("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]+[-]", mod_pep, ignore.case = TRUE)
modifC <- gregexpr("[-][ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]+$", mod_pep, ignore.case = TRUE)

listmod <- vector("list", NROW(modif))
for (i in (1:NROW(modif))){
  listmod[[i]] <- vector("character",4)
  for (j in (1:NROW(modif[[i]]))){
    if(modif[[i]][[j]] > -1){
      listmod[[i]][[j]] <- substr(mod_pep[[i]], modif[[i]][[j]], modif[[i]][[j]]+attr(modif[[i]],"match.length")[[j]]-1)
    }
  }
}

listmodN <- vector("list", NROW(modifN))
for (i in (1:NROW(modifN))){
  listmodN[[i]] <- vector("character",1)
  for (j in (1:NROW(modifN[[i]]))){
    if(modifN[[i]][[j]] > -1){
      listmodN[[i]][[j]] <- substr(mod_pep[[i]], modifN[[i]][[j]], modifN[[i]][[j]]+attr(modifN[[i]],"match.length")[[j]]-1)
    }
  }
}

listmodC <- vector("list", NROW(modifC))
for (i in (1:NROW(modifC))){
  listmodC[[i]] <- vector("character",1)
  for (j in (1:NROW(modifC[[i]]))){
    if(modifC[[i]][[j]] > -1){
      listmodC[[i]][[j]] <- substr(mod_pep[[i]], modifC[[i]][[j]], modifC[[i]][[j]]+attr(modifC[[i]],"match.length")[[j]]-1)
    }
  }
}

pepUnique[, listMod := listmod]
pepUnique <- pepUnique[, list(modseq, listMod)]
pepUnique[, mod1 := as.character(lapply(listmod, "[[", 1))]
pepUnique[, mod2 := as.character(lapply(listmod, "[[", 2))]
pepUnique[, mod3 := as.character(lapply(listmod, "[[", 3))]
pepUnique[, mod4 := as.character(lapply(listmod, "[[", 4))]
pepUnique[, modN := as.character(lapply(listmodN, "[[", 1))]
pepUnique[, modC := as.character(lapply(listmodC, "[[", 1))]

setkey(pepUnique, modseq)
setkey(result_filtered, modseq)
result_filtered <- pepUnique[result_filtered]

# For some reason Percolator rejects all the modified Q
result_filtered <- result_filtered[(mod1=="M<15.994915>" | mod1=="" | mod1=="C<57.021464>" | mod1=="Q<-17>") &
                                     modN=="0-"]

result_filtered[substr(modseq, 1, 2) == "0-", elude_sequence := substr(modseq, 3,nchar(modseq))]
result_filtered[, elude_sequence := gsub("<","[",elude_sequence)]
result_filtered[, elude_sequence := gsub(">","]",elude_sequence)]
result_filtered[, elude_sequence := gsub("\\.",",",elude_sequence)]


# Calling ELUDE to predict a normalized retention time
setkey(result_filtered, lab)
onelab <- result_filtered["lab1"]

setkey(onelab, elude_sequence, fileOrigin)
pepUnique <- unique(onelab)[, list(elude_sequence, rtsec)]
setkey(pepUnique,elude_sequence)
pepUnique[, q50 := quantile(rtsec, probs = 0.5), by = elude_sequence]
pepCalibration <- unique(pepUnique)[, list(elude_sequence, q50)]

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
write.table(pepUnique[, list(elude_sequence, rtsec)], file=paste0(projectPath, "/data/ELUDE/CPTAC_Peptides.txt"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(pepCalibration, file=paste0(projectPath, "/data/ELUDE/CPTAC_Calibration.txt"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)


verbFlag <- "-v"
verbLevel <- "5"

trainFlag <- "-t"
calibrationData <- "/data/ELUDE/CPTAC_Calibration.txt"
calibrationData <- shQuote(paste0(projectPath, calibrationData))
testFlag <- "-e"
testData <- "/data/ELUDE/CPTAC_Peptides.txt"
testData <- shQuote(paste0(projectPath, testData))

testRTFlag <- "-g"
ltsAdjustmentFlag <- "-c"
ltsAdjustmentCoverage <- "0.5"

loadModelFlag <- "-l"
loadModel <- "/data/ELUDE/library/modelHydrophil.model"
loadModel <- shQuote(paste0(projectPath, loadModel))

savePredictFlag <- "-o"
savePredict <- "/data/ELUDE/predictionsCPTAC.out"
savePredict <- shQuote(paste0(projectPath, savePredict))

ignoreNewTestPTMFlag <- "-p"


eludeCall <- "elude"
system2(eludeCall, args = c(verbFlag, verbLevel, trainFlag, calibrationData, testFlag, testData, testRTFlag,
                            loadModelFlag, loadModel,
                            savePredictFlag, savePredict, ignoreNewTestPTMFlag),
        stdout=TRUE, stderr = TRUE)

system2(eludeCall, args = c(verbFlag, verbLevel, testFlag, testData, testRTFlag,
                             ltsAdjustmentFlag, ltsAdjustmentCoverage, 
                             loadModelFlag, loadModel,
                             savePredictFlag, savePredict, ignoreNewTestPTMFlag),
        stdout=TRUE, stderr = TRUE)



results <- data.table(read.table(file=paste0(projectPath, "/data/ELUDE/predictionsCPTAC.out"), header = TRUE, sep = "\t"))
setkey(results, Peptide, Observed_RT)
setkey(onelab, elude_sequence, rtsec)

# ggplot(results, aes(Predicted_RT)) + geom_histogram(aes(y = ..density..), binwidth = 0.05)

test <- results[onelab]

result_filtered <- results[result_filtered]

save(result_filtered, file =  paste0(projectPath,"/data/CPTAC_predicted.RData"), compression_level=1)


