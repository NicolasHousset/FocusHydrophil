
inputPath <- "//uv2522.ugent.be/compomics/Andrea/CPTAC/"
outputPath <- "//uv2522.ugent.be/compomics/Nicolas/CPTAC/"

# Test Command
# fileName <- inputList[[2]]

joinFiles <- function(fileName,...){
  inputFile1 <- paste0(fileName, ".mgf")
  
  
  mgfDT <- fread(input = paste0(inputPath, inputFile1), sep="\n", header = FALSE)
  
  
  rtDT <- data.table(mgfDT[grepl("TITLE=", V1)][, sub("TITLE=", "", V1)])
  rtDT[, rtsec := mgfDT[grepl("RTINSECONDS=", V1)][, sub("RTINSECONDS=", "",V1)]]
  rtDT[, pepmass := mgfDT[grepl("PEPMASS=", V1)][, sub("PEPMASS=", "",V1)]]
  # Some very rare cases where there is no intensity recorded
  rtDT <- rtDT[grepl(" ", pepmass)]
  rtDT[, splitLine := strsplit(pepmass, " ")]
  rtDT[, pepMZ := lapply(splitLine, "[[", 1)]
  rtDT[, MS1_Intensity := lapply(splitLine, "[[", 2)]
  rtDT[, spectrum_id := V1]
  rtDT[, V1 := NULL]
  rtDT <- rtDT[, list(spectrum_id, rtsec, pepMZ, MS1_Intensity)]
  rtDT[, rtsec := as.double(rtsec)]
  rtDT[, pepMZ := as.double(pepMZ)]
  rtDT[, MS1_Intensity := as.double(MS1_Intensity)]
  
  inputFile2 <- paste0(fileName,"_SvenSPyeast.dat.MASCOT")
  
  mascotDT <- fread(input = paste0(inputPath, inputFile2), header = TRUE, sep = "\t")
  
  mascotDT <- mascotDT[grepl("_rank1_", paste0(spectrum_id, "_"))]
  mascotDT[, spectrum_id := sub("_rank1", "", spectrum_id)]
  mascotDT[, query := as.character(query)]
  
  setkey(rtDT, spectrum_id)
  setkey(mascotDT, spectrum_id)
  
  outputDT <- rtDT[mascotDT]
  
  
  inputFile3 <- paste0(fileName, "_SvenSPyeast.dat.PERCOLATOR.csv")
  
  percolatorDT <- data.table(read.table(file = paste0(inputPath, inputFile3),
                                        header = FALSE, row.names = NULL,
                                        sep = "\t", fill = TRUE))
  
  if(NCOL(percolatorDT)==6){
    percolatorDT <- percolatorDT[, list(V1,V2,V3,V4,V5,V6)]
    # multiple proteins make new lines with a different structure
    percolatorDT[, multi_pro := !grepl("query", V1)]
    percolatorDT[1:(NROW(percolatorDT)-1), multi_pro := percolatorDT[2:NROW(percolatorDT), multi_pro]]
    percolatorDT <- percolatorDT[grepl("query", V1)]
  }else{
    percolatorDT <- percolatorDT[, list(V1,V2,V3,V4,V5,V6,V7)]
    percolatorDT <- percolatorDT[grepl("query", V1)]
    percolatorDT[, multi_pro := (V7!="")]
  }
  percolatorDT <- percolatorDT[, list(V1,V2,V3,V4,V5,V6, multi_pro)]
  percolatorDT[, query := as.character(sub("query:","",sub(";rank:1","",V1)))] # extracting the query number
  
  list_var=list("PSMId","score_percolator","q_value_percolator","pep_percolator","peptide","protein")
  for(i in 1:6){
    eval(parse(text=paste0("percolatorDT[,", list_var[i], ":= V", i, "]")))
    eval(parse(text=paste0("percolatorDT[,V", i, ":= NULL]")))
  }
  
  
  
  setkey(outputDT, query)
  setkey(percolatorDT, query)
  
  test <- percolatorDT[outputDT]
  
  outputDT <- percolatorDT[outputDT]
  outputDT[, fileOrigin := fileName]
  
  return(outputDT)
}


inputList <- list(
  "20080311_CPTAC6_07_6A005",
  "20080311_CPTAC6_10_6B019",
  "20080311_CPTAC6_13_6C012",
  "20080311_CPTAC6_16_6D014",
  "20080311_CPTAC6_19_6E010",
  "20080313_CPTAC6_07_6A005",
  "20080313_CPTAC6_10_6B019",
  "20080313_CPTAC6_13_6C012",
  "20080313_CPTAC6_16_6D014",
  "20080313_CPTAC6_19_6E010",
  "20080315_CPTAC6_07_6A005",
  "20080315_CPTAC6_10_6B019",
  "20080315_CPTAC6_13_6C012",
  "20080315_CPTAC6_16_6D014",
  "20080315_CPTAC6_19_6E010",
  "mam_042408o_CPTAC_study6_6A018",
  "mam_042408o_CPTAC_study6_6B011",
  "mam_042408o_CPTAC_study6_6C008",
  "mam_042408o_CPTAC_study6_6D004",
  "mam_042408o_CPTAC_study6_6E004",
  "mam_050108o_CPTAC_study6_6A018",
  "mam_050108o_CPTAC_study6_6B011",
  "mam_050108o_CPTAC_study6_6C008",
  "mam_050108o_CPTAC_study6_6D004",
  "mam_050108o_CPTAC_study6_6E004",
  "mam_050108o_CPTAC_study6_6A018_080504183404",
  "mam_050108o_CPTAC_study6_6B011_080504231912",
  "mam_050108o_CPTAC_study6_6C008_080505040419",
  "mam_050108o_CPTAC_study6_6D004_080505084927",
  "mam_050108o_CPTAC_study6_6E004_080505133441",
  "Orbi2_study6a_W080314_6B007_yeast_S48_ft8_pc",
  "Orbi2_study6a_W080314_6C001_yeast_S48_ft8_pc",
  "Orbi2_study6a_W080314_6D007_yeast_S48_ft8_pc",
  "Orbi2_study6a_W080314_6E008_yeast_S48_ft8_pc",
  "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_01",
  "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_01",
  "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_01",
  "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_01",
  "Orbi2_study6b_W080321_6A013_yeast_S48_ft8_pc_02",
  "Orbi2_study6b_W080321_6B007_yeast_S48_ft8_pc_02",
  "Orbi2_study6b_W080321_6C001_yeast_S48_ft8_pc_02",
  "Orbi2_study6b_W080321_6D007_yeast_S48_ft8_pc_02",
  "Orbi2_study6b_W080321_6E008_yeast_S48_ft8_pc_02"
)

result <- data.table(NULL)

for(fileName in inputList){
  print(fileName)
  result <- rbind(result, joinFiles(fileName))
}

save(list = c("result"), file=paste0(outputPath, "CPTAC_processed.RData"),
     compress = "gzip", compression_level = 1)
write.csv(result, file=paste0(outputPath, "CPTAC_processed.csv"))