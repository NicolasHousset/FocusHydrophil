library(reshape2)
library(data.table)
## Code to generate the plots associated with the BBC poster

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"

plotPath <- "/plot/MappingPredGrad_2"
load(file= paste0(projectPath, "/data/CPTAC_predicted.RData"))

# First step
# Predicting on QC, get the 3 quartiles of observation and prediction for each lab, fit the linear model
# Apply the linear model to each sample (A to E)

globalResults[, q25 := quantile(rtsec, probs = 0.25), by = list(lab, rep, sample)]
globalResults[, q50 := quantile(rtsec, probs = 0.5), by = list(lab, rep, sample)]
globalResults[, q75 := quantile(rtsec, probs = 0.75), by = list(lab, rep, sample)]

setkey(globalResults, lab, rep, sample)
lm_part1 <- globalResults[, list(lab, rep, sample, q25, q50, q75)]
setkey(lm_part1, lab, rep, sample)
lm_part1 <- unique(lm_part1)

lm_part1 <- data.table(melt(data = lm_part1, 
                            id.vars = c("lab", "rep", "sample"), 
                            measure.vars = c("q25","q50","q75"),
                            variable.name = "quartile",
                            value.name = "rtquartile"))

globalResults[, q25 := quantile(ELUDE_RT_NOLTS, probs = 0.25), by = list(lab, rep, sample)]
globalResults[, q50 := quantile(ELUDE_RT_NOLTS, probs = 0.5), by = list(lab, rep, sample)]
globalResults[, q75 := quantile(ELUDE_RT_NOLTS, probs = 0.75), by = list(lab, rep, sample)]

setkey(globalResults, lab, rep, sample)
lm_part2 <- globalResults[, list(lab, rep, sample, q25, q50, q75)]
setkey(lm_part2, lab, rep, sample)
lm_part2 <- unique(lm_part2)

lm_part2 <- data.table(melt(data = lm_part2, 
                            id.vars = c("lab", "rep", "sample"), 
                            measure.vars = c("q25","q50","q75"),
                            variable.name = "quartile",
                            value.name = "predquartile"))

setkey(lm_part1, lab, rep, sample, quartile)
setkey(lm_part2, lab, rep, sample, quartile)

lm_melted <- lm_part2[lm_part1]
setkey(lm_melted, lab, rep, sample)

list_mapping <- unique(lm_melted[, list(lab, rep, sample)])
setkey(list_mapping, lab, rep, sample)
# Fitting a linear model of quartile(obs)~quartile(pred) for each lab.rep.sample combination
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    for(k in c("1_QC2","2_A","3_B","4_C","5_D","6_E")){
      if(NROW(lm_melted[list(i,j,k)]) > 1){
        linear_model <- lm(rtquartile~predquartile, data=lm_melted[list(i,j,k)])
        list_mapping[list(i,j,k), intercept := linear_model[[1]][[1]]]
        list_mapping[list(i,j,k), slope := linear_model[[1]][[2]]]    
      }
    }
  }
}

# Here I decide on which sample to base the calibration
# By changing which sample is used for calibration, maybe it could be considered as cross-validation ?
setkey(globalResults, lab, rep)
k <- "1_QC2"
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    globalResults[list(i,j), rt_mapped := list_mapping[list(i,j,k)][, intercept] + list_mapping[list(i,j,k)][, slope] * ELUDE_RT_NOLTS]
  }
}


# Second step
# Count the number of (non modified based sequence) peptides per protein (per lab per rep per sample)
# Filter out proteins non identified by at least 2 unique peptides

bbc_data <- globalResults[index_rt2==1][ups==FALSE]
setkey(bbc_data, lab, rep, sample, protein, peptide, rtsec) # In the case of two or more identical peptides carrying different modifications, we take into account the earliest rt
setkey(bbc_data, lab, rep, sample, protein,peptide)
bbc_data <- unique(bbc_data)
# protein has been converted to factor, necessary to first convert it back to character
bbc_data[, str_protein := as.character(protein)]
list_protein_global <- data.table(NULL)
setkey(bbc_data, lab, rep, sample)

for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    for(k in c("1_QC2","2_A","3_B","4_C","5_D","6_E")){
      list_protein <- data.frame(summary(bbc_data[list(i,j,k)][, factor(str_protein)], maxsum=10000))
      names(list_protein)[c(1)] <- c("count_unique_pep")
      list_protein <- data.table(list_protein, keep.rownames = TRUE)
      if(!(list_protein[1, rn]=="NA's")){ # Sometimes the sample can be missing
        list_protein[, lab := i]
        list_protein[, rep := j]
        list_protein[, sample := k]
        list_protein <- list_protein[count_unique_pep > 1] # proteins identified by only one unique pep are excluded
        list_protein_global <- rbind(list_protein_global, list_protein)
      }
    }
  }
}

setkey(bbc_data, lab, rep, sample, protein)
setkey(list_protein_global, lab, rep, sample, rn)

bbc_data <- bbc_data[list_protein_global]

# DEPRECATED
# Third step
# Count the number of protein identifications (replicate based) per lab
# Filter out proteins non identified in 2 replicates out of 3
# This step is deprecated since I do the analysis sample by sample

setkey(bbc_data, lab, rep, protein)
bbc_data2 <- unique(bbc_data)
# Once again, converting back protein from factor to character
bbc_data2[, str_protein := as.character(protein)]
list_protein_2_global <- data.table(NULL)
setkey(bbc_data2, lab)
for(i in c("lab1","lab2","lab3")){
  list_protein_2 <- data.frame(summary(bbc_data2[i][,factor(str_protein)], maxsum = 10000))
  names(list_protein_2)[c(1)] <- c("count_rep_id")
  list_protein_2 <- data.table(list_protein_2, keep.rownames = TRUE)
  list_protein_2 <- list_protein_2[count_rep_id > 1] # proteins identified in only one replicate are excluded
  list_protein_2[, lab := i]
  list_protein_2_global <- rbind(list_protein_2_global, list_protein_2)
}
setkey(bbc_data, lab, protein)
setkey(list_protein_2_global, lab, rn)
bbc_data <- bbc_data[list_protein_2_global]


# Fourth step
# For each lab/rep combination, count the number of proteins identified at different gradient
# To determine the min and the max, go for mapping (-1.2) and mapping (1.2) ?
# In any case divide the interval in 100 bins of x seconds

# This step is hard
# I think I want to compute the rt value of the 2nd most hydrophilic peptide of the protein
# Also I want to compute the borders of the "effective LC range" for each lab (based on rep2), and deduce my gradient levels.



bbc_data[, str_key:=paste0(lab,rep,sample,protein)]
setkey(bbc_data, str_key, rtsec)
list_key <- summary(bbc_data[, factor(str_key)], maxsum = 40000)
# This takes a while, it can probably be further optimised
# convenient_vector was overkill so it's already a slight improvement
# This piece of code computes for every time a protein is identified the order in which comes each of its peptides
# Then, by zooming on the 2nd peptides, we can immediately estimate the amount of proteins identified for every level of gradient
for(i in 1:NROW(list_key)){
  print(i)
  bbc_data[labels(list_key[i]), index_protein := 1:list_key[[i]]]
}

setkey(bbc_data, lab, rep, sample, rtsec)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    for(k in c("1_QC2","2_A","3_B","4_C","5_D","6_E")){
      bbc_data[list(i,j,k), nbPepIdentified := 1:NROW(bbc_data[list(i,j,k)])]
    }
  }
}
bbc_data[, grp_rt := ceiling((rtsec+100)/800)*800+500] # It will show the middle of the retention time interval
setkey(bbc_data, lab, rep, sample, grp_rt)
bbc_data[, error_95_LTS_glob := quantile(abs(Elude_RT_LTS95 - rtsec), probs = 0.95), by = list(lab, rep, sample)]
bbc_data[, error_95_LTS_grp := quantile(abs(Elude_RT_LTS95 - rtsec), probs = 0.95), by = list(lab, rep, sample, grp_rt)]
bbc_data[, error_95_QUARTS_glob := quantile(abs(rt_mapped - rtsec), probs = 0.95), by = list(lab, rep, sample)]
bbc_data[, error_95_QUARTS_grp := quantile(abs(rt_mapped - rtsec), probs = 0.95), by = list(lab, rep, sample, grp_rt)]


bbc_data[, error_median := quantile(rt_mapped - rtsec, probs = 0.50), by = list(lab, rep, sample, grp_rt)]
bbc_data[, error_median_LTS := quantile(Elude_RT_LTS95 - rtsec, probs = 0.50), by = list(lab, rep, sample, grp_rt)]

error_graph <- unique(bbc_data)[, list(lab, rep, sample, grp_rt, error_95_LTS_glob, error_95_LTS_grp, error_95_QUARTS_glob, error_95_QUARTS_grp, error_median, error_median_LTS)]
# Removing the calibrating data
error_graph <- error_graph[sample != "1_QC2"]
# That might be my first graph
setkey(error_graph, lab, rep, sample, grp_rt)

# Ok, I want to remove the legend, and more explicit variable names
plotPath <- "/plot/MappingPredGrad_2"
png(filename = paste0(projectPath,plotPath,"/1_ErrorByRT_QUARTS.png"),
    width = 1600, height = 900, units = "px")
ggplot(error_graph, aes(grp_rt, error_95_QUARTS_grp)) + geom_boxplot(aes(group=as.character(grp_rt)), outlier.size = 0) + geom_point(aes(color = lab), size = 3) + xlim(500,6300) + ylim(0,900) + facet_grid(. ~ lab) + theme(text = element_text(size = 36),
                                                                                                                                                                                                                    panel.background = element_blank(),
                                                                                                                                                                                                                    panel.grid = element_blank(),
                                                                                                                                                                                                                    legend.position = "none",
                                                                                                                                                                                                                    axis.title.y = element_text(angle=270, size = 45),
                                                                                                                                                                                                                    axis.title.x = element_text(size = 45)) +                                                                                                                                                                                                                   
  xlab("Retention time in seconds")+
  ylab("Absolute prediction error (95%)")
dev.off()

png(filename = paste0(projectPath,plotPath,"/2_ErrorByRT_LTS.png"),
    width = 800, height = 800, units = "px")
ggplot(error_graph, aes(grp_rt, error_95_LTS_grp)) + geom_boxplot(aes(group=as.character(grp_rt)), outlier.size = 0) + geom_point(aes(color = lab)) + xlim(500,6300) + ylim(0,900) + facet_grid(. ~ lab) + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank(),axis.title.y = element_text(angle=270, size = 45),
                                                                                                                                                                                                                 axis.title.x = element_text(size = 45))+
  xlab("Retention time in seconds")+
  ylab("Number of proteins identified")
dev.off()

setkey(error_graph, lab, rep, sample)
ggplot(unique(error_graph), aes(lab, error_95_QUARTS_glob)) + geom_boxplot(aes(group=lab), outlier.size = 0) + geom_point(aes(color = lab)) + ylim(0,900) 
ggplot(unique(error_graph), aes(lab, error_95_LTS_glob)) + geom_boxplot(aes(group=lab), outlier.size = 0) + geom_point(aes(color = lab)) + ylim(0,900)

a <- data.table(melt(data = unique(error_graph), 
                id.vars = c("lab", "rep", "sample"), 
                measure.vars = c("error_95_QUARTS_glob","error_95_LTS_glob"),
                variable.name = "type",
                value.name = "error"))
View(a)

# The conclusion of this graph is obvious: LTS performs better than QUARTS.
png(filename = paste0(projectPath,plotPath,"/3_GlobError_QUARTS_VS_LTS.png"),
    width = 800, height = 800, units = "px")
ggplot(a, aes(x=type, y=error)) + geom_boxplot() + ylim(0,800) + facet_grid(. ~ lab) + theme(text = element_text(size = 14), panel.background = element_blank(), panel.grid = element_blank())
dev.off()


bbc_data_2<- bbc_data[index_protein==2]
setkey(bbc_data_2, lab, rep, sample, rtsec)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    for(k in c("1_QC2","2_A","3_B","4_C","5_D","6_E")){
      bbc_data_2[list(i,j,k), nbProtIdentified := 1:NROW(bbc_data_2[list(i,j,k)])]
    }
  }
}

ggplot(bbc_data[list("lab2")], aes(rtsec, nbPepIdentified)) + geom_point() + xlim(0,6250) + facet_grid(rep ~ sample)
png(filename = paste0(projectPath,plotPath,"/A_protsExp.png"),
    width = 800, height = 800, units = "px")
ggplot(bbc_data_2[list("lab2","rep1","1_QC2")], aes(rtsec, nbProtIdentified)) + geom_point() + xlim(900,8100) + facet_grid(rep ~ sample) + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()

# This is my second graph
png(filename = paste0(projectPath,plotPath,"/4_ProteomeCoverageByGradient.png"),
    width = 1600, height = 900, units = "px")
ggplot(bbc_data_2[list("lab2")], aes(rtsec, nbProtIdentified)) + geom_point() + xlim(900,6100) + facet_grid(rep ~ sample) + theme(text = element_text(size = 36), panel.background = element_blank(), panel.grid = element_blank(),axis.title.y = element_text(angle=270, size = 45), axis.title.x = element_text(size = 45)) +
  xlab("Retention time in seconds")+
  ylab("Number of proteins identified")
dev.off()

View(bbc_data_2)
# Nearly done, but no more energy

bbc_data_2[, q975 := quantile(abs(rt_mapped - rtsec), probs = 0.950), by = list(lab, rep, rtsec.f)]
ggplot(bbc_data_2[list("lab2","rep3")], aes(rtsec, q975)) + geom_point()