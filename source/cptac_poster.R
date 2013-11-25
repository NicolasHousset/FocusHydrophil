## Code to generate the plots associated with the BBC poster

projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"

plotPath <- "/plot/MappingPredGrad"
load(file= paste0(projectPath, "/data/CPTAC_predicted.RData"))

# First step
# Predicting on rep2, get the 3 quartiles of observation and prediction for each lab, fit the linear model
# Apply the linear model to rep1, rep2 and rep3

globalResults[, q25 := quantile(rtsec, probs = 0.25), by = list(lab, rep)]
globalResults[, q50 := quantile(rtsec, probs = 0.5), by = list(lab, rep)]
globalResults[, q75 := quantile(rtsec, probs = 0.75), by = list(lab, rep)]

setkey(globalResults, lab, rep)
lm_part1 <- globalResults[, list(lab, rep, q25, q50, q75)]
setkey(lm_part1, lab, rep)
lm_part1 <- unique(lm_part1)

lm_part1 <- data.table(melt(data = lm_part1, 
                            id.vars = c("lab", "rep"), 
                            measure.vars = c("q25","q50","q75"),
                            variable.name = "quartile",
                            value.name = "rtquartile"))

globalResults[, q25 := quantile(ELUDE_RT_NOLTS, probs = 0.25), by = list(lab, rep)]
globalResults[, q50 := quantile(ELUDE_RT_NOLTS, probs = 0.5), by = list(lab, rep)]
globalResults[, q75 := quantile(ELUDE_RT_NOLTS, probs = 0.75), by = list(lab, rep)]

setkey(globalResults, lab, rep)
lm_part2 <- globalResults[, list(lab, rep, q25, q50, q75)]
setkey(lm_part2, lab, rep)
lm_part2 <- unique(lm_part2)

lm_part2 <- data.table(melt(data = lm_part2, 
                            id.vars = c("lab", "rep"), 
                            measure.vars = c("q25","q50","q75"),
                            variable.name = "quartile",
                            value.name = "predquartile"))

setkey(lm_part1, lab, rep, quartile)
setkey(lm_part2, lab, rep, quartile)

lm_melted <- lm_part2[lm_part1]
setkey(lm_melted, lab, rep)

list_mapping <- unique(lm_melted[, list(lab, rep)])
setkey(list_mapping, lab, rep)
# Fitting a linear model of quartile(obs)~quartile(pred) for each lab.rep combination 
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    linear_model <- lm(rtquartile~predquartile, data=lm_melted[list(i,j)])
    list_mapping[list(i,j), intercept := linear_model[[1]][[1]]]
    list_mapping[list(i,j), slope := linear_model[[1]][[2]]]    
  }
}

# Here I decide on which replicate to base the calibration
# By changing which replicate is used for calibration, maybe it could be considered as cross-validation ?
setkey(globalResults, lab, rep)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    globalResults[list(i,j), rt_mapped := list_mapping[list(i,"rep2")][, intercept] + list_mapping[list(i,"rep2")][, slope] * ELUDE_RT_NOLTS]
  }
}


# Second step
# Count the number of (non modified based sequence) peptides per protein (per lab per rep)
# Filter out proteins non identified by at least 2 unique peptides

bbc_data <- globalResults[index_rt2==1][ups==FALSE]
setkey(bbc_data, lab, rep, peptide)
bbc_data <- unique(bbc_data)
# protein has been converted to factor, necessary to first convert it back to character
bbc_data[, str_protein := as.character(protein)]
list_protein_global <- data.table(NULL)
setkey(bbc_data, lab, rep)

for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    list_protein <- data.frame(summary(bbc_data[list(i,j)][, factor(str_protein)], maxsum=10000))
    names(list_protein)[c(1)] <- c("count_unique_pep")
    list_protein <- data.table(list_protein, keep.rownames = TRUE)
    list_protein <- list_protein[count_unique_pep > 1]
    list_protein[, lab := i]
    list_protein[, rep := j]
    list_protein_global <- rbind(list_protein_global, list_protein)
  }
}


setkey(bbc_data, lab, rep, protein)
setkey(list_protein_global, lab, rep, rn)

bbc_data <- bbc_data[list_protein_global]
# Third step
# Count the number of protein identifications (replicate based) per lab
# Filter out proteins non identified in 2 replicates out of 3

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
  list_protein_2 <- list_protein_2[count_rep_id > 1]
  list_protein_2[, lab := i]
  list_protein_2_global <- rbind(list_protein_2_global, list_protein_2)
}
setkey(bbc_data, lab, protein)
setkey(list_protein_2_global, lab, rn)
test <- bbc_data[list_protein_2_global]


# Fourth step
# For each lab/rep combination, count the number of proteins identified at different gradient
# To determine the min and the max, go for mapping (-1.2) and mapping (1.2) ?
# In any case divide the interval in 100 bins of x seconds