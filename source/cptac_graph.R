projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
projectPath <- "/mnt/compomics/Nicolas/R_Projects/FocusHydrophil"

plotPath <- "/plot/MappingPredGrad"
load(file= paste0(projectPath, "/data/CPTAC_predicted.RData"))

setkey(globalResults, exp)


# Trying to map the predicted RT and what we observe
graph <- ggplot(globalResults["lab1rep1"], aes(ELUDE_RT_NOLTS)) + geom_histogram(aes(y = ..density..), binwidth = 0.04)

png(filename = paste0(projectPath,plotPath,"/1_GradPred.png"),
                      width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()

graph <- ggplot(globalResults["lab1rep1"], aes(rtsec)) + geom_histogram(aes(y =..density..), binwidth = 100)
png(filename = paste0(projectPath,plotPath,"/2_GradExp.png"),
    width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()


graph <- ggplot(globalResults["lab1rep1"], aes(ELUDE_RT_NOLTS, rtsec)) + geom_point(alpha = 1/4) + xlim(-1.5,1.5) + ylim(0,8000)
png(filename = paste0(projectPath,plotPath,"/3_RawPredExp.png"),
    width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()


graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())





graph <- ggplot(globalResults, aes(ELUDE_RT_NOLTS)) + geom_density(aes(colour = lab)) + xlim(-1.5,1.5) + facet_grid(. ~ rep)
graph



# Let's try to automate the process
# A starting point, using the three quartiles


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

for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    linear_model <- lm(rtquartile~predquartile, data=lm_melted[list(i,j)])
    list_mapping[list(i,j), intercept := linear_model[[1]][[1]]]
    list_mapping[list(i,j), slope := linear_model[[1]][[2]]]    
  }
}

setkey(globalResults, lab, rep)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    globalResults[list(i,j), rt_mapped := list_mapping[list(i,"rep2")][, intercept] + list_mapping[list(i,"rep2")][, slope] * ELUDE_RT_NOLTS]
  }
}

## Example of what the previous for is computing
## setkey(result_filtered, lab)
## result_filtered["lab1", rt_mapped := (-1.17 + rtsec * 0.0004194)]
## result_filtered["lab2", rt_mapped := (-1.2 + rtsec * 0.0004)]
## result_filtered["lab3", rt_mapped := (-1.7 + rtsec * 0.00045)]

globalResults[, diff := Elude_RT_LTS50 - rtsec]
setkey(globalResults, lab, rep, rtsec)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    # Who said Lisp was the parenthesis language ?
    globalResults[list(i,j), centile := ceiling((1:NROW(globalResults[list(i,j)]))*100/NROW(globalResults[list(i,j)]))]
  }
}
setkey(globalResults, lab, rep, centile)
globalResults[, q975 := quantile(diff, probs = 0.975), by = list(lab,rep,centile)]
globalResults[, q025 := quantile(diff, probs = 0.025), by = list(lab,rep,centile)]
globalResults[, q500 := quantile(diff, probs = 0.500), by = list(lab,rep,centile)]
globalResults[, q250 := quantile(diff, probs = 0.250), by = list(lab,rep,centile)]
globalResults[, q750 := quantile(diff, probs = 0.750), by = list(lab,rep,centile)]


graphDS <- data.table(melt(data = globalResults, 
                           id.vars = c("lab", "rep", "elude_sequence","index_rt2", "centile","rtsec"), 
                           measure.vars = c("q025","q250","q500","q750","q975"),
                           variable.name = "quantile",
                           value.name = "error"))

png(filename = paste0(projectPath,plotPath,"/4_LTS_Way_rtsec.png"),
    width = 800, height = 800, units = "px")
ggplot(graphDS, aes(rtsec, error, colour = quantile)) + geom_point(alpha=(1)) + xlim(0,7500)+ ylim(-1250,800) + facet_grid(lab ~ rep)
dev.off()

setkey(graphDS, lab, rep, quantile, centile)
graphDS <- unique(graphDS)

png(filename = paste0(projectPath,plotPath,"/4_LTS_Way_centile.png"),
    width = 800, height = 800, units = "px")
ggplot(graphDS, aes(centile, error, colour = quantile)) + geom_point(alpha=(1)) + xlim(0,100)+ ylim(-1250,800) + facet_grid(lab ~ rep)
dev.off()


globalResults[, diff := rt_mapped - rtsec]
setkey(globalResults, lab, rep, rtsec)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    # Who said Lisp was the parenthesis language ?
    globalResults[list(i,j), centile := ceiling((1:NROW(globalResults[list(i,j)]))*100/NROW(globalResults[list(i,j)]))]
  }
}
setkey(globalResults, lab, rep, centile)
globalResults[, q975 := quantile(diff, probs = 0.975), by = list(lab,rep,centile)]
globalResults[, q025 := quantile(diff, probs = 0.025), by = list(lab,rep,centile)]
globalResults[, q500 := quantile(diff, probs = 0.500), by = list(lab,rep,centile)]
globalResults[, q250 := quantile(diff, probs = 0.250), by = list(lab,rep,centile)]
globalResults[, q750 := quantile(diff, probs = 0.750), by = list(lab,rep,centile)]

ggplot(globalResults, aes(centile, q975)) + geom_point() + xlim(0,100) + ylim(300,1000) + facet_grid(lab ~ rep)

graphDS_2 <- data.table(melt(data = globalResults, 
                             id.vars = c("lab", "rep", "elude_sequence","index_rt2", "centile","rtsec"), 
                             measure.vars = c("q025","q250","q500","q750","q975"),
                             variable.name = "quantile",
                             value.name = "error"))

globalResults[,quantile(q500, probs = 0.500), by = list(lab,rep)]
png(filename = paste0(projectPath,plotPath,"/5_New_Way_rtsec.png"),
    width = 800, height = 800, units = "px")
ggplot(graphDS_2, aes(rtsec, error, colour = quantile)) + geom_point(alpha=(1)) + xlim(0,7500)+ ylim(-1250,800) + facet_grid(lab ~ rep)
dev.off()

setkey(graphDS_2, lab, rep, quantile, centile)
graphDS_2 <- unique(graphDS_2)

png(filename = paste0(projectPath,plotPath,"/5_New_Way_centile.png"),
    width = 800, height = 800, units = "px")
ggplot(graphDS_2, aes(centile, error, colour = quantile)) + geom_point(alpha=(1)) + xlim(0,100)+ ylim(-1250,800) + facet_grid(lab ~ rep)
dev.off()

ggplot(result_filtered, aes(centile, q500)) + geom_point(alpha=(1/2), aes(colour=rep)) + xlim(-150,150) + ylim(-1, 0.5)+ facet_grid(lab ~ rep) 