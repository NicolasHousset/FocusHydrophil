projectPath <- "C:/Users/Nicolas Housset/Documents/R_Projects/FocusHydrophil"
plotPath <- "/plot/MappingPredGrad"
load(file= paste0(projectPath, "/data/CPTAC_predicted.RData"))

setkey(result_filtered, exp)


# Trying to map the predicted RT and what we observe
graph <- ggplot(result_filtered["lab1rep1"], aes(Predicted_RT)) + geom_histogram(aes(y = ..density..), binwidth = 0.04)

png(filename = paste0(projectPath,plotPath,"/1_GradPred.png"),
                      width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()

graph <- ggplot(result_filtered["lab1rep1"], aes(rtsec)) + geom_histogram(aes(y = ..density..), binwidth = 100)
png(filename = paste0(projectPath,plotPath,"/2_GradExp.png"),
    width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()


graph <- ggplot(result_filtered["lab1rep1"], aes(Predicted_RT, rtsec)) + geom_point(alpha = 1/4) + xlim(-1.5,1.5) + ylim(0,6000)
png(filename = paste0(projectPath,plotPath,"/3_RawPredExp.png"),
    width = 800, height = 800, units = "px")
graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())
dev.off()


graph + theme(text = element_text(size = 30), panel.background = element_blank(), panel.grid = element_blank())





graph <- ggplot(result_filtered, aes(Predicted_RT)) + geom_density(aes(colour = lab)) + xlim(-1.5,1.5) + facet_grid(. ~ rep)
graph



# Let's try to automate the process
# A starting point, using the three quartiles

result_filtered[, rt25 := quantile(rtsec, probs = 0.25), by = list(lab, rep)]
result_filtered[, rt50 := quantile(rtsec, probs = 0.5), by = list(lab, rep)]
result_filtered[, rt75 := quantile(rtsec, probs = 0.75), by = list(lab, rep)]
result_filtered[, pred25 := quantile(Predicted_RT, probs = 0.25), by = list(lab, rep)]
result_filtered[, pred50 := quantile(Predicted_RT, probs = 0.5), by = list(lab, rep)]
result_filtered[, pred75 := quantile(Predicted_RT, probs = 0.75), by = list(lab, rep)]

result_filtered[, q25 := quantile(rtsec, probs = 0.25), by = list(lab, rep)]
result_filtered[, q50 := quantile(rtsec, probs = 0.5), by = list(lab, rep)]
result_filtered[, q75 := quantile(rtsec, probs = 0.75), by = list(lab, rep)]

setkey(result_filtered, lab, rep)
lm_part1 <- result_filtered[, list(lab, rep, q25, q50, q75)]
setkey(lm_part1, lab, rep)
lm_part1 <- unique(lm_part1)

lm_part1 <- data.table(melt(data = lm_part1, 
                            id.vars = c("lab", "rep"), 
                            measure.vars = c("q25","q50","q75"),
                            variable.name = "quartile",
                            value.name = "rtquartile"))


result_filtered[, q25 := quantile(Predicted_RT, probs = 0.25), by = list(lab, rep)]
result_filtered[, q50 := quantile(Predicted_RT, probs = 0.5), by = list(lab, rep)]
result_filtered[, q75 := quantile(Predicted_RT, probs = 0.75), by = list(lab, rep)]

setkey(result_filtered, lab, rep)
lm_part2 <- result_filtered[, list(lab, rep, q25, q50, q75)]
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
list_mapping[list(i,j), intercept := linear_model[[1]][[1]]]
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    linear_model <- lm(predquartile~rtquartile, data=lm_melted[list(i,j)])
    list_mapping[list(i,j), intercept := linear_model[[1]][[1]]]
    list_mapping[list(i,j), slope := linear_model[[1]][[2]]]    
  }
}

setkey(result_filtered, lab, rep)
for(i in c("lab1","lab2","lab3")){
  for(j in c("rep1","rep2","rep3")){
    result_filtered[list(i,j), rt_mapped := list_mapping[list(i,j)][, intercept] + list_mapping[list(i,j)][, slope] * rtsec]
  }
}

## Example of what the previous for is computing
## setkey(result_filtered, lab)
## result_filtered["lab1", rt_mapped := (-1.17 + rtsec * 0.0004194)]
## result_filtered["lab2", rt_mapped := (-1.2 + rtsec * 0.0004)]
## result_filtered["lab3", rt_mapped := (-1.7 + rtsec * 0.00045)]


result_filtered[, diff := Predicted_RT - rt_mapped]

result_filtered[, centile := ceiling(rt_mapped * 50)]
result_filtered[, meanError := mean(diff), by = list(lab,rep,centile)]
result_filtered[, q975 := quantile(diff, probs = 0.975), by = list(lab,rep,centile)]
result_filtered[, q025 := quantile(diff, probs = 0.025), by = list(lab,rep,centile)]
result_filtered[, q500 := quantile(diff, probs = 0.500), by = list(lab,rep,centile)]
result_filtered[, q250 := quantile(diff, probs = 0.250), by = list(lab,rep,centile)]
result_filtered[, q750 := quantile(diff, probs = 0.750), by = list(lab,rep,centile)]

ggplot(result_filtered, aes(centile, q975-q025)) + geom_point() + xlim(-75,75)+ ylim(0,1) + facet_grid(lab ~ rep)

graphDS <- data.table(melt(data = result_filtered, 
                           id.vars = c("lab", "rep", "Peptide","index_rt2", "centile"), 
                           measure.vars = c("q025","q250","q500","q750","q975"),
                           variable.name = "quantile",
                           value.name = "error"))

setkey(graphDS, lab, rep, quantile, centile)
graphDS <- unique(graphDS)

ggplot(graphDS, aes(centile, error, colour = quantile)) + geom_point(alpha=(1)) + xlim(-75,75)+ ylim(-0.6,0.3) + facet_grid(lab ~ rep)

ggplot(result_filtered, aes(centile, q500)) + geom_point(alpha=(1/2), aes(colour=rep)) + xlim(-150,150) + ylim(-1, 0.5)+ facet_grid(lab ~ rep) 