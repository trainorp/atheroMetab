############ Prereqs ############
## Start always run:
options(stringsAsFactors = FALSE, scipen = 900)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")

## End always run

# Import WCMC data:
df1 <- readxl::read_xlsx("Data/targetedData_20181203.xlsx")
metabKey <- readxl::read_xlsx("Data/metabKey_20181123.xlsx")

# Import cohort phenotypes:
oxPLDF <- read.csv("~/gdrive/Athro/oxPL6/wide_data_20150529.csv")
oxPLDF$ptid <- as.character(oxPLDF$ptid)
phenoDF <- oxPLDF %>% select(ptid, group = migroup2, coarseGroup = group, oldMIGroup = migroup)
phenoDF$group[is.na(phenoDF$group)] <- "sCAD"
phenoDF$group[phenoDF$group=="Type 1"] <- "Thrombotic MI"
phenoDF$group[phenoDF$group=="Type 2"] <- "Non-Thrombotic MI"

phenoDF$oldGroup <- phenoDF$group
for(i in 1:nrow(phenoDF)){
  if(!is.na(phenoDF$oldMIGroup[i]) & phenoDF$group[i] == "Thrombotic MI") phenoDF$oldGroup[i] <- NA
  if(phenoDF$group[i] == "Indeterminate") phenoDF$oldGroup[i] <- NA
}

# Import Untargeted data and key:
untarKey <- read.csv("../Data/metabolite_key2.csv")
untarDF1 <- read.csv("../Data/scaled.csv")

# Mapping:
idMap <- readxl::read_xlsx("Data/idMap.xlsx")

############ Data processing ############
# QC data:
df1$samp <- gsub("Biorec_preDeFilippis0", "QC", df1$samp)
df1$samp <- gsub("Biorec_postDeFilippis0", "QC", df1$samp)
qcDF <- df1[grepl("QC", df1$samp),]
qcDFL <- qcDF %>% gather(key = "metabID", value = "Concentration", -samp)

# Blank data:
df1$samp <- gsub("MtdBlank_preDeFilippis0", "Blank", df1$samp)
df1$samp <- gsub("MtdBlank_postDeFilippis0", "Blank", df1$samp)

# Phenotype data:
df1$ptid <- gsub("Blank[[:digit:]]", "", gsub("QC[[:digit:]]", "",
                                          str_split(df1$samp, "-", simplify = TRUE)[,1]))
df1$timept <- str_split(df1$samp, "-", simplify = TRUE)[,2]
phenoDF <- df1 %>% select(samp, ptid, timept) %>% left_join(phenoDF)
df1 <- df1 %>% left_join(phenoDF)

# Imputation and normalization / scaling:
which(is.na(df1[,names(df1) %in% metabKey$metabID]), arr.ind = TRUE)

impFun <- function(x){
  x2 <- x[x>0]
  x[x <= 0] <- min(x2) / 2
  return(x)
}
df2 <- df1
df2[,names(df2) %in% metabKey$metabID] <- apply(df2[,names(df2) %in% metabKey$metabID], 2, impFun)
df2[,names(df2) %in% metabKey$metabID] <- log2(df2[,names(df2) %in% metabKey$metabID])

# With QC and blanks removed:
df1b <- df1 %>% filter(!is.na(group))
df2b <- df2 %>% filter(!is.na(group))

############ QC samples throughout run ############
qcMeans <- qcDFL %>% group_by(metabID) %>% 
  summarize(meanConc = mean(Concentration), sdConc = sd(Concentration), CV = sdConc / meanConc * 100)
qcMeansECDF <- ecdf(qcMeans$meanConc)
qcMeans$prob <- qcMeansECDF(qcMeans$meanConc)

qcPlotDF <- qcDFL %>% filter(metabID %in% c("m22","m40","m37","m34","m10","m43","m5",
                                           "m25","m44","m50","m30"))
qcPlotDF$metabID <- factor(qcPlotDF$metabID, levels = c("m22","m40","m37","m34","m10","m43","m5",
                                                   "m25","m44","m50","m30"))
qcPlotDF <- qcPlotDF %>% arrange(metabID)
qcPlotDF$Metabolite <- metabKey$Metabolite[match(qcPlotDF$metabID, metabKey$metabID)]
qcPlotDF$Metabolite <- factor(qcPlotDF$Metabolite, levels = unique(qcPlotDF$Metabolite))

#png(file = "Plots/qcPlot1.png", height = 4, width = 7, units = "in", res = 300)
ggplot(qcPlotDF,aes(x = samp,y = log(Concentration), group = Metabolite, color = Metabolite)) + 
  geom_point() + geom_line() + theme_bw() + xlab("Sample")
#dev.off()

# CV's:
qcMeans <- qcMeans %>% left_join(metabKey)
qcMeans <- qcMeans %>% arrange(CV)
qcMeans$Metabolite <- factor(qcMeans$Metabolite, levels = qcMeans$Metabolite)

# png(file = "Plots/qcPlot2.png", height = 5*1.5, width = 4*1.5, units = "in", res = 300)
ggplot(qcMeans, aes(x = Metabolite, y = CV)) + geom_point() + geom_hline(yintercept = 25, color = "darkred", lty = 2) +
  coord_flip() + theme_bw()
# dev.off()

# write.csv(qcMeans, file = "Results/QCCVs.csv", row.names = FALSE)

goodMetabs <- qcMeans$metabID[qcMeans$CV < 25] 

rm(qcDF, qcDFL, qcMeans, qcPlotDF, qcMeansECDF)

############ Summary statistics ############
summaryFun2 <- function(metab){
  tab1 <- df1b %>% select(group, timept, x = metab) %>% 
    group_by(group, timept) %>% summarize(mean = mean(x, na.rm = TRUE),sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
      IQR = IQR(x, na.rm = TRUE))
  tab1$metabID <- metab
  return(tab1)
}
summaryDF <- do.call("rbind", lapply(metabKey$metabID, function(metab) summaryFun2(metab)))

# Separate summaries by timepoint: 
summaryDFT0 <- summaryDF %>% filter(timept == "T0") %>% select(-timept)
summaryDFTFU <- summaryDF %>% filter(timept == "TFU") %>% select(-timept)
names(summaryDFT0)[!names(summaryDFT0) %in% c("group", "metabID")] <-
  paste(names(summaryDFT0)[!names(summaryDFT0) %in% c("group", "metabID")], "T0", sep = "_")
names(summaryDFTFU)[!names(summaryDFTFU) %in% c("group", "metabID")] <-
  paste(names(summaryDFTFU)[!names(summaryDFTFU) %in% c("group", "metabID")], "TFU", sep = "_")
summaryDF <- summaryDFT0 %>% left_join(summaryDFTFU)

# Make pretty:
summaryDF$meanSD_T0 <- paste0(format(summaryDF$mean_T0, digits = 2, nsmall = 2, trim = TRUE),
                            "+/-", format(summaryDF$sd_T0, digits = 2, nsmall = 2, trim = TRUE))
summaryDF$medianIQR_T0 <- paste0(format(summaryDF$median_T0, digits = 2, nsmall = 2, trim = TRUE),
                               "+/-", format(summaryDF$IQR_T0, digits = 2, nsmall = 2, trim = TRUE),
                               " [", format(summaryDF$min_T0, digits = 2, nsmall = 2, trim = TRUE),
                               ", ", format(summaryDF$max_T0, digits = 2, nsmall = 2, trim = TRUE), "]")
summaryDF$meanSD_TFU <- paste0(format(summaryDF$mean_TFU, digits = 2, nsmall = 2, trim = TRUE),
                            "+/-", format(summaryDF$sd_TFU, digits = 2, nsmall = 2, trim = TRUE))
summaryDF$medianIQR_TFU <- paste0(format(summaryDF$median_TFU, digits = 2, nsmall = 2, trim = TRUE),
                               "+/-", format(summaryDF$IQR_TFU, digits = 2, nsmall = 2, trim = TRUE),
                               " [", format(summaryDF$min_TFU, digits = 2, nsmall = 2, trim = TRUE),
                               ", ", format(summaryDF$max_TFU, digits = 2, nsmall = 2, trim = TRUE), "]")
summaryDF <- summaryDF %>% select(metabID, group, meanSD_T0, meanSD_TFU, medianIQR_T0, medianIQR_TFU)
summaryDF <- metabKey %>% select(metabID, Metabolite) %>% left_join(summaryDF)
summaryDF$group <- factor(summaryDF$group, levels = c("sCAD", "Non-Thrombotic MI",
                      "Indeterminate", "Thrombotic MI"))
summaryDF$metabID <- factor(summaryDF$metabID, levels = metabKey$metabID)
summaryDF <- summaryDF %>% arrange(metabID, group)
# write.csv(summaryDF, file = "Results/MetabSummaryDF.csv", row.names = FALSE)

rm(summaryDFTFU, summaryDFT0, summaryDF)

############ Untargeted vs MRM ############
# Make comparison data frame from targeted data:
compDF <- df2 %>% gather(key = "metabID", value = "conc", -samp, -(ptid:oldGroup))
compDF$conc <- as.numeric(compDF$conc)
tempSamp <- str_split(compDF$samp, "-", simplify = TRUE)[,1:2]
tempSamp <- paste(tempSamp[,1], tempSamp[,2], sep = "-")
compDF$samp <- tempSamp
rm(tempSamp)

# Make a new sample name in untargeted data for matching:
untarDF1$samp <- paste(untarDF1$ptid, gsub("-", "", untarDF1$timepoint), sep = "-")
compDF$relAbund <- NA
for(i in 1:nrow(compDF)){
  tempSamp <- untarDF1$samp[match(compDF$samp[i], untarDF1$samp)]
  tempMetab <- idMap$untargeted[match(compDF$metabID[i], idMap$targeted)]
  if(!is.na(tempSamp) & !is.na(tempMetab)){
    compDF$relAbund[i] <- untarDF1[untarDF1$samp == tempSamp, names(untarDF1) == tempMetab]
  }
}
compDF <- compDF[!is.na(compDF$relAbund),]

# Plot & correlation:
rDF <- data.frame(metabID = unique(compDF$metabID), r = NA)
plotList <- list()
for(i in 1:nrow(rDF)){
  rDF$r[i] <- cor(x = compDF[compDF$metabID == rDF$metabID[i], "conc"],
      y = compDF[compDF$metabID == rDF$metabID[i], "relAbund"], method = "spearman")
  
  name1 <- paste(rDF$metabID[i],
               ": ", metabKey$`Full Name, synonym`[metabKey$metabID == rDF$metabID[i]])
  name2 <- untarKey$biochemical[match(idMap$untargeted[match(rDF$metabID[i], idMap$targeted)],
        untarKey$id)]
  
  # Label location:
  minX <- min(compDF[compDF$metabID == rDF$metabID[i], "conc"])
  maxX <- max(compDF[compDF$metabID == rDF$metabID[i], "conc"])
  xLoc <- minX + (maxX - minX) / 9
  minY <- min(log2(compDF[compDF$metabID == rDF$metabID[i], "relAbund"]))
  maxY <- max(log2(compDF[compDF$metabID == rDF$metabID[i], "relAbund"]))
  yLoc <- minY + 8 * (maxY - minY) / 9
  
  plotList[[i]] <- ggplot(compDF %>% filter(metabID == rDF$metabID[i]), aes(x = conc, y = log2(relAbund))) + 
    geom_point() + theme_bw() + xlab("Log(Concentration)") + 
    ylab("Log(scaled abund)") + ggtitle(paste(name1, name2, sep = ";\n")) +
    annotate("text", x = xLoc, y = yLoc, label = paste("r = ", formatC(rDF$r[i], digits = 3)))
}
# png(file = "Plots/rel2AbsCor.png", height = 67, width = 12, res = 100, units = "in")
gridExtra::grid.arrange(grobs = plotList, ncol = 3)
# dev.off()

rm(plotList, compDF)

############ Correlation between metabolites ############
metabCor <- cor(df2b[, names(df2b) %in% metabKey$metabID], method = "spearman")
metabCorT0 <- cor(df2b[grepl("T0", df2b$samp), names(df2b) %in% metabKey$metabID], method = "spearman")
metabCorTFU <- cor(df2b[grepl("TFU", df2b$samp), names(df2b) %in% metabKey$metabID], method = "spearman")
rownames(metabCor) <- colnames(metabCor) <- rownames(metabCorT0) <- colnames(metabCorT0) <- rownames(metabCorTFU) <- 
  colnames(metabCorTFU) <- metabKey$Metabolite
# write.csv(metabCor, file = "Results/metabCor.csv", row.names = FALSE)

# Correlation plots:
# png(file = "Plots/absCor.png", height = 5, width = 5, units = "in", res = 300)
corrplot::corrplot(metabCor, order = "hclust", tl.cex = .4)
# dev.off()

# png(file = "Plots/absCorT0.png", height = 5, width = 5, units = "in", res = 300)
corrplot::corrplot(metabCorT0, order = "hclust", tl.cex = .4)
# dev.off()

# png(file = "Plots/absCorTFU.png", height = 5, width = 5, units = "in", res = 300)
corrplot::corrplot(metabCorTFU, order = "hclust", tl.cex = .4)
# dev.off()

############ Heatmaps ############
# Make DF's
df2bT0 <- as.data.frame(df2b[grepl("T0", df2b$samp), names(df2b) %in% metabKey$metabID])
df2bT0names <- df2b[grepl("T0", df2b$samp), match(c("group", "ptid"), names(df2b))]
rownames(df2bT0) <- paste(df2bT0names$group, df2bT0names$ptid, sep = "_")
colnames(df2bT0) <- metabKey$Metabolite[match(colnames(df2bT0), metabKey$metabID)]

df2bTFU <- as.data.frame(df2b[grepl("TFU", df2b$samp), names(df2b) %in% metabKey$metabID])
df2bTFUnames <- df2b[grepl("TFU", df2b$samp), match(c("group", "ptid"), names(df2b))]
rownames(df2bTFU) <- paste(df2bTFUnames$group, df2bTFUnames$ptid, sep = "_")
colnames(df2bTFU) <- metabKey$Metabolite[match(colnames(df2bTFU), metabKey$metabID)]

# Make plots
# png(file = "Plots/heatmapT0.png", height = 6, width = 6, units = "in", res = 300)
heatmap(scale(as.matrix(df2bT0), scale = FALSE), cexRow = .5 * (0.2 + 1/log10(nrow(df2bT0))), 
        cexCol = .5 * (0.2 + 1/log10(ncol(df2bT0))))
# dev.off()

# png(file = "Plots/heatmapTFU.png", height = 6, width = 6, units = "in", res = 300)
heatmap(scale(as.matrix(df2bTFU), scale = FALSE), cexRow = .5 * (0.2 + 1/log10(nrow(df2bTFU))), 
        cexCol = .5 * (0.2 + 1/log10(ncol(df2bTFU))))
# dev.off()

############ Change score linear model ############
# Prepare data:
df2DfromB <- df2b %>% select(-samp, -coarseGroup, -oldMIGroup, -oldGroup) %>% 
  gather(key = "metabID", value = "value", -(ptid:group))
df2DfromB <- df2DfromB %>% spread(key = timept, value = value)
df2DfromB$d <- df2DfromB$T0 - df2DfromB$TFU
df2DfromB <- cbind(df2DfromB, psych::dummy.code(df2DfromB$group))

# JAGS model:
rJAGSModel <- "
model{
  # Likelihood:
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], invVar)
    mu[i] <- beta[1] + beta[2] * NonThrombMI[i] + beta[3] * Indeterminate[i] +
      beta[4] * ThrombMI[i]
  }
  
  # Prior for beta:
  for(j in 1:4){
    beta[j] ~ dnorm(0, 0.0001)
  }
  
  # Prior for the variance:
  invVar ~ dgamma(0.01, 0.01)
  sigma <- 1 / sqrt(invVar)
}"

# Run the sampler for each metabolite
samp3 <- list()
for(metab in metabKey$metabID){
  df2DfromBTemp <- df2DfromB %>% filter(metabID==metab & !is.na(TFU) & !is.na(T0))
  model <- rjags::jags.model(textConnection(rJAGSModel),
                           data = list(Y = df2DfromBTemp$d, ThrombMI = df2DfromBTemp$`Thrombotic MI`,
                                     NonThrombMI = df2DfromBTemp$`Non-Thrombotic MI`,
                                     Indeterminate = df2DfromBTemp$Indeterminate, n = nrow(df2DfromBTemp)))
  
  update(model, 10000) # Burnin for 10000 samples
  samp <- rjags::coda.samples(model, variable.names = c("beta","sigma"), n.iter = 20000)
  samp2 <- data.frame(metab = metab,
                    sCAD = as.numeric(samp[[1]][,"beta[1]"]),
                    T2 = as.numeric(samp[[1]][,"beta[2]"] + samp[[1]][,"beta[1]"]),
                    Ind = as.numeric(samp[[1]][,"beta[3]"] + samp[[1]][,"beta[1]"]),
                    T1 = as.numeric(samp[[1]][,"beta[4]"] + samp[[1]][,"beta[1]"]),
                    T1vssCAD = as.numeric(samp[[1]][,"beta[4]"]),
                    T1vsT2 = as.numeric(samp[[1]][,"beta[4]"] - samp[[1]][,"beta[2]"]),
                    T1vsInd = as.numeric(samp[[1]][,"beta[4]"] - samp[[1]][,"beta[3]"]))
  samp3[[metab]] <- samp2
}
samp3 <- do.call("rbind", samp3)
rownames(samp3) <- NULL

# Wide to long:
samp4 <- samp3 %>% gather(key = effect, value = value, -metab)

# Summarize:
ciFun <- function(x){
  quants <- quantile(x, probs = c(.025, .975))
  quants2 <- HDInterval::hdi(x)
  quants3 <- HDInterval::hdi(x, credMass = .9)
  
  data.frame(mode = MCMCglmm::posterior.mode(x), mean = mean(x), median = median(x), lq = as.numeric(quants[1]),
             uq = as.numeric(quants[2]), lq2 = as.numeric(quants2[1]), uq2 = as.numeric(quants2[2]),
             lq3 = as.numeric(quants3[1]), uq3 = as.numeric(quants3[2]))
}
samp4Sum <- samp4 %>% group_by(metab, effect) %>% do(ciFun(.$value))

# Join metabolite names:
sampSum <- samp4Sum
sampSum$Est <- paste0(sapply(samp4Sum$mode, FUN = function(x) format(x, digits = 2, nsmall = 2))," (",
      sapply(samp4Sum$lq2, FUN = function(x) format(x, digits = 2, nsmall = 2)), ", ",
      sapply(samp4Sum$uq2, FUN = function(x) format(x, digits = 2, nsmall = 2)), ")")
sampSum$Est2 <- paste0(sapply(samp4Sum$mode, FUN = function(x) format(x, digits = 2, nsmall = 2))," (",
                      sapply(samp4Sum$lq3, FUN = function(x) format(x, digits = 2, nsmall = 2)), ", ",
                      sapply(samp4Sum$uq3, FUN = function(x) format(x, digits = 2, nsmall = 2)), ")")
sampSum <- sampSum %>% select(metabID = metab, effect, Est, Est2)
sampSum <- metabKey %>% select(metabID, Metabolite, `Full Name, synonym`) %>% 
  right_join(sampSum)

# Make nice:
sampSum <- sampSum %>% select(-Est) %>% spread(key = effect, value = Est2)
sampSum <- sampSum %>% select(metabID, Metabolite, Name = `Full Name, synonym`, sCAD, T2,
                            Ind, T1, T1vssCAD, T1vsT2, T1vsInd)
# write.csv(sampSum, file = "Results/changeModelSum.csv", row.names = FALSE)

# Make HDI plots:
hdiPlots <- metabKey %>% select(metabID, Metabolite, `Full Name, synonym`) %>% 
  right_join(samp4Sum %>% filter(effect %in% c("T1", "Ind", "T2", "sCAD")), by = c("metabID" = "metab"))
hdiPlots <- hdiPlots %>% group_by(Metabolite) %>% mutate(meanMode = mean(mode))
hdiPlots <- hdiPlots %>% arrange(desc(meanMode))
hdiPlots$Metabolite <- factor(hdiPlots$Metabolite, levels = unique(hdiPlots$Metabolite))

hdiPlots$effect <- factor(hdiPlots$effect)
levels(hdiPlots$effect)[match(c("T1", "Ind", "T2", "sCAD"), levels(hdiPlots$effect))] <- 
  c("Thrombotic MI","Indeterminate MI", "Non-Thrombotic MI", "Stable CAD")
hdiPlots$effect <- factor(hdiPlots$effect, levels = c("Thrombotic MI","Indeterminate MI", "Non-Thrombotic MI", "Stable CAD"))

p1 <- ggplot(hdiPlots %>% filter(Metabolite %in% levels(hdiPlots$Metabolite)[1:11]), aes(x = effect, y = mode, color = effect)) + 
  geom_point() + geom_errorbar(aes(ymin = lq3, ymax = uq3), width = .2) +
  geom_hline(yintercept = 0, lty = 2) + coord_flip() + facet_wrap(~Metabolite, nrow = 1) + 
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")

p2 <- ggplot(hdiPlots %>% filter(Metabolite %in% levels(hdiPlots$Metabolite)[12:22]), aes(x = effect, y = mode, color = effect)) + 
  geom_point() + geom_errorbar(aes(ymin = lq3, ymax = uq3), width = .2) +
  geom_hline(yintercept = 0, lty = 2) + coord_flip() + facet_wrap(~Metabolite, nrow = 1) + 
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")

p3 <- ggplot(hdiPlots %>% filter(Metabolite %in% levels(hdiPlots$Metabolite)[23:33]), aes(x = effect, y = mode, color = effect)) + 
  geom_point() + geom_errorbar(aes(ymin = lq3, ymax = uq3), width = .2) +
  geom_hline(yintercept = 0, lty = 2) + coord_flip() + facet_wrap(~Metabolite, nrow = 1) + 
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")

p4 <- ggplot(hdiPlots %>% filter(Metabolite %in% levels(hdiPlots$Metabolite)[34:44]), aes(x = effect, y = mode, color = effect)) + 
  geom_point() + geom_errorbar(aes(ymin = lq3, ymax = uq3), width = .2) +
  geom_hline(yintercept = 0, lty = 2) + coord_flip() + facet_wrap(~Metabolite, nrow = 1) + 
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")

p5 <- ggplot(hdiPlots %>% filter(Metabolite %in% levels(hdiPlots$Metabolite)[45:55]), aes(x = effect, y = mode, color = effect)) + 
  geom_point() + geom_errorbar(aes(ymin = lq3, ymax = uq3), width = .2) +
  geom_hline(yintercept = 0, lty = 2) + coord_flip() + facet_wrap(~Metabolite, nrow = 1) + 
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")

# png(file = "Plots/changeModelSum.png", height = 8 * 1.5, width = 14 * 1.5, units = "in", res = 300)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, layout_matrix = matrix(c(1:5), nrow = 5))
# dev.off()

# Temp save:
rm(samp, samp2, samp3, samp4, samp4Sum, sampSum, model, p1, p2, p3, p4, p5, hdiPlots)
# save.image("working_20190113.RData")

############ T0 Bayesian model selection ############
load("working_20190113.RData")

rJAGSModel3 <- "
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i] <- exp(beta0[r] + sum(beta[r,1:p] * X[i,1:p]))
    }
    # Likelihood:
    y[i] ~ dcat(pi[1:nGrps,i])
    logdensi[i] <- logdensity.cat(y[i], pi[1:nGrps,i])
  }
  logdens <- sum(logdensi)
  
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, tau2)
  }
  for(j in 1:p){
    delta[1,j] ~ dbern(prob)
    for(r in 2:nGrps){
      gamma[r,j] ~ dnorm(0, tau)
      beta[r,j] <- gamma[r,j] * delta[1,j]
    }
  }
  nVarsInc <- sum(delta[1,1:p])
  prob ~ dbeta(5, 100)
  tau ~ dgamma(2, 2)
  tau2 ~ dgamma(2, 2)
  SD <- sqrt(1 / tau)
}"

df2bT0 <- df2b %>% filter(timept=="T0" & group != "Indeterminate")
df2bT0$group <- factor(df2bT0$group, levels = c("Thrombotic MI", "Non-Thrombotic MI", "sCAD"))
y <- as.numeric(as.factor(df2bT0$group))
nGrps <- length(unique(y))
X <- scale(df2bT0[,grepl("m\\d", names(df2bT0))])
X <- X[,colnames(X) %in% goodMetabs]
p <- dim(X)[2]
n <- dim(X)[1]

# Wrapper function for parallel:
parWrapper<-function(seedIter){
  # Model definition:
  model <- rjags::jags.model(file = textConnection(rJAGSModel3),
             inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seedIter),
             data = list(y = y, X = X, p = p, n = n, nGrps = nGrps), n.chains = 1, n.adapt = 1000)
  codaSamples<-tryCatch({
      codaSamples <- rjags::coda.samples(model,
           variable.names = c("logdens", "prob", "nVarsInc", "delta", "tau", "SD", "beta0", "beta"),
           n.iter = 100000, thin = 10)
    },error=function(e){
      codaSamples <- e
      return(codaSamples)
    }
  )
  codaSamples[[1]]
}

# Create the cluster and export needed variables/data:
nChains <- 8
cl <- parallel::makeCluster(8)
parallel::clusterExport(cl, list("y", "X", "p", "n", "nGrps", "rJAGSModel3"))

# Do the distributed MCMC sampling:
ptm <- proc.time()
samps <- parallel::clusterApply(cl, 1:nChains, parWrapper)
proc.time() - ptm
parallel::stopCluster(cl)

# Save result:
# save.image("working_20190121.RData")
load("working_20190121.RData")

# Wide to long:
chainDF <- data.frame()
for(i in 1:nChains){
  chainDFTemp <- as.data.frame(as.matrix(samps[[i]])) %>% gather(key="parameter")
  chainDFTemp$chain <- i
  chainDFTemp$iter <- 1:nrow(as.matrix(samps[[i]]))
  chainDF <- rbind(chainDF, chainDFTemp)
}
chainDF$chain <- factor(chainDF$chain)

# Higher parameters:
chainDFHigher <- data.frame()
for(i in 1:nChains){
  chainDFTemp <- as.data.frame(as.matrix(samps[[i]]))
  chainDFTemp <- chainDFTemp[,colnames(chainDFTemp) %in% c("logdens", "nVarsInc", "prob", "tau", "SD")]
  chainDFTemp$chain <- i
  chainDFHigher <- rbind(chainDFHigher, chainDFTemp)
}

# Process MCMC samples:
chainDF$paramType <- str_split(chainDF$parameter, "\\[|\\,", simplify = TRUE)[,1]
chainDF <- chainDF[!chainDF$paramType %in% c("logdens", "nVarsInc", "prob"),]
chainDF$i <- gsub("\\]", "", str_split(chainDF$parameter, "\\[|\\,", simplify = TRUE)[,2])
chainDF$j <- gsub("\\]", "", str_split(chainDF$parameter, "\\[|\\,", simplify = TRUE)[,3])
goodMetabDF <- data.frame(metabID = colnames(X), j = as.character(1:length(colnames(X))))
chainDF <- chainDF %>% left_join(goodMetabDF)

# Filter out un-neaded due to level coding:
chainDF <- chainDF[!(chainDF$paramType == "beta" & chainDF$i == "1"),]
chainDF$group <- "Non-Thrombotic MI"
chainDF$group[chainDF$i == 3] <- "sCAD"

# Add annotation data:
chainDF <- chainDF %>% left_join(metabKey)
rm(samps)

############ Variable selection analysis ############
rJAGSModel3Priors <- "
model{
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, 0.05)
  }
  for(j in 1:p){
    delta[1,j] ~ dbern(prob)
    for(r in 2:nGrps){
      gamma[r,j] ~ dnorm(0, tau)
      beta[r,j] <- gamma[r,j] * delta[1,j]
    }
  }
  nVarsInc <- sum(delta[1,1:p])
  prob ~ dbeta(5, 100)
  tau ~ dgamma(1, 1)
  SD <- sqrt(1 / tau)
}"
priorModel <- rjags::jags.model(file = textConnection(rJAGSModel3Priors),
                           inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 333),
                           data = list(p = p, nGrps = nGrps), n.chains = 10, n.adapt = 10000)
priorSamples <- rjags::coda.samples(priorModel,
  variable.names = c("prob", "nVarsInc","gamma", "delta", "tau", "SD", "beta0", "beta"), n.iter = 400000, thin = 10)
priorSamples <- do.call(rbind, priorSamples)

# Plot prior distributions:
priorProb <- data.frame(x = priorSamples[,"prob"])
priorTau <- data.frame(x = priorSamples[,"tau"])
priorTau$variance <- 1 / priorTau$x
priorTau$sd <- sqrt(priorTau$var)
priorGamma <- data.frame(x = c(priorSamples[,grepl("gamma", colnames(priorSamples))])[1:1e6])

p1 <- ggplot(priorProb, aes(x = x, y = ..density..)) + 
  geom_histogram(fill = "grey60", color = "black", bins = 50) + 
  theme_bw() + xlab("Value") + ylab("Density") + ggtitle(expression(paste("Prior for ",delta))) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(priorTau, aes(x = x, y = ..density..)) + 
  geom_histogram(fill = "grey60", color = "black", bins = 50) + 
  theme_bw() + xlab("Value") + ylab("Density") + xlim(-1, 10) + ggtitle(expression(paste("Prior for ", tau))) + 
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(priorTau, aes(x = log(variance), y = ..density..)) + 
  geom_histogram(fill = "grey60", color = "black", bins = 50) + 
  theme_bw() + xlab("Value") + ylab("Density") + ggtitle(expression(paste("Prior for ", log(sigma ^ 2)))) + 
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(priorGamma, aes(x = x, y = ..density..)) + 
  geom_histogram(fill = "grey60", color = "black", bins = 50) + 
  theme_bw() + xlab("Value") + ylab("Density") + xlim(-15, 15) + ggtitle(expression(paste("Prior for ", gamma))) + 
  theme(plot.title = element_text(hjust = 0.5))

# png(file="Plots/priorDists.png",height=6,width=7,res=300,units="in")
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# dev.off()

rm(priorModel, priorSamples, priorProb, priorTau, priorGamma, p1, p2, p3, p4)

# Plot some chains:
set.seed(3)
chainSamp <- sample(1:10, 5)
p1 <- ggplot(chainDF %>% filter(metabID == "m21" & paramType == "beta" & group == "Non-Thrombotic MI" & chain %in% chainSamp),
           aes(x = iter, y = value, color = chain)) + geom_line() + 
  theme_bw() + xlim(500, 1500) + ylim(-7, 1) + xlab("Iteration") + ylab(expression(beta))
p2 <- ggplot(chainDF %>% filter(metabID == "m21" & paramType == "beta" & group == "Non-Thrombotic MI" & chain %in% chainSamp),
           aes(x = iter, y = value, color = chain)) + geom_line() + 
  theme_bw() + xlim(500, 600) + ylim(-5, 1) + xlab("Iteration") + ylab(expression(beta))
p3 <- ggplot(chainDF %>% filter(metabID == "m21" & paramType == "beta" & group == "Non-Thrombotic MI" & chain %in% chainSamp),
           aes(x = value, fill = chain, y = ..density..)) + 
  geom_histogram(binwidth = .35, alpha = 0.55, position = "identity", color = rgb(0, 0, 0, .4)) + xlab(expression(beta)) +
  ylab("Density") + xlim(-8, 3) + theme_bw()

# png(file = "Plots/cortisolBetaChain.png", height = 6, width = 9, res = 300, units = "in")
lmat <- rbind(c(1, 1, 1, 1), c(2, 2, 3, 3))
gridExtra::grid.arrange(p1, p2, p3, layout_matrix = lmat)
# dev.off()

p1 <- ggplot(chainDF %>% filter(metabID == "m21" & paramType == "beta" & group == "Non-Thrombotic MI"),
           aes(x = value, y = ..density..)) + geom_histogram(bins = 35, color = "black", fill = "grey60") + 
   theme_bw() + ggtitle(metabKey$Metabolite[metabKey$metabID == "m21"]) + ylim(0, 1) +
  xlim(-7, 2.5) + theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(chainDF %>% filter(metabID == "m16" & paramType == "beta" & group == "Non-Thrombotic MI"),
           aes(x = value, y = ..density..)) + geom_histogram(bins = 35, color = "black", fill = "grey60") + 
  theme_bw() + ggtitle(metabKey$Metabolite[metabKey$metabID == "m16"]) + ylim(0, 3.25) +
  xlim(-7, 2.5) + theme(plot.title = element_text(hjust = 0.5))
p3 <- ggplot(chainDF %>% filter(metabID == "m32" & paramType == "beta" & group == "Non-Thrombotic MI"),
           aes(x = value, y = ..density..)) + geom_histogram(bins = 35, color = "black", fill = "grey60") + 
  theme_bw() + ggtitle(metabKey$Metabolite[metabKey$metabID == "m32"]) + ylim(0, 4.25) +
  xlim(-4, 4) + theme(plot.title = element_text(hjust = 0.5))
p4 <- ggplot(chainDF %>% filter(metabID == "m54" & paramType == "beta" & group == "Non-Thrombotic MI"),
           aes(x = value, y = ..density..)) + geom_histogram(bins = 35, color = "black", fill = "grey60") + 
  theme_bw() + ggtitle(metabKey$Metabolite[metabKey$metabID == "m54"]) + ylim(0, 7) +
  xlim(-2.5, 2.5) + theme(plot.title = element_text(hjust = 0.5))

# png(file = "Plots/betaHistograms.png", height = 7, width = 8, res = 300, units = "in")
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# dev.off()

# ACF:
acfDF <- chainDF %>% filter(metabID == "m21" & paramType == "beta" & group == "Non-Thrombotic MI" & chain == 1)
# png(file = "Plots/MCMCChainACF.png", height = 5, width = 6, res = 300, units = "in")
acf(acfDF$value, na.action = na.omit, main = "ACF for Cortisol Non-Thrombotic Effect")
# dev.off()

# Likelihood:
rf <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
hexbin::hexbinplot(logdens ~ SD, data = chainDFHigher, colramp = rf)
hexbin::hexbinplot(logdens ~ nVarsInc, data = chainDFHigher, colramp = rf)

# Generate parameter summary:
chainParamSum <- chainDF %>% group_by(parameter, paramType, metabID, Metabolite, group) %>% 
  summarize(mean = mean(value))

cpsTemp <- chainParamSum %>% filter(paramType == "delta") %>% arrange(desc(mean))
cpsTemp2 <- cpsTemp %>% arrange(mean)
cpsTemp2$Metabolite <- factor(cpsTemp2$Metabolite, levels = cpsTemp2$Metabolite)

# png(file = "Plots/SVSSPosteriorMean.png", height = 7, width = 8, res = 300, units = "in")
ggplot(cpsTemp2, aes(x = Metabolite, y = mean))+
  geom_point() + ylab(expression(paste("Posterior Mean of ", delta))) + theme_bw() + coord_flip() 
# dev.off()

rm(chainDF, chainDFHigher, chainDFTemp, p1, p2, p3, p4, acfDF)
save.image(file = "working_20190223.RData")

############ T0 Bayesian model fitting ############
# LOH 
load(file = "working_20190223.RData")

metabInclude <- cpsTemp$metabID[1:15]
rJAGSModel2 <- "
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i] <- exp(beta0[r] + sum(beta[r,1:p] * X[i,1:p]))
    }
    # Likelihood:
    y[i] ~ dcat(pi[1:nGrps,i])
    logdensi[i] <- logdensity.cat(y[i], pi[1:nGrps,i])
  }
  logdens <- sum(logdensi)
  
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, tau2)
  }
  for(j in 1:p){
    for(r in 2:nGrps){
      beta[r,j] ~ dnorm(0, tau)
    }
  }
  tau ~ dgamma(2,2)
  tau2 ~ dgamma(2,2)
  SD <- sqrt(1 / tau)
}"

tropT0 <- oxPLDF %>% select(ptid, tropT0) %>% mutate(logTrop = log(tropT0 + .0001)) %>% select(-tropT0)
df2bT0 <- df2bT0 %>% left_join(tropT0)
X <- scale(df2bT0[,names(df2bT0) %in% c(metabInclude, "logTrop")])
# X<-scale(df2bT0[,names(df2bT0) %in% metabInclude])
p <- dim(X)[2]

# Model priors:
rJAGSModel2Priors <- "
model{
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, tau2)
  }
  for(j in 1:p){
    for(r in 2:nGrps){
      beta[r,j] ~ dnorm(0, tau)
    }
  }
  tau ~ dgamma(2, 2)
  tau2 ~ dgamma(2, 2)
  SD <- sqrt(1 / tau)
}"
priorModel<-rjags::jags.model(file = textConnection(rJAGSModel2Priors),
                         data = list(p = p, nGrps = nGrps), n.chains = 6, n.adapt = 10000)
codaSamplesPriorModel <- rjags::coda.samples(priorModel,
              variable.names = c("tau", "SD", "beta0", "beta"), n.iter = 100000, thin = 10)
codaSamplesPriorModel <- as.data.frame(do.call("rbind", codaSamplesPriorModel))
codaSamplesPriorModel <- codaSamplesPriorModel[,!grepl("beta\\[1,", names(codaSamplesPriorModel))]
codaSamplesPriorModel$`beta0[1]` <- NULL

ggplot(codaSamplesPriorModel, aes(x = `beta[2,1]`)) + geom_histogram(bins = 60, color = "black", fill = "grey60") +
  theme_bw()
ggplot(codaSamplesPriorModel, aes(x = `beta0[2]`)) + geom_histogram(bins = 60, color = "black", fill = "grey60") +
  theme_bw()
ggplot(codaSamplesPriorModel, aes(x = tau)) + geom_histogram(bins = 60, color = "black", fill = "grey60") +
  theme_bw()
ggplot(codaSamplesPriorModel, aes(x = SD)) + geom_histogram(bins = 60, color = "black", fill = "grey60") +
  theme_bw()


# Sample for cross-validation
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
ptm <- proc.time()
cvList <- foreach(i = 1:nrow(X), .inorder = FALSE) %dopar% {
  X2 <- X[-i,]
  n <- dim(X2)[1]
  set.seed(33333)
  model <- rjags::jags.model(file = textConnection(rJAGSModel2),
                           data = list(y = y, X = X2, p = p, n = n, nGrps = nGrps), n.chains = 1, n.adapt = 1000)
  
  codaSamples <- rjags::coda.samples(model,
                                   variable.names = c("logdens", "tau", "SD", "beta0", "beta"), 
                                   n.iter = 1000, thin = 10)
  
  # Make into one MCMC chain:
  codaSamples <- as.data.frame(do.call("rbind", codaSamples))
  
  # Calculate group probabilities from LOO-CV posteriors
  groupProbsList <- list()
  for(j in 1:nrow(codaSamples)){
    groupExp <- matrix(0, nrow = 1, ncol = 3)
    colnames(groupExp) <- levels(as.factor(df2bT0$group))
    for(g in 2:3){
      betaVars <- paste0("beta[", g, ",", 1:p, "]")
      codaSampBeta <- codaSamples[j, match(betaVars, colnames(codaSamples))]
      codaSampBeta0 <- codaSamples[j, match(paste0("beta0[", g, "]"), colnames(codaSamples))]
      groupExp[1, g] <- exp(codaSampBeta0 + X[i,] %*% t(codaSampBeta))
    }
    groupExp[,1] <- 1
    groupProbs <- groupExp / apply(groupExp, 1, sum)
    groupProbs <- data.frame(ptid = df2bT0$ptid[i], groupProbs)
    groupProbsList[[j]] <- groupProbs
  }
  groupProbsDF <- do.call("rbind", groupProbsList)
  groupProbsDF
}
proc.time()-ptm
stopCluster(cl)

# LOH

save.image(file="working_20190223b.RData")


# Not LOO model:
n<-dim(X)[1]
set.seed(3333333)
model<-rjags::jags.model(file=textConnection(rJAGSModel2),
                         data=list(y=y,X=X,p=p,n=n,nGrps=nGrps),n.chains=6,n.adapt=1000)
codaSamplesOneModel<-rjags::coda.samples(model,
                                         variable.names=c("logdens","tau","SD","beta0","beta"),n.iter=10000,thin=10)
codaSamplesOneModel<-as.data.frame(do.call("rbind",codaSamplesOneModel))

codaSOMParamQuant<-data.frame()
for(colName in colnames(codaSamplesOneModel)){
  codaSOMParamQuant<-rbind(codaSOMParamQuant,t(as.matrix(quantile(codaSamplesOneModel[,colName],probs=seq(0,1,.1)))))
}
codaSOMParamQuant$param<-colnames(codaSamplesOneModel)

############ T0 Bayesian model prediction from LOO-CV ############
# Combind sets of chains:
codaSamples<-do.call("rbind",codaSamples)

# Calculate group probabilities from LOO-CV posteriors
groupExpList<-groupProbsList<-list()
for(i in 1:nrow(codaSamples)){
  groupExp<-matrix(0,nrow=1,ncol=3)
  colnames(groupExp)<-levels(as.factor(df2bT0$group))
  for(g in 2:3){
    betaVars<-paste0("beta[",g,",",1:p,"]")
    codaSampBeta<-codaSamples[i,match(betaVars,colnames(codaSamples))]
    codaSampBeta0<-codaSamples[i,match(paste0("beta0[",g,"]"),colnames(codaSamples))]
    groupExp[1,g]<-exp(codaSampBeta0 + X[codaSamples$lo[i],] %*% t(codaSampBeta))
  }
  groupExp[,1]<-1
  groupProbs<-groupExp/apply(groupExp,1,sum)
  groupProbs<-data.frame(ptid=df2bT0$ptid[codaSamples$lo[i]],groupProbs)
  
  # Add to lists:
  groupExp<-as.data.frame(groupExp)
  groupExp$iter<-i
  groupProbs$iter<-i
  groupExpList[[i]]<-groupExp
  groupProbsList[[i]]<-groupProbs
  if(i %% 2880==0){
    print(i/nrow(groupProbs)*100)
  }
}
groupExp<-do.call("rbind",groupExpList)
groupProbs<-do.call("rbind",groupProbsList)

# Multinomial loss:
groupProbs<-df2bT0 %>% select(ptid,group) %>% left_join(groupProbs)
groupProbs$ind<-match(gsub("-",".",gsub(" ",".",groupProbs$group)),names(groupProbs))
groupProbs$multLoss<-NA
for(i in 1:nrow(groupProbs)){
  groupProbs$multLoss[i]<-groupProbs[i,groupProbs$ind[i]]
  if(i %% 2880==0){
    print(i/nrow(groupProbs)*100)
  }
}
groupProbs$multLoss<-unlist(groupProbs$multLoss)
groupProbs$multLoss<-(-log(groupProbs$multLoss))
save.image(file="working_20190223c.RData")
load(file="working_20190223c.RData")

# CV-estimated confusion matrix:
groupProbsSum<-groupProbs %>% select(ptid,Thrombotic.MI,Non.Thrombotic.MI,sCAD) %>% group_by(ptid) %>% 
  summarize(Thrombotic.MI=median(Thrombotic.MI),Non.Thrombotic.MI=median(Non.Thrombotic.MI),
            sCAD=median(sCAD))
groupProbsSum$predGroup<-c("Thrombotic MI","Non-Thrombotic MI","sCAD")[apply(groupProbsSum[,
                                                                                           c("Thrombotic.MI","Non.Thrombotic.MI","sCAD")],1,which.max)]
groupProbsSum<-df2bT0 %>% select(ptid,group) %>% left_join(groupProbsSum,by="ptid")
xtabs(~group+predGroup,data=groupProbsSum)
prop.table(xtabs(~group+predGroup,data=groupProbsSum),2)

# Plots
groupProbsL<-groupProbs %>% select(-multLoss,-group,-ind,-predGroup) %>% 
  gather(key="Group",value="Probability",-iter,-ptid)

ggplot(groupProbsL %>% filter(ptid=="2003"),aes(x=Probability,y=..density..,fill=Group)) + 
  geom_histogram(binwidth=.04,alpha=0.5,position="identity",color="black") + theme_bw()
ggplot(groupProbsL %>% filter(ptid=="2010"),aes(x=Probability,y=..density..,fill=Group)) + 
  geom_histogram(binwidth=.04,alpha=0.5,position="identity",color="black") + theme_bw()

ggplot(groupProbsL %>% filter(ptid=="2003"),aes(x=Probability,y=..density..)) + 
  geom_histogram(binwidth=.04,alpha=0.5,position="identity",color="black") + 
  facet_grid(~Group) + theme_bw()
ggplot(groupProbsL %>% filter(ptid=="2010"),aes(x=Probability,y=..density..)) + 
  geom_histogram(binwidth=.04,alpha=0.5,position="identity",color="black") + 
  facet_grid(~Group) + theme_bw()

############ GLM-net variable selection ############
df2bT0<-oxPLDF %>% select(ptid,tropT0) %>% right_join(df2bT0)

# CV selection of alpha:
cvDF1<-expand.grid(rep=1:20,alpha=seq(0,1,.05),mis=NA,deviance=NA)
for(i in 1:nrow(cvDF1)){
  set.seed(cvDF1$rep[i])
  folds<-sample(rep(seq(10),length=nrow(df2bT0)))
  cvFit<-glmnet::cv.glmnet(x=as.matrix(df2bT0[,names(df2bT0) %in% metabKey$metabID]),
                           y=df2bT0$group,alpha=cvDF1$alpha[i],type.multinomial="grouped",
                           family="multinomial",type.measure="class",foldid=folds)
  cvDF1$mis[i]<-cvFit$cvm[cvFit$lambda==cvFit$lambda.min]
  cvFit2<-glmnet::cv.glmnet(x=as.matrix(df2bT0[,names(df2bT0) %in% metabKey$metabID]),
                           y=df2bT0$group,alpha=cvDF1$alpha[i],type.multinomial="grouped",
                           family="multinomial",foldid=folds)
  cvDF1$deviance[i]<-cvFit2$cvm[cvFit2$lambda==cvFit2$lambda.min]
  print(i)
}
cvDF1Sum<-cvDF1 %>% group_by(alpha) %>% summarize(mis=mean(mis),deviance=mean(deviance))
cvDF1Sum<-cvDF1Sum %>% mutate(misRank=rank(mis,ties.method='first'),
                         devianceRank=rank(deviance,ties.method='first'),
                         rankSum=misRank+devianceRank)

# Set alpha:
alpha<-.10
# CV selection of lambda:
cvDF2<-data.frame(rep=1:20,minLambdaMis=NA,minLambdaDeviance=NA)
for(i in 1:nrow(cvDF2)){
  set.seed(cvDF2$rep[i])
  folds<-sample(rep(seq(10),length=nrow(df2bT0)))
  cvFit<-glmnet::cv.glmnet(x=as.matrix(df2bT0[,names(df2bT0) %in% metabKey$metabID]),
                           y=df2bT0$group,alpha=alpha,type.multinomial="grouped",
                           family="multinomial",type.measure="class",foldid=folds)
  cvDF2$minLambdaMis[i]<-cvFit$lambda.min
  cvFit2<-glmnet::cv.glmnet(x=as.matrix(df2bT0[,names(df2bT0) %in% metabKey$metabID]),
                            y=df2bT0$group,alpha=alpha,type.multinomial="grouped",
                            family="multinomial",foldid=folds)
  cvDF2$minLambdaDeviance[i]<-cvFit2$lambda.min
  print(i)
}
cvDF2Sum<-cvDF2 %>% summarize(lambdaMis=mean(minLambdaMis),lambdaDeviance=mean(minLambdaDeviance))
lambda<-cvDF2Sum$lambdaMis

png(file="Plots/eNetCVMis.png",height=5,width=6,units="in",res=300)
plot(cvFit)
dev.off()

png(file="Plots/eNetCVDev.png",height=5,width=6,units="in",res=300)
plot(cvFit2)
dev.off()

# "Final model"
eNetModel<-glmnet::glmnet(x=as.matrix(df2bT0[,names(df2bT0) %in% metabKey$metabID]),
               y=df2bT0$group,alpha=alpha,type.multinomial="grouped",lambda=lambda,
               family="multinomial")
eNetModelCoef<-coef(eNetModel)
eNetCoefNames<-names(eNetModelCoef)
eNetModelCoef<-do.call("cbind",lapply(eNetModelCoef,as.matrix))
eNetModelCoef<-as.data.frame(eNetModelCoef)
names(eNetModelCoef)<-eNetCoefNames
eNetModelCoef<-eNetModelCoef[eNetModelCoef[,1]!=0,]
eNetModelCoef$metabID<-rownames(eNetModelCoef)

# Join with metabolite key:
eNetModelCoef<-as.data.frame(right_join(metabKey %>% select(Metabolite,metabID), eNetModelCoef))

# Comparison to ref:
eNetModelCoef$ThrombvsNon<-eNetModelCoef[,"Thrombotic MI"]-eNetModelCoef[,"Non-Thrombotic MI"]
eNetModelCoef$ThrombvsSCAD<-eNetModelCoef[,"Thrombotic MI"]-eNetModelCoef[,"sCAD"]
write.csv(eNetModelCoef,"Results/eNetCoefs.csv",row.names=FALSE)

# Predictions:
set.seed(333)
looList<-sample(1:nrow(df2bT0))
looENetDF<-data.frame()
for(i in 1:length(looList)){
  eNetModel<-glmnet::glmnet(x=as.matrix(df2bT0[-i,names(df2bT0) %in% metabKey$metabID]),
                y=df2bT0$group[-i],alpha=alpha,type.multinomial="grouped",lambda=lambda,
                family="multinomial")
  eNetPred<-predict(eNetModel,newx=as.matrix(df2bT0[i,names(df2bT0) %in% metabKey$metabID]),
                    type="response")[,,1]
  eNetPred<-cbind(df2bT0[i,c("ptid","group")],t(as.matrix(eNetPred)))
  looENetDF<-rbind(looENetDF,eNetPred)
}
looENetDF$pred<-c("Thrombotic MI","Non-Thrombotic MI","sCAD")[apply(looENetDF[,c("Thrombotic MI",
                                            "Non-Thrombotic MI","sCAD")],1,which.max)]
write.csv(looENetDF,"Results/looENet.csv",row.names=FALSE)

# Confusion analysis:
ENetConf<-xtabs(~group+pred,data=looENetDF)
write.csv(ENetConf,"Results/ENetConf.csv")

looENetDF$nT1<-ifelse(looENetDF$group!="Thrombotic MI","nT1","T1")
looENetDF$predNT1<-ifelse(looENetDF$pred!="Thrombotic MI","nT1","T1")
ENetConfNT1<-xtabs(~nT1+predNT1,data=looENetDF)
write.csv(ENetConfNT1,"Results/ENetConfNT1.csv")

############ RF ############
# Formula without troponin:
form0<-as.formula(paste0("group~",paste(metabKey$metabID,collapse="+")))
# Formula with troponin:
form1<-as.formula(paste0("group~",paste(metabKey$metabID,collapse="+"),"+tropT0"))

# Ensemble without troponin:
rF0<-randomForest::randomForest(form0,data=df2bT0,ntree=10000,importance=TRUE,
                           strata=df2bT0,sampsize=12)
rF0Imp<-randomForest::importance(rF0)

# Ensemble with troponin:
rF1<-randomForest::randomForest(form1,data=df2bT0,ntree=10000,importance=TRUE,
                           strata=df2bT0,sampsize=12)
rF1Imp<-randomForest::importance(rF1)
rF1Imp<-as.data.frame(rF1Imp)
rF1Imp$metabID<-rownames(rF1Imp)
rF1Imp<-rF1Imp %>% left_join(metabKey)
