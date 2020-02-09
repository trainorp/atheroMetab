############ Prereqs ############
## Start always run:
options(stringsAsFactors = FALSE, scipen = 900)
library(tidyverse)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive", "~/gdrive")

setwd(paste0(baseDir, "/AthroMetab/WCMC"))

# Import WCMC data:
df1 <- readxl::read_xlsx("Data/targetedData_20181203.xlsx")
df1$m55 <- NULL
metabKey <- readxl::read_xlsx("Data/metabKey_20181123.xlsx")
metabKey <- metabKey[!metabKey$metabID == "m55",]

# Import cohort phenotypes:
oxPLDF <- read.csv(paste0(baseDir, "/Athro/oxPL6/wide_data_20150529.csv"))
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
qcMeans <- qcDFL %>% group_by(metabID) %>% summarize(meanConc = mean(Concentration))
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

############ PCA ############
pcaDF <- as.data.frame(df2b)
rownames(pcaDF) <- paste0(gsub(" ", "", pcaDF$group), "_", pcaDF$timept, "_", pcaDF$ptid)
pcaDF <- pcaDF[,sapply(names(pcaDF), function(x) grepl("m\\d", x))]
pca1 <- prcomp(pcaDF, center = TRUE, scale. = TRUE)

pcaDF2 <- as.data.frame(pca1$x)
pcaDF2$samp <- rownames(pcaDF2)
pcaDF2$group <- str_split(pcaDF2$samp, "_", simplify = TRUE)[,1]
pcaDF2$timepoint <- str_split(pcaDF2$samp, "_", simplify = TRUE)[,2]
pcaDF2$groupTime <- paste0(pcaDF2$group, pcaDF2$timepoint)
pcaDF2$groupTime <- factor(pcaDF2$groupTime, levels = c("ThromboticMIT0", "Non-ThromboticMIT0", "sCADT0",
                                                        "ThromboticMITFU", "Non-ThromboticMITFU", "sCADTFU"))
ggplot(pcaDF2 %>% filter(group != "Indeterminate"), aes(x = PC1, y = PC2, color = groupTime)) + 
  geom_point() + theme_bw() + scale_color_manual(values = c("red", "blue", "black", 
                                                            rgb(255, 229, 229, 255, maxColorValue = 255), rgb(229, 229, 255, 255, maxColorValue = 255), "grey80"))

############ Summary statistics ############
summaryFun2 <- function(metab){
  tab1 <- df1b %>% select(group, timept, x = metab) %>% 
    group_by(group, timept) %>% summarize(mean = mean(x, na.rm = TRUE),sd = sd(x, na.rm = TRUE),
                                          min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),median = median(x, na.rm = TRUE),
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
rm(summaryDFT0, summaryDFTFU)

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
#write.csv(summaryDF,file="Results/MetabSummaryDF.csv",row.names=FALSE)

############ Correlation between metabolites ############
# Make this separate for 
metabCor <- cor(df2b[,names(df2b) %in% metabKey$metabID], method = "spearman")
rownames(metabCor) <- colnames(metabCor) <- metabKey$Metabolite
# write.csv(metabCor, file = "Results/metabCor.csv", row.names = FALSE)

# Correlation plot:
# png(file = "Plots/absCor.png", height = 5, width = 5, units = "in", res = 300)
corrplot::corrplot(metabCor, order = "hclust", tl.cex = .4)
# dev.off()

############ Linear model ############
metabKey$metabID

t1 <- df2b %>% filter(!group =="Indeterminate") %>% select(ptid, group, timept, m10) %>% 
  spread(key = "timept", value = "m10")
t1$d <- t1$T0 - t1$TFU
t1$group <- factor(t1$group, levels = c("Thrombotic MI", "Non-Thrombotic MI", "sCAD"))

t2 <- df2b %>% filter(!group =="Indeterminate")
t2$group <- factor(t2$group, levels = c("Thrombotic MI", "Non-Thrombotic MI", "sCAD"))

lm0 <- lm(T0 ~ TFU + group, data = t1)
lm1 <- lm(T0 ~ TFU * group, data = t1)
lm2 <- lm(d ~ group, data = t1)
lm3 <- lm(m10 ~ group * timept, data = t2)

emLm0 <- emmeans::emmeans(lm0, ~ group)
pairs(emLm0, adjust = "none")

emLm1 <- emmeans::emmeans(lm1, ~ group)
pairs(emLm1, adjust = "none")

emLm2 <- emmeans::emmeans(lm2, ~ group)
pairs(emLm2, adjust = "none")

# write.csv(t1, file = "t1.csv", row.names = FALSE, na = "")

############ Change Score Bayesian Linear model ############
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
    beta[j] ~ dnorm(0, 0.000001)
  }
  
  # Prior for the variance:
  invVar ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(invVar)
}"

# Prepare data:
df2DfromB <- df2b %>% select(-samp, -coarseGroup, -oldMIGroup, -oldGroup) %>% 
  gather(key = "metabID", value = "value", -(ptid:group))
df2DfromB <- df2DfromB %>% spread(key = timept, value = value)
df2DfromB$d <- df2DfromB$T0 - df2DfromB$TFU
df2DfromB <- cbind(df2DfromB, psych::dummy.code(df2DfromB$group))

# Run the sampler for each metabolite
samp3 <- list()
for(metab in metabKey$metabID){
  df2DfromBTemp <- df2DfromB %>% filter(metabID==metab & !is.na(TFU) & !is.na(T0))
  model <- rjags::jags.model(textConnection(rJAGSModel),
                             data = list(Y = df2DfromBTemp$d, ThrombMI = df2DfromBTemp$`Thrombotic MI`,
                                         NonThrombMI = df2DfromBTemp$`Non-Thrombotic MI`,
                                         Indeterminate = df2DfromBTemp$Indeterminate, n = nrow(df2DfromBTemp)),
                             n.chains = 5)
  
  update(model, 10000) # Burnin for 10000 samples
  samp <- rjags::coda.samples(model, variable.names = c("beta","sigma"), n.iter = 20000, thin = 5)
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
  
  data.frame(mean = mean(x), median = median(x), lq = as.numeric(quants[1]),
             uq = as.numeric(quants[2]), lq2 = as.numeric(quants2[1]), uq2 = as.numeric(quants2[2]),
             lq3 = as.numeric(quants3[1]), uq3 = as.numeric(quants3[2]))
}
samp4Sum <- samp4 %>% group_by(metab, effect) %>% do(ciFun(.$value))

# Join metabolite names:
sampSum <- samp4Sum
sampSum$Est <- paste0(sapply(samp4Sum$median, FUN = function(x) format(x, digits = 2, nsmall = 2))," (",
                      sapply(samp4Sum$lq2, FUN = function(x) format(x, digits = 2, nsmall = 2)), ", ",
                      sapply(samp4Sum$uq2, FUN = function(x) format(x, digits = 2, nsmall = 2)), ")")
sampSum$Est2 <- paste0(sapply(samp4Sum$median, FUN = function(x) format(x, digits = 2, nsmall = 2))," (",
                       sapply(samp4Sum$lq3, FUN = function(x) format(x, digits = 2, nsmall = 2)), ", ",
                       sapply(samp4Sum$uq3, FUN = function(x) format(x, digits = 2, nsmall = 2)), ")")
sampSum <- sampSum %>% select(metabID = metab, effect, Est, Est2)
sampSum <- metabKey %>% select(metabID, Metabolite, `Full Name, synonym`) %>% 
  right_join(sampSum)

# Make nice:
sampSum <- sampSum %>% select(-Est) %>% spread(key = effect, value = Est2)
sampSum <- sampSum %>% select(metabID, Metabolite, Name = `Full Name, synonym`, sCAD, T2,
                              Ind, T1, T1vssCAD, T1vsT2, T1vsInd)
# write.csv(sampSum, file = paste0("Results/changeModelSum", "_", gsub("-", "", Sys.Date()), ".csv"),
#                                  row.names = FALSE)

############ Change Score Bayesian Linear model Figure ############
samp4Sum <- metabKey %>% select(metabID, Metabolite, Name = `Full Name, synonym`) %>% 
  right_join(samp4Sum, by = c("metabID" = "metab"))
metabOrder <- samp4Sum %>% filter(effect == "T1") %>% arrange(mean) %>% select(Metabolite) %>% as.data.frame()
samp4Sum$Metabolite <- factor(samp4Sum$Metabolite, levels = metabOrder[, 1])
samp4Sum$effect[samp4Sum$effect == "T1"] <- "Thrombotic MI"
samp4Sum$effect[samp4Sum$effect == "T2"] <- "Non-Thrombotic MI"
samp4Sum$effect <- factor(samp4Sum$effect)
samp4Sum$effect <- factor(samp4Sum$effect, levels = levels(samp4Sum$effect)[c(7, 1:6)])

# png(file = paste0("Plots/changeModelSum", "_", gsub("-", "", Sys.Date()), ".png"),
#     height = 5*1.5, width = 7*1.5, units = "in", res = 300)
ggplot(data = subset(samp4Sum, effect == "Non-Thrombotic MI"), 
       aes(x = Metabolite, y = mean, ymin = lq3, ymax = uq3)) +
  geom_point(aes(color = effect)) +
  geom_errorbar(aes(ymin = lq3, ymax = uq3, col = effect), width = .6, cex = .75) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(data = subset(samp4Sum, effect == "Thrombotic MI"), aes(color = effect)) +
  geom_errorbar(data = subset(samp4Sum, effect == "Thrombotic MI"), 
                aes(ymin = lq3, ymax = uq3, col = effect), width = 0.6, cex = .75) + 
  scale_color_manual(values = c("#377EB8", "#E41A1C")) + 
  coord_flip() + theme_bw() + labs(y = "Relative Change from Quiescent to Event", color = "MI Type")
# dev.off()

save.image("working_20200205.RData")
load("working_20200205.RData")

############ T0 Variable importance ############
df2T0 <- as.data.frame(df2[df2$timept == "T0" & df2$group %in% c("sCAD", "Thrombotic MI", "Non-Thrombotic MI"), ])
df2T0 <- df2T0[, names(df2T0) %in% c("ptid", "group") | grepl("m\\d", names(df2T0))]
rownames(df2T0) <- df2T0$ptid
df2T0$ptid <- NULL
df2T0$group <- factor(df2T0$group)
set.seed(33)
rF1 <- randomForest::randomForest(group ~ ., data = df2T0, ntree = 10000, importance = TRUE)
rF1Imp <- as.data.frame(randomForest::importance(rF1))
rF1Imp$metabID <- rownames(rF1Imp)
rF1Imp <- metabKey %>% select(metabID, Metabolite, `Full Name, synonym`) %>% 
  right_join(rF1Imp) %>% as.data.frame()

rF1ImpDel <- rF1Imp[order(rF1Imp$MeanDecreaseAccuracy, decreasing = TRUE),][21:nrow(rF1Imp), "metabID"]
df2T0b <- df2T0[, !names(df2T0) %in% rF1ImpDel]
set.seed(333)
rF2 <- randomForest::randomForest(group ~ ., data = df2T0b, ntree = 1000, importance = TRUE)

############ T0 Bayesian model selection ############
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

df2T0$group <- factor(df2T0$group, levels = c("Thrombotic MI", "Non-Thrombotic MI", "sCAD"))
y <- as.numeric(as.factor(df2T0$group))
nGrps <- length(unique(y))
X <- scale(df2T0[,grepl("m\\d", names(df2T0))])
p <- dim(X)[2]
n <- dim(X)[1]

# Wrapper function for parallel:
parWrapper<-function(seedIter){
  # Model definition:
  model <- rjags::jags.model(file = textConnection(rJAGSModel3),
                             inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seedIter),
                             data = list(y = y, X = X, p = p, n = n, nGrps = nGrps), n.chains = 1, n.adapt = 2000)
  codaSamples<-tryCatch({
    codaSamples <- rjags::coda.samples(model,
                                       variable.names = c("logdens", "prob", "nVarsInc", "delta", "tau", "SD", "beta0", "beta"),
                                       n.iter = 20000, thin = 20)
  },error=function(e){
    codaSamples <- e
    return(codaSamples)
  }
  )
  codaSamples[[1]]
}

nChains <- 7
cl <- parallel::makeCluster(4)
parallel::clusterExport(cl, list("y", "X", "p", "n", "nGrps", "rJAGSModel3"))

# Do the distributed MCMC sampling:
ptm <- proc.time()
samp <- parallel::clusterApply(cl, 1:nChains, parWrapper)
proc.time() - ptm
parallel::stopCluster(cl)

# Make chain matricies:
# samp<-samp[-5]
# nChains<-nChains-1
# samps<-lapply(samp,function(x) x[[1]])
samps <- samp
chainMatrix1 <- as.matrix(samps[[1]])

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
chainDF$metabID <- paste0("m", chainDF$j)

# Filter out un-neaded due to level coding:
chainDF <- chainDF[!(chainDF$paramType == "beta" & chainDF$i == "1"),]
chainDF$group <- "Non-Thrombotic MI"
chainDF$group[chainDF$i == 3] <- "sCAD"

# Add annotation data:
chainDF <- chainDF %>% left_join(metabKey)

# Generate parameter summary:
chainParamSum <- chainDF %>% group_by(parameter, paramType, metabID, Metabolite, group) %>% 
  summarize(mean = mean(value))
cpsTemp <- chainParamSum %>% filter(paramType == "delta") %>% arrange(desc(mean))
cpsTemp2 <- cpsTemp %>% arrange(mean)
cpsTemp2$Metabolite <- factor(cpsTemp2$Metabolite, levels = cpsTemp2$Metabolite)

png(file = paste0("Plots/SVSSPosteriorMean", "_", gsub("-", "", Sys.Date()), ".png"),
                  height = 7, width = 8, res = 300, units = "in")
ggplot(cpsTemp2, aes(x = Metabolite, y = mean))+
  geom_point()+ ylab(expression(paste("Posterior Mean of ", delta))) + theme_bw() + coord_flip() 
dev.off()

# Delta parameter only:
chainDFDelta <- chainDF[chainDF$paramType == "delta", c("Metabolite", "value", "chain", "iter")]
chainDFDelta$chainIter <- paste0(chainDFDelta$chain, "_", chainDFDelta$iter)

deltaMatrix <- matrix(NA, nrow = length(unique(chainDFDelta$Metabolite)), 
                      ncol = length(unique(chainDFDelta$Metabolite)))
rownames(deltaMatrix) <- colnames(deltaMatrix) <- unique(chainDFDelta$Metabolite)
# Calculate Jaccard similarity:
for(i in 1:nrow(deltaMatrix)){
  for(j in 1:ncol(deltaMatrix)){
    tDF1 <- data.frame(x = chainDFDelta$value[chainDFDelta$Metabolite == rownames(deltaMatrix)[i]], 
                       y = chainDFDelta$value[chainDFDelta$Metabolite == colnames(deltaMatrix)[j]])
    deltaMatrix[i,j] <- jaccard::jaccard(tDF1$x, tDF1$y)
  }
}
hclustDeltaMatrix <- hclust(as.dist(as.matrix(deltaMatrix)))
cutree(hclustDeltaMatrix, h = .28)
deltaMatrix2 <- deltaMatrix[hclustDeltaMatrix$order, hclustDeltaMatrix$order]
diag(deltaMatrix2) <- NA
image(deltaMatrix2)
heatmap(deltaMatrix2)

save.image("working_20200205b.RData")
load("working_20200205b.RData")

############ T0 model ############
metabInclude <- cpsTemp$metabID[1:15]
rJAGSModel2 <- "
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i] <- exp(beta0[r] + sum(beta[r, 1:p] * X[i, 1:p]))
    }
    # Likelihood:
    y[i] ~ dcat(pi[1:nGrps, i])
    logdensi[i] <- logdensity.cat(y[i], pi[1:nGrps, i])
  }
  logdens <- sum(logdensi)
  
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, 0.01)
  }
  for(j in 1:p){
    for(r in 2:nGrps){
      beta[r,j] ~ dnorm(0, tau)
    }
  }
  tau ~ dgamma(2,1)
  SD <- sqrt(1 / tau)
}"

runjags::load.runjagsmodule()
rJAGSModel2b <- "
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i] <- exp(beta0[r] + sum(beta[r, 1:p] * X[i, 1:p]))
    }
    # Likelihood:
    y[i] ~ dcat(pi[1:nGrps, i])
    logdensi[i] <- logdensity.cat(y[i], pi[1:nGrps, i])
  }
  logdens <- sum(logdensi)
  
  # Priors:
  beta0[1] <- 0
  for(j in 1:p){
    beta[1,j] <- 0
  }
  for(r in 2:nGrps){
    beta0[r] ~ dnorm(0, 0.01)
  }
  for(j in 1:p){
    for(r in 2:nGrps){
      lambda[r,j] ~ dhalfcauchy(1)
      tau[r,j] ~ dhalfcauchy(lambda[r,j])
      beta[r,j] ~ dnorm(0, 1 / (tau[r,j]^2))
    }
  }
}"

# Sample for cross-validation
X <- scale(df2T0[,names(df2T0) %in% c(metabInclude)])
p <- dim(X)[2]
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
ptm <- proc.time()
codaSamples <- foreach(i=1:nrow(X), .inorder=FALSE) %dopar% {
  X2 <- X[-i,]
  y2 <- y[-i]
  n <- dim(X2)[1]
  set.seed(33333)
  model <- rjags::jags.model(file = textConnection(rJAGSModel2b),
            data = list(y = y2, X = X2, p = p, n = n, nGrps = nGrps), n.chains = 6, n.adapt = 1000)
  
  codaSamples <- rjags::coda.samples(model,
               variable.names = c("logdens", "tau", "SD", "beta0", "beta"), n.iter = 10000, thin = 10)
  
  # Make into one MCMC chain:
  codaSamples<-as.data.frame(do.call("rbind",codaSamples))
  codaSamples$lo<-i
  codaSamples
}
proc.time()-ptm
stopCluster(cl)
save.image('working_20200205c.RData')

load('working_20200205c.RData')
looPreds <- data.frame()
for(j in 1:length(codaSamples)){
  loo <- codaSamples[[j]]
  looPi<- data.frame(iter = 1:nrow(loo), T2 = NA, sCAD = NA)
  for(i in 1:nrow(loo)){
    looBeta <- as.matrix(loo[i, grepl("beta\\[", colnames(loo))])
    looBeta <- matrix(c(looBeta), ncol = 3, byrow = TRUE)
    looBeta2 <- t(looBeta) %*% t(t(X[j,]))
    T2 <- exp(looBeta2[2])
    sCAD <- exp(looBeta2[3])
    denom <- 1 + T2 + sCAD
    looPi$T2[i] <- T2 / denom
    looPi$sCAD[i] <- sCAD / denom
  }
  looPi$T1 <- 1 - (looPi$T2 + looPi$sCAD)
  looProp <- prop.table(table(apply(looPi[, -1], 1, which.max)))
  nam1 <- c("Non-Thrombotic MI", "sCAD", "Thrombotic MI")[match(names(looProp), c("1", "2", "3"))]
  looProp <- as.numeric(looProp)
  names(looProp) <- nam1
  looPred <- as.data.frame(t(looProp))
  looPred$ptid <- rownames(X)[j]
  looPreds <- bind_rows(looPreds, looPred)
  print(j)
}
df2T0$ptid <- rownames(df2T0)
looPreds <- looPreds %>% left_join(df2T0 %>% select(ptid, group))
looPreds$pred <- colnames(looPreds)[apply(looPreds[,1:3], 1, which.max)]
xtabs(~group + pred, data = looPreds)

model <- rjags::jags.model(file = textConnection(rJAGSModel2),
                         data = list(y = y, X = X, p = p, n = n, nGrps = nGrps), n.chains = 6, n.adapt = 1000)

codaSamples <- rjags::coda.samples(model,
                    variable.names = c("logdens", "tau", "SD", "beta0", "beta"), n.iter = 10000, thin = 10)

tropT0<-oxPLDF %>% select(ptid,tropT0) %>% mutate(logTrop=log(tropT0+.0001)) %>% select(-tropT0)
df2bT0<-df2bT0 %>% left_join(tropT0)
X<-scale(df2bT0[,names(df2bT0) %in% c(metabInclude,"logTrop")])
# X<-scale(df2bT0[,names(df2bT0) %in% metabInclude])
cat("p is ",p,"\n")