############ Prereqs ############
## Start always run:
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")
## End always run

# Import WCMC data:
df1<-readxl::read_xlsx("Data/targetedData_20181203.xlsx")
metabKey<-readxl::read_xlsx("Data/metabKey_20181123.xlsx")

# Import cohort phenotypes:
oxPLDF<-read.csv("~/gdrive/Athro/oxPL6/wide_data_20150529.csv")
oxPLDF$ptid<-as.character(oxPLDF$ptid)
phenoDF<-oxPLDF %>% select(ptid,group=migroup2,coarseGroup=group,oldMIGroup=migroup)
phenoDF$group[is.na(phenoDF$group)]<-"sCAD"
phenoDF$group[phenoDF$group=="Type 1"]<-"Thrombotic MI"
phenoDF$group[phenoDF$group=="Type 2"]<-"Non-Thrombotic MI"

phenoDF$oldGroup<-phenoDF$group
for(i in 1:nrow(phenoDF)){
  if(!is.na(phenoDF$oldMIGroup[i]) & phenoDF$group[i]=="Thrombotic MI") phenoDF$oldGroup[i]<-NA
  if(phenoDF$group[i]=="Indeterminate") phenoDF$oldGroup[i]<-NA
}

# Import Untargeted data and key:
untarKey<-read.csv("../Data/metabolite_key2.csv")
untarDF1<-read.csv("../Data/scaled.csv")

# Mapping:
idMap<-readxl::read_xlsx("Data/idMap.xlsx")

############ Data processing ############
# QC data:
df1$samp<-gsub("Biorec_preDeFilippis0","QC",df1$samp)
df1$samp<-gsub("Biorec_postDeFilippis0","QC",df1$samp)
qcDF<-df1[grepl("QC",df1$samp),]
qcDFL<-qcDF %>% gather(key="metabID",value="Concentration",-samp)

# Blank data:
df1$samp<-gsub("MtdBlank_preDeFilippis0","Blank",df1$samp)
df1$samp<-gsub("MtdBlank_postDeFilippis0","Blank",df1$samp)

# Phenotype data:
df1$ptid<-gsub("Blank[[:digit:]]","",gsub("QC[[:digit:]]","",
                                          str_split(df1$samp,"-",simplify=TRUE)[,1]))
df1$timept<-str_split(df1$samp,"-",simplify=TRUE)[,2]
phenoDF<-df1 %>% select(samp,ptid,timept) %>% left_join(phenoDF)
df1<-df1 %>% left_join(phenoDF)

# Imputation and normalization / scaling:
which(is.na(df1[,names(df1) %in% metabKey$metabID]),arr.ind=TRUE)

impFun<-function(x){
  x2<-x[x>0]
  x[x<=0]<-min(x2)/2
  return(x)
}
df2<-df1
df2[,names(df2) %in% metabKey$metabID]<-apply(df2[,names(df2) %in% metabKey$metabID],2,impFun)
df2[,names(df2) %in% metabKey$metabID]<-log2(df2[,names(df2) %in% metabKey$metabID])

# With QC and blanks removed:
df1b<-df1 %>% filter(!is.na(group))
df2b<-df2 %>% filter(!is.na(group))

############ QC samples throughout run ############
qcMeans<-qcDFL %>% group_by(metabID) %>% summarize(meanConc=mean(Concentration))
qcMeansECDF<-ecdf(qcMeans$meanConc)
qcMeans$prob<-qcMeansECDF(qcMeans$meanConc)

qcPlotDF<-qcDFL %>% filter(metabID %in% c("m22","m40","m37","m34","m10","m43","m5",
                                           "m25","m44","m50","m30"))
qcPlotDF$metabID<-factor(qcPlotDF$metabID,levels=c("m22","m40","m37","m34","m10","m43","m5",
                                                   "m25","m44","m50","m30"))
qcPlotDF<-qcPlotDF %>% arrange(metabID)
qcPlotDF$Metabolite<-metabKey$Metabolite[match(qcPlotDF$metabID,metabKey$metabID)]
qcPlotDF$Metabolite<-factor(qcPlotDF$Metabolite,levels=unique(qcPlotDF$Metabolite))

# png(file="Plots/qcPlot1.png",height=4,width=7,units="in",res=300)
ggplot(qcPlotDF,aes(x=samp,y=log(Concentration),group=Metabolite,color=Metabolite)) + 
  geom_point() + geom_line() + theme_bw() + xlab("Sample")
# dev.off()

############ Summary statistics ############
summaryFun2<-function(metab){
  tab1<-df1b %>% select(group,timept,x=metab) %>% 
    group_by(group,timept) %>% summarize(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),
      min=min(x,na.rm=TRUE),max=max(x,na.rm=TRUE),median=median(x,na.rm=TRUE),
      IQR=IQR(x,na.rm=TRUE))
  tab1$metabID<-metab
  return(tab1)
}
summaryDF<-do.call("rbind",lapply(metabKey$metabID,function(metab) summaryFun2(metab)))

# Separate summaries by timepoint: 
summaryDFT0<-summaryDF %>% filter(timept=="T0") %>% select(-timept)
summaryDFTFU<-summaryDF %>% filter(timept=="TFU") %>% select(-timept)
names(summaryDFT0)[!names(summaryDFT0) %in% c("group","metabID")]<-
  paste(names(summaryDFT0)[!names(summaryDFT0) %in% c("group","metabID")],"T0",sep="_")
names(summaryDFTFU)[!names(summaryDFTFU) %in% c("group","metabID")]<-
  paste(names(summaryDFTFU)[!names(summaryDFTFU) %in% c("group","metabID")],"TFU",sep="_")
summaryDF<-summaryDFT0 %>% left_join(summaryDFTFU)
rm(summaryDFT0,summaryDFTFU)

# Make pretty:
summaryDF$meanSD_T0<-paste0(format(summaryDF$mean_T0,digits=2,nsmall=2,trim=TRUE),
                            "+/-",format(summaryDF$sd_T0,digits=2,nsmall=2,trim=TRUE))
summaryDF$medianIQR_T0<-paste0(format(summaryDF$median_T0,digits=2,nsmall=2,trim=TRUE),
                               "+/-",format(summaryDF$IQR_T0,digits=2,nsmall=2,trim=TRUE),
                               " [",format(summaryDF$min_T0,digits=2,nsmall=2,trim=TRUE),
                               ", ",format(summaryDF$max_T0,digits=2,nsmall=2,trim=TRUE),"]")
summaryDF$meanSD_TFU<-paste0(format(summaryDF$mean_TFU,digits=2,nsmall=2,trim=TRUE),
                            "+/-",format(summaryDF$sd_TFU,digits=2,nsmall=2,trim=TRUE))
summaryDF$medianIQR_TFU<-paste0(format(summaryDF$median_TFU,digits=2,nsmall=2,trim=TRUE),
                               "+/-",format(summaryDF$IQR_TFU,digits=2,nsmall=2,trim=TRUE),
                               " [",format(summaryDF$min_TFU,digits=2,nsmall=2,trim=TRUE),
                               ", ",format(summaryDF$max_TFU,digits=2,nsmall=2,trim=TRUE),"]")
summaryDF<-summaryDF %>% select(metabID,group,meanSD_T0,meanSD_TFU,medianIQR_T0,medianIQR_TFU)
summaryDF<-metabKey %>% select(metabID,Metabolite) %>% left_join(summaryDF)
summaryDF$group<-factor(summaryDF$group,levels=c("sCAD","Non-Thrombotic MI",
                      "Indeterminate","Thrombotic MI"))
summaryDF$metabID<-factor(summaryDF$metabID,levels=metabKey$metabID)
summaryDF<-summaryDF %>% arrange(metabID,group)
# write.csv(summaryDF,file="Results/MetabSummaryDF.csv",row.names=FALSE)

############ Untargeted vs MRM ############
# Make comparison data frame from targeted data:
compDF<-df2 %>% gather(key="metabID",value="conc",-samp)
compDF$conc<-as.numeric(compDF$conc)
tempSamp<-str_split(compDF$samp,"-",simplify=TRUE)[,1:2]
tempSamp<-paste(tempSamp[,1],tempSamp[,2],sep="-")
compDF$samp<-tempSamp
rm(tempSamp)

# Make a new sample name in untargeted data for matching:
untarDF1$samp<-paste(untarDF1$ptid,gsub("-","",untarDF1$timepoint),sep="-")
compDF$relAbund<-NA
for(i in 1:nrow(compDF)){
  tempSamp<-untarDF1$samp[match(compDF$samp[i],untarDF1$samp)]
  tempMetab<-idMap$untargeted[match(compDF$metabID[i],idMap$targeted)]
  if(!is.na(tempSamp) & !is.na(tempMetab)){
    compDF$relAbund[i]<-untarDF1[untarDF1$samp==tempSamp,names(untarDF1)==tempMetab]
  }
}
compDF<-compDF[!is.na(compDF$relAbund),]

# Plot & correlation:
rDF<-data.frame(metabID=unique(compDF$metabID),r=NA)
plotList<-list()
for(i in 1:nrow(rDF)){
  rDF$r[i]<-cor(x=compDF[compDF$metabID==rDF$metabID[i],"conc"],
      y=compDF[compDF$metabID==rDF$metabID[i],"relAbund"],method="spearman")
  
  name1<-paste(rDF$metabID[i],
               ": ",metabKey$`Full Name, synonym`[metabKey$metabID==rDF$metabID[i]])
  name2<-untarKey$biochemical[match(idMap$untargeted[match(rDF$metabID[i],idMap$targeted)],
        untarKey$id)]
  
  # Label location:
  minX<-min(compDF[compDF$metabID==rDF$metabID[i],"conc"])
  maxX<-max(compDF[compDF$metabID==rDF$metabID[i],"conc"])
  xLoc<-minX+1*(maxX-minX)/9
  minY<-min(log2(compDF[compDF$metabID==rDF$metabID[i],"relAbund"]))
  maxY<-max(log2(compDF[compDF$metabID==rDF$metabID[i],"relAbund"]))
  yLoc<-minY+8*(maxY-minY)/9
  
  plotList[[i]]<-ggplot(compDF %>% filter(metabID==rDF$metabID[i]),aes(x=conc,log2(relAbund))) + 
    geom_point() + theme_bw() + xlab("Log(Concentration)") + 
    ylab("Log(scaled abund)") + ggtitle(paste(name1,name2,sep=";\n")) +
    annotate("text",x=xLoc,y=yLoc,label=paste("r = ",formatC(rDF$r[i],digits=3)))
}
# png(file="Plots/rel2AbsCor.png",height=67,width=12,res=100,units="in")
# gridExtra::grid.arrange(grobs=plotList,ncol=3)
# dev.off()

############ Correlation between metabolites ############
metabCor<-cor(df2b[,names(df2b) %in% metabKey$metabID],method="spearman")
rownames(metabCor)<-colnames(metabCor)<-metabKey$Metabolite
# write.csv(metabCor,file="Results/metabCor.csv",row.names=FALSE)

# Correlation plot:
# png(file="Plots/absCor.png",height=5,width=5,units="in",res=300)
corrplot::corrplot(metabCor,order="hclust",tl.cex=.4)
# dev.off()

############ Change score linear model ############
# Prepare data:
df2DfromB<-df2b %>% select(-samp,-coarseGroup,-oldMIGroup,-oldGroup) %>% 
  gather(key="metabID",value="value",-(ptid:group))
df2DfromB<-df2DfromB %>% spread(key=timept,value=value)
df2DfromB$d<-df2DfromB$T0-df2DfromB$TFU
df2DfromB<-cbind(df2DfromB,psych::dummy.code(df2DfromB$group))

# JAGS model:
rJAGSModel<-"
model{
  # Likelihood:
  for(i in 1:n){
    Y[i]~dnorm(mu[i],invVar)
    mu[i]<-beta[1]+beta[2]*NonThrombMI[i]+beta[3]*Indeterminate[i]+
      beta[4]*ThrombMI[i]
  }
  
  # Prior for beta:
  for(j in 1:4){
    beta[j]~dnorm(0,0.0001)
  }
  
  # Prior for the variance:
  invVar~dgamma(0.01,0.01)
  sigma<-1/sqrt(invVar)
}"

# Run the sampler for each metabolite
samp3<-list()
for(metab in metabKey$metabID){
  df2DfromBTemp<-df2DfromB %>% filter(metabID==metab & !is.na(TFU) & !is.na(T0))
  model<-rjags::jags.model(textConnection(rJAGSModel),
                           data=list(Y=df2DfromBTemp$d,ThrombMI=df2DfromBTemp$`Thrombotic MI`,
                                     NonThrombMI=df2DfromBTemp$`Non-Thrombotic MI`,
                                     Indeterminate=df2DfromBTemp$Indeterminate,n=nrow(df2DfromBTemp)))
  
  update(model,10000); # Burnin for 10000 samples
  samp<-rjags::coda.samples(model, variable.names=c("beta","sigma"),n.iter=20000)
  samp2<-data.frame(metab=metab,
                    sCAD=as.numeric(samp[[1]][,"beta[1]"]),
                    T2=as.numeric(samp[[1]][,"beta[2]"]+samp[[1]][,"beta[1]"]),
                    Ind=as.numeric(samp[[1]][,"beta[3]"]+samp[[1]][,"beta[1]"]),
                    T1=as.numeric(samp[[1]][,"beta[4]"]+samp[[1]][,"beta[1]"]),
                    T1vssCAD=as.numeric(samp[[1]][,"beta[4]"]),
                    T1vsT2=as.numeric(samp[[1]][,"beta[4]"]-samp[[1]][,"beta[2]"]),
                    T1vsInd=as.numeric(samp[[1]][,"beta[4]"]-samp[[1]][,"beta[3]"]))
  samp3[[metab]]<-samp2
}
samp3<-do.call("rbind",samp3)
rownames(samp3)<-NULL

# Wide to long:
samp4<-samp3 %>% gather(key=effect,value=value,-metab)

# Summarize:
ciFun<-function(x){
  quants<-quantile(x,probs=c(.025,.975))
  data.frame(mean=mean(x),median=median(x),lq=as.numeric(quants[1]),
             uq=as.numeric(quants[2]))
}
samp4Sum<-samp4 %>% group_by(metab,effect) %>% do(ciFun(.$value))

# Join metabolite names:
sampSum<-samp4Sum
sampSum$Est<-paste0(sapply(samp4Sum$median,FUN=function(x) format(x,digits=2,nsmall=2))," (",
      sapply(samp4Sum$lq,FUN=function(x) format(x,digits=2,nsmall=2)),", ",
      sapply(samp4Sum$uq,FUN=function(x) format(x,digits=2,nsmall=2)),")")
sampSum<-sampSum %>% select(metabID=metab,effect,Est)
sampSum<-metabKey %>% select(metabID,Metabolite,`Full Name, synonym`) %>% 
  right_join(sampSum)

# Make nice:
sampSum<-sampSum %>% spread(key=effect,value=Est)
sampSum<-sampSum %>% select(metabID,Metabolite,Name=`Full Name, synonym`,sCAD,T2,
                            Ind,T1,T1vssCAD,T1vsT2,T1vsInd)
# write.csv(sampSum,file="Results/changeModelSum.csv",row.names=FALSE)

# Temp save:
save.image("working_20190113.RData")

############ T0 Bayesian model fitting ############
load("working_20190113.RData")

rJAGSModel2<-"
data{
  for(j in 1:p){
    meanX[j]<-mean(x[,j])
    sdX[j]<-sd(x[,j])
    for(i in 1:n){
      scaledX[i,j]<-(x[i,j]-meanX[j])/sdX[j]
    }
  }
}
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      explambda[r,i]<-exp(scaledBeta0[r]+sum(scaledBeta[r,1:p]*scaledX[i,1:p]))
    }
    y[i]~dcat(explambda[1:nGrps,i])
  }
  scaledBeta0[1]<-0
  for(j in 1:p){
    scaledBeta[1,j]<-0
  }
  for(r in 2:nGrps){
    scaledBeta0[r]~dnorm(0,0.01)
    for(j in 1:p){
      scaledBeta[r,j]~dnorm(0,.01)
    }
  }
  for(r in 1:nGrps){
    beta[r,1:p]<-scaledBeta[r,1:p]/sdX[1:p]
    beta0[r]<-scaledBeta0[r]-sum(scaledBeta[r,1:p]*meanX[1:p]/sdX[1:p])
  }
}"
df2bT0<-df2b %>% filter(timept=="T0" & group!="Indeterminate")
df2bT0$group<-factor(df2bT0$group,levels=c("Thrombotic MI","Non-Thrombotic MI","sCAD"))
y<-as.numeric(as.factor(df2bT0$group))
x<-df2bT0[,names(df2bT0)%in%c("m10","m11","m12","m13","m21","m22","m26","m33","m45")]
p<-dim(x)[2]
n<-dim(x)[1]
nGrps<-length(unique(y))
model<-rjags::jags.model(file=textConnection(rJAGSModel2),
                         data=list(y=y,x=x,p=p,n=n,nGrps=nGrps),
                         n.chains=1)

update(model,10000); # Burnin for 10000 samples
samp<-rjags::coda.samples(model,
                          variable.names=c("beta0","beta","scaledBeta0","scaledBeta"),
                          n.iter=20000)

samp<-as.matrix(samp[[1]])
acf(samp[,"beta[3,1]"][1:10000])
plot(1:10000,samp[,"beta[2,1]"][1:10000],type="l")

############ T0 Bayesian model prediction ############
# Calculate group probabilities for one iteration of Gibbs sampler
groupExp<-matrix(NA,nrow=nrow(x),ncol=3)
colnames(groupExp)<-levels(as.factor(df2bT0$group))
for(g in 1:3){
  betaVars<-paste0("beta[",g,",",1:9,"]")
  sampBeta<-samp[,match(betaVars,colnames(samp))]
  sampBeta0<-samp[,match(paste0("beta0[",g,"]"),colnames(samp))]
  groupExp[,g]<-exp(sampBeta0[1] + as.matrix(x) %*% sampBeta[1,])
}
apply(groupExp,1,sum)
groupProbs<-groupExp/apply(groupExp,1,sum)
groupProbs<-cbind(ptid=df2bT0$ptid,groupProbs)

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
