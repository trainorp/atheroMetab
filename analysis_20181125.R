############ Prereqs ############
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")

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

png(file="Plots/qcPlot1.png",height=4,width=7,units="in",res=300)
ggplot(qcPlotDF,aes(x=samp,y=log(Concentration),group=Metabolite,color=Metabolite)) + 
  geom_point() + geom_line() + theme_bw() + xlab("Sample")
dev.off()

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
write.csv(summaryDF,file="Results/MetabSummaryDF.csv",row.names=FALSE)

############ Change from baseline linear model ############
rJagsModel<-"model{
  # Likelihood:
  for(i in 1:n){
    Y[i]~dnorm(mu[i],invVar)
    mu[i]<-beta[1]+beta[2]*Baseline[i]+beta[3]*ThrombMI[i]
  }
  
  # Prior for beta:
  beta[1]~dnorm(2794,0.0001)
  for(j in 2:3){
    beta[j]~dnorm(0,0.0001)
  }
  
  # Prior for the variance:
  invVar~dgamma(0.01,0.01)
  sigma<-1/sqrt(invVar)
}"
dc<-psych::dummy.code(df1b$group)
dc<-dc[,4]

model<-rjags::jags.model(textConnection(rJagsModel),
                         data=list(Y=df1b$m1,ThrombMI=dc,n=nrow(df1b)))

update(model,10000); # Burnin for 10000 samples
samp<-rjags::coda.samples(model, variable.names=c("beta","sigma"),n.iter=20000)

summary(samp)

############ BEST ############
priors<-list(muM=0,muSD=2)
na.omit(df2$m1[df1$group=="Thrombotic MI"])
bestM1<-BEST::BESTmcmc(na.omit(df2$m1[df1$group=="Thrombotic MI"]),
                    na.omit(df2$m1[df1$group=="Non-Thrombotic MI"]),
                    priors=NULL, parallel=FALSE)

# BANOVA
coefs2<-list()
for(var1 in metabKey$metabID){
  df2Temp<-as.data.frame(df2[!(is.na(df2$group)),
                             names(df2) %in% c(var1,"group","ptid")])
  df2Temp$group<-factor(df2Temp$group,
    levels=c("Thrombotic MI","sCAD","Non-Thrombotic MI","Indeterminate"))
  df2$ptid<-factor(df2$ptid)
  form1<-as.formula(paste0(var1,"~1"))
  banova1<-BANOVA::BANOVA.Normal(form1,~group,data=df2Temp,id=df2Temp$ptid)
  coefs1<-banova1$coef.tables$coeff_table
  coefs1$var<-var1
  coefs2[[var1]]<-coefs1
}
coefs2<-do.call("rbind",coefs2)
coefs2$group<-str_split(rownames(coefs2),"  :  ",simplify=TRUE)[,2]
rownames(coefs2)<-NULL

############ T0 Analysis ############
ggplot(df2 %>% filter(group %in% c("Thrombotic MI","Non-Thrombotic MI","sCAD")),
       aes(x=timept,y=log(m34),group=ptid,color=group)) + 
  geom_line() + geom_point() + theme_bw()

# Random Forests
df2c<-df2b %>% filter(group %in% c("Thrombotic MI","Non-Thrombotic MI","sCAD"))
df2c$group<-factor(df2c$group)
df2c<-oxPLDF %>% select(ptid,tropT0) %>% full_join(df2c)
form0<-as.formula(paste0("group~",paste(metabKey$metabID,collapse="+")))
form1<-as.formula(paste0("group~",paste(metabKey$metabID,collapse="+"),"+tropT0"))
randomForest::randomForest(form0,data=df2c %>% filter(timept=="T0"),ntree=10000,
                           strata=df2c$group,sampsize=15)
randomForest::randomForest(form0,data=df2c %>% filter(timept=="T0"),ntree=10000)
randomForest::randomForest(form1,data=df2c %>% filter(timept=="T0"),ntree=10000,
                           strata=df2c$group,sampsize=11)

