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
png(file="Plots/rel2AbsCor.png",height=67,width=12,res=100,units="in")
gridExtra::grid.arrange(grobs=plotList,ncol=3)
dev.off()

############ Correlation between metabolites ############
metabCor<-cor(df2b[,names(df2b) %in% metabKey$metabID],method="spearman")
rownames(metabCor)<-colnames(metabCor)<-metabKey$Metabolite
write.csv(metabCor,file="Results/metabCor.csv",row.names=FALSE)

# Correlation plot:
png(file="Plots/absCor.png",height=5,width=5,units="in",res=300)
corrplot::corrplot(metabCor,order="hclust",tl.cex=.4)
dev.off()

############ Change score linear model ############
# Prepare data:
df2DfromB<-df2b %>% select(-samp,-coarseGroup,-oldMIGroup,-oldGroup) %>% 
  gather(key="metabID",value="value",-(ptid:group))
df2DfromB<-df2DfromB %>% spread(key=timept,value=value)
df2DfromB$d<-df2DfromB$T0-df2DfromB$TFU
df2DfromB<-cbind(df2DfromB,psych::dummy.code(df2DfromB$group))

# JAGS model:
rJAGSModel<-"model{
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
ptm<-proc.time()
samp4Sum<-samp4 %>% group_by(metab,effect) %>% do(ciFun(.$value))
proc.time()-ptm

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
write.csv(sampSum,file="Results/changeModelSum.csv",row.names=FALSE)

############ T0 Bayesian model ############
rJAGSModel2<-"
data{
  for(j in 1:p){
    meanX[j]<-mean(x[,j])
    sdX[j]<-sd(x[,j])
    for(i in 1:n){
      zx[i,j]<-(x[i,j]-meanX[j])/sdX[j]
    }
  }
}
model{
  for(i in 1:n){
    y[i]~dcat(explambda[1:nGrps,i])
    for(r in 1:nGrps){
      explambda[r,i]<-exp(zbeta0[r]+sum(zbeta[r,1:p]*zx[i,1:p]))
    }
  }
  zbeta0[1]<-0
  for(j in 1:p){
    zbeta[1,j]<-0
  }
  for(r in 2:nGrps){
    zbeta0[r]~dnorm(0,0.0001)
    for(j in 1:p){
      zbeta[r,j]~dnorm(0,.0001)
    }
  }
  for(r in 1:nGrps){
    beta[r,1:p]<-zbeta[r,1:p]/sdX[1:p]
    beta0[r]<-zbeta0[r]-sum(zbeta[r,1:p]*meanX[1:p]/sdX[1:p])
  }
}"
y<-as.numeric(as.factor(df2b$group))
x<-df2b[,names(df2b)%in%c("m1","m3","m12","m27")]
p<-dim(x)[2]
n<-dim(x)[1]
nGrps<-length(unique(y))
model<-rjags::jags.model(file=textConnection(rJAGSModel2),
                         data=list(y=y,x=x,p=p,n=n,nGrps=nGrps),
                         n.chains=4)

update(model,10000); # Burnin for 10000 samples
samp<-rjags::coda.samples(model,variable.names=c("beta0" ,  "beta" ,  
                                                 "zbeta0" , "zbeta"),n.iter=20000)

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

