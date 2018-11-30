############ Prereqs ############
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")

# Import WCMC data:
df1<-readxl::read_xlsx("targetedData_20181123.xlsx")
metabKey<-readxl::read_xlsx("metabKey_20181123.xlsx")

# Import cohort phenotypes:
oxPLDF<-read.csv("~/gdrive/Athro/oxPL6/wide_data_20150529.csv")
phenoDF<-oxPLDF %>% select(ptid,group=migroup2,coarseGroup=group,oldMIGroup=migroup)
phenoDF$group[is.na(phenoDF$group)]<-"sCAD"
phenoDF$group[phenoDF$group=="Type 1"]<-"Thrombotic MI"
phenoDF$group[phenoDF$group=="Type 2"]<-"Non-Thrombotic MI"
phenoDF$ptid<-as.character(phenoDF$ptid)

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

png(file="qcPlot1.png",height=4,width=7,units="in",res=300)
ggplot(qcPlotDF,aes(x=samp,y=log(Concentration),group=Metabolite,color=Metabolite)) + 
  geom_point() + geom_line() + theme_bw() + xlab("Sample")
dev.off()

############ Summary statistics ############
summaryFun<-function(x){
  m1<-c(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),min=min(x,na.rm=TRUE),
        max=max(x,na.rm=TRUE),median=median(x,na.rm=TRUE),IQR=IQR(x,na.rm=TRUE))
  # Need to change after imputation
  return(m1)
}
summaryFun2<-function(metab){
  tab1<-as.data.frame(do.call("rbind",by(c(df1[,metab]),df1[,"group"],summaryFun)))
  tab1$metabID<-metab
  tab1$group<-rownames(tab1)
  return(tab1)
}
summaryDF<-do.call("rbind",lapply(metabKey$metabID,summaryFun2))

############ T0 Analysis ############
ggplot(df1 %>% filter(group %in% c("Thrombotic MI","Non-Thrombotic MI","sCAD")),
       aes(x=timept,y=log(m34),group=ptid,color=group)) + 
  geom_line() + geom_point() + theme_bw()

df2<-df1 %>% filter(group %in% c("Thrombotic MI","Non-Thrombotic MI","sCAD"))
form1<-as.formula(paste0("group~",paste(metabKey$metabID,collapse="+")))
df2$group<-factor(df2$group)
randomForest::randomForest(form1,data=df2 %>% filter(timept=="T0"),na.action=na.omit)

df3<-df1 %>% filter(oldGroup %in% c("Thrombotic MI","Non-Thrombotic MI","sCAD"))
form2<-as.formula(paste0("oldGroup~",paste(metabKey$metabID,collapse="+")))
df3$oldGroup<-factor(df3$oldGroup)
randomForest::randomForest(form2,data=df3 %>% filter(timept=="T0"),na.action=na.omit)
