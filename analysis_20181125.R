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
phenoDF<-oxPLDF %>% select(ptid,group=migroup2,coarseGroup=group)
phenoDF$group[is.na(phenoDF$group)]<-"sCAD"
phenoDF$group[phenoDF$group=="Type 1"]<-"Thrombotic MI"
phenoDF$group[phenoDF$group=="Type 2"]<-"Non-Thrombotic MI"
phenoDF$ptid<-as.character(phenoDF$ptid)

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
df1<-df1 %>% select(-samp)

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

############ Da ############

