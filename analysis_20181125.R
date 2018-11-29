############ Prereqs ############
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")
df1<-readxl::read_xlsx("targetedData_20181123.xlsx")
metabKey<-readxl::read_xlsx("metabKey_20181123.xlsx")

############ Data processing ############
# QC data:
df1$samp<-gsub("Biorec_preDeFilippis0","QC",df1$samp)
df1$samp<-gsub("Biorec_postDeFilippis0","QC",df1$samp)
df1qc<-df1[grepl("QC",df1$samp),]
df1qcL<-df1qc %>% gather(key="metabID",value="Concentration",-samp)

# Blank data:
df1$samp<-gsub("MtdBlank_preDeFilippis0","Blank",df1$samp)
df1$samp<-gsub("MtdBlank_postDeFilippis0","Blank",df1$samp)

# Phenotype data:
df1$ptid<-gsub("Blank[[:digit:]]","",gsub("QC[[:digit:]]","",
                                          str_split(df1$samp,"-",simplify=TRUE)[,1]))
df1$timept<-str_split(df1$samp,"-",simplify=TRUE)[,2]
phenoDF<-df1 %>% select(samp,ptid,timept)

############ QC samples throughout run ############
qcMeans<-df1qcL %>% group_by(metabID) %>% summarize(meanConc=mean(Concentration))
qcMeansECDF<-ecdf(qcMeans$meanConc)
qcMeans$prob<-qcMeansECDF(qcMeans$meanConc)

qcPlotDF<-df1qcL %>% filter(metabID %in% c("m22","m40","m37","m34","m10","m43","m5",
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
