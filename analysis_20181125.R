############ Prereqs ############
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/gdrive/AthroMetab/WCMC")
df1<-readxl::read_xlsx("targetedData_20181123.xlsx")
