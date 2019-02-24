############ Prereqs ############
## Start always run:
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

load(file="working_20190223.RData")
# Sample for cross-validation
library(doParallel)
cl<-makeCluster(3)
registerDoParallel(cl)
ptm<-proc.time()
codaSamples<-foreach(i=1:nrow(X),.inorder=FALSE) %dopar% {
  X2<-X[-i,]
  n<-dim(X2)[1]
  set.seed(33333)
  model<-rjags::jags.model(file=textConnection(rJAGSModel2),
                           data=list(y=y,X=X2,p=p,n=n,nGrps=nGrps),n.chains=6,n.adapt=1000)
  
  codaSamples<-rjags::coda.samples(model,
                                   variable.names=c("logdens","tau","SD","beta0","beta"),n.iter=10000,thin=10)
  
  # Make into one MCMC chain:
  codaSamples<-as.data.frame(do.call("rbind",codaSamples))
  codaSamples$lo<-i
  codaSamples
}
proc.time()-ptm
stopCluster(cl)

save.image(file="working_20190223b.RData")