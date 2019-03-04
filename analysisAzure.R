############ Prereqs ############
## Start always run:
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

load(file="working_20190223.RData")

metabInclude<-cpsTemp$metabID[1:25]
rJAGSModel2<-"
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i]<-exp(beta0[r]+sum(beta[r,1:p]*X[i,1:p]))
    }
    # Likelihood:
    y[i]~dcat(pi[1:nGrps,i])
    logdensi[i]<-logdensity.cat(y[i],pi[1:nGrps,i])
  }
  logdens<-sum(logdensi)
  
  # Priors:
  beta0[1]<-0
  for(j in 1:p){
    beta[1,j]<-0
  }
  for(r in 2:nGrps){
    beta0[r]~dnorm(0,0.01)
  }
  for(j in 1:p){
    for(r in 2:nGrps){
      beta[r,j] ~ dnorm(0,tau)
    }
  }
  tau~dgamma(3,1)
  SD<-sqrt(1/tau)
}"

X<-scale(df2bT0[,names(df2bT0) %in% metabInclude])
p<-dim(X)[2]

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