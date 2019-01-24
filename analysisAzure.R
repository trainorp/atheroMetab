############ Prereqs ############
## Start always run:
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

setwd("~/atheroMetab")

## End always run

############ T0 Bayesian model selection ############
load("working_20190113.RData")
rJAGSModel3<-"
model{
  for(i in 1:n){
    for(r in 1:nGrps){
      pi[r,i]<-exp(beta0[r]+sum(beta[r,1:p]*X[i,1:p]))
    }
    # Likelihood:
    y[i]~dcat(pi[1:nGrps,i])
  }
  
  # Priors:
  beta0[1]<-0
  for(j in 1:p){
    beta[1,j]<-0
  }
  for(r in 2:nGrps){
    beta0[r]~dnorm(0,0.01)
  }
  for(j in 1:p){
    delta[1,j] ~ dbern(prob)
    for(r in 2:nGrps){
      gamma[r,j] ~ dnorm(0,tau)
      beta[r,j] <- gamma[r,j]*delta[1,j]
    }
  }
  nVarsInc<-sum(delta[1,1:p])
  prob~dbeta(5,100)
  tau~dgamma(1,1)
  SD<-sqrt(1/tau)
}"

df2bT0<-df2b %>% filter(timept=="T0" & group!="Indeterminate")
df2bT0$group<-factor(df2bT0$group,levels=c("Thrombotic MI","Non-Thrombotic MI","sCAD"))
y<-as.numeric(as.factor(df2bT0$group))
nGrps<-length(unique(y))
X<-scale(df2bT0[,grepl("m\\d",names(df2bT0))])
p<-dim(X)[2]
n<-dim(X)[1]

# Wrapper function for parallel:
parWrapper<-function(seedIter){
  # Model definition:
  model<-rjags::jags.model(file=textConnection(rJAGSModel3),
                           inits=list(.RNG.name="base::Wichmann-Hill",.RNG.seed=seedIter),
                           data=list(y=y,X=X,p=p,n=n,nGrps=nGrps),n.chains=1,n.adapt=1000)
  tryCatch({
    codaSamples<-rjags::coda.samples(model,
                                     variable.names=c("prob","nVarsInc","delta","tau","SD","beta0","beta"),
                                     n.iter=50000,thin=10)
  }, error=function(e){
    codaSamples<-e
  }
  )
  codaSamples
}

# Create the cluster and export needed variables/data:
nChains<-6
cl<-parallel::makeCluster(nChains)
parallel::clusterExport(cl,list("y","X","p","n","nGrps","rJAGSModel3"))

# Do the distributed MCMC sampling:
ptm<-proc.time()
samp<-parallel::clusterApply(cl,1:nChains,parWrapper)
proc.time()-ptm
parallel::stopCluster(cl)

# Save result:
save.image("working_20190121.RData")