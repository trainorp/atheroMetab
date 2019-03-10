############ Prereqs ############
## Start always run:
options(stringsAsFactors=FALSE,scipen=900)
library(pillar)
library(tidyverse)

load(file="working_20190223.RData")

metabInclude<-cpsTemp$metabID[1:15]
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
  tau~dgamma(2,1)
  SD<-sqrt(1/tau)
}"

tropT0<-oxPLDF %>% select(ptid,tropT0) %>% mutate(logTrop=log(tropT0+.0001)) %>% select(-tropT0)
df2bT0<-df2bT0 %>% left_join(tropT0)
X<-scale(df2bT0[,names(df2bT0) %in% c(metabInclude,"logTrop")])
# X<-scale(df2bT0[,names(df2bT0) %in% metabInclude])
p<-dim(X)[2]
cat("p is ",p,"\n")

# Sample for cross-validation
library(doParallel)
cl<-makeCluster(4)
registerDoParallel(cl)
ptm<-proc.time()
codaSamples<-foreach(i=1:nrow(X),.inorder=FALSE) %dopar% {
  X2<-X[-i,]
  y<-y[-i]
  # tempGroupInd<-df2bT0$group[-1]
  # resampled<-sample(which(tempGroupInd=="Non-Thrombotic MI"),4)
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


# Not LOO model:
n<-dim(X)[1]
set.seed(3333333)
model<-rjags::jags.model(file=textConnection(rJAGSModel2),
                         data=list(y=y,X=X,p=p,n=n,nGrps=nGrps),n.chains=6,n.adapt=1000)
codaSamplesOneModel<-rjags::coda.samples(model,
                                         variable.names=c("logdens","tau","SD","beta0","beta"),n.iter=10000,thin=10)
codaSamplesOneModel<-as.data.frame(do.call("rbind",codaSamplesOneModel))

codaSOMParamQuant<-data.frame()
for(colName in colnames(codaSamplesOneModel)){
  codaSOMParamQuant<-rbind(codaSOMParamQuant,t(as.matrix(quantile(codaSamplesOneModel[,colName],probs=seq(0,1,.1)))))
}
codaSOMParamQuant$param<-colnames(codaSamplesOneModel)

############ T0 Bayesian model prediction from LOO-CV ############
# Combind sets of chains:
codaSamples<-do.call("rbind",codaSamples)

# Calculate group probabilities from LOO-CV posteriors
groupExpList<-groupProbsList<-list()
for(i in 1:nrow(codaSamples)){
  groupExp<-matrix(0,nrow=1,ncol=3)
  colnames(groupExp)<-levels(as.factor(df2bT0$group))
  for(g in 2:3){
    betaVars<-paste0("beta[",g,",",1:p,"]")
    codaSampBeta<-codaSamples[i,match(betaVars,colnames(codaSamples))]
    codaSampBeta0<-codaSamples[i,match(paste0("beta0[",g,"]"),colnames(codaSamples))]
    groupExp[1,g]<-exp(codaSampBeta0 + X[codaSamples$lo[i],] %*% t(codaSampBeta))
  }
  groupExp[,1]<-1
  groupProbs<-groupExp/apply(groupExp,1,sum)
  groupProbs<-data.frame(ptid=df2bT0$ptid[codaSamples$lo[i]],groupProbs)
  
  # Add to lists:
  groupExp<-as.data.frame(groupExp)
  groupExp$iter<-i
  groupProbs$iter<-i
  groupExpList[[i]]<-groupExp
  groupProbsList[[i]]<-groupProbs
  if(i %% 2880==0){
    print(i/nrow(groupProbs)*100)
  }
}
groupExp<-do.call("rbind",groupExpList)
groupProbs<-do.call("rbind",groupProbsList)

# Multinomial loss:
groupProbs<-df2bT0 %>% select(ptid,group) %>% left_join(groupProbs)
groupProbs$ind<-match(gsub("-",".",gsub(" ",".",groupProbs$group)),names(groupProbs))
groupProbs$multLoss<-NA
for(i in 1:nrow(groupProbs)){
  groupProbs$multLoss[i]<-groupProbs[i,groupProbs$ind[i]]
  if(i %% 2880==0){
    print(i/nrow(groupProbs)*100)
  }
}
groupProbs$multLoss<-unlist(groupProbs$multLoss)
groupProbs$multLoss<-(-log(groupProbs$multLoss))
save.image(file="working_20190223c.RData")

# CV-estimated confusion matrix:
groupProbsSum<-groupProbs %>% select(ptid,Thrombotic.MI,Non.Thrombotic.MI,sCAD) %>% group_by(ptid) %>% 
  summarize(Thrombotic.MI=median(Thrombotic.MI),Non.Thrombotic.MI=median(Non.Thrombotic.MI),
            sCAD=median(sCAD))
groupProbsSum$predGroup<-c("Thrombotic MI","Non-Thrombotic MI","sCAD")[apply(groupProbsSum[,
                                                                                           c("Thrombotic.MI","Non.Thrombotic.MI","sCAD")],1,which.max)]
groupProbsSum<-df2bT0 %>% select(ptid,group) %>% left_join(groupProbsSum,by="ptid")
xtabs(~group+predGroup,data=groupProbsSum)