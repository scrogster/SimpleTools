#===================================================

Fit.NmixModel<- function(Z, nq) {
  
  modelFilename = 'NMixModel.txt'
  cat('
      model{
      
      for(i in 1:M){
      # ecological process model  
      N[i] ~ dpois(lambda[i]) 
      log(lambda[i])<- b
      }
      for(i in 1:n) {
      C[i] ~ dbin(p[i], N[site[i]])
      logit(p[i])<- a       
      }
      
      a ~ dnorm(0,0.01)I(-20,20) 
      b ~ dnorm(0,0.01)I(-20,20) 
      
      TotalN<- sum(N[])
      PopEst<- exp(b)*nq
      }
      ', fill=TRUE, file=modelFilename)
  
  Z<- Z[,1:10]
  M = dim(Z)[1]
  J = dim(Z)[2]
  n = M*J
  site = rep(1:M,J)
  C<- c(Z)
  
  data <- list(M = M, n=n, site=site, C=C, nq=nq)
  
  # arguments for bugs()
  
  params = list('a','b','TotalN','PopEst')
  inits = function() {
    list(a=0, b=0, N=apply(Z,1,max)+1)
    #list(a=rep(0,5), b=rep(0,5),u=10,prior.p=0.5)
  }
  
  
  # call to bugs()
  library(R2OpenBUGS)
  fit = bugs(data, inits, params, model.file=modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1,summary.only=T)
  
  est<- fit$stats[,c(1:5)]
  est
}

#--------------------------------------

Fit.RNModel<- function(Z, nq) {
  
  modelFilename = 'RNmodel1.txt'
  cat('
      model{
      
      for(i in 1:M){
      # ecological process model
      
      N[i] ~ dpois(lambda[i]) 
      log(lambda[i])<- b
      mu[i] <- 1-pow(1-r[i],N[i])
      logit(r[i])<- a
      Y[i] ~ dbin(mu[i],J)
      }
      a ~ dnorm(0,0.01)I(-20,20) 
      b ~ dnorm(0,0.01)I(-20,20) 
      
      TotalN<- sum(N[])
      PopEst<- exp(b)*nq
      }
      ', fill=TRUE, file=modelFilename)
  
  # simulated monitoring data: Sign and cameras
  
  
  Y<- Z
  Y[Y>0]<- 1 
  Y<- apply(Y, 1, sum)
  M = dim(Z)[1]
  J = dim(Z)[2]
  
  data <- list(M = M, J=J, Y=Y, nq=nq)
  
  
  # arguments for bugs()
  
  params = list('a','b','TotalN','PopEst')
  inits = function() {
    list(a=0, b=0, N=rep(1,M))
  }
  
  
  # call to bugs()
  library(R2OpenBUGS)
  fit2 = bugs(data, inits, params, model.file=modelFilename,
              n.chains=3, n.iter=60000, n.burnin=50000, n.thin=1,summary.only=T)
  
  
  est<- fit2$stats[,c(1:5)]
  est
}

#==================================================

Fit.CompMon<- function(caps, PopMax) {
  # model specification in WinBUGS
  # model specification in WinBUGS
  modelFilename = 'composite-monitor.txt'
  cat('
      model {
      # declining catchability model
      for(i in 1:K) {
      
      # removal part
      n[i] ~ dbin(p1[i], N[i]) # removal trapping 
      cams[i] ~ dbin(p2[i],ncams[i]) # cameras
      sign[i] ~ dbin(p3[i],nsign[i]) # sign
      N[i+1] <- N[i] - n[i]
      p2[i] <- 1-pow(1-r2[i],N[i])
      p3[i] <- 1-pow(1-r3[i],N[i])
      
      logit(p1[i])<- a1 + b1*ntraps[i]
      logit(r2[i])<- a2 
      logit(r3[i])<- a3  
      
      }
      
      # prior on initial population size N0
      
      N[1] <- round(exp(u))
      nsum<- log(sum(n[]))
      u ~ dunif(nsum,PopMax)
      a1 ~ dunif(-20,10)
      a2 ~ dunif(-20,10) 
      a3 ~ dunif(-20,10)
      b1 ~ dunif(-10,10)
      
      Nresid<- N[1] - sum(n[])
      }
      
      ', fill=TRUE, file=modelFilename)
  
  # import data
  rem<- caps$removed
  
  # simulated monitoring data: Sign and cameras
  K<- nrow(rem)
  n<- as.vector(rem[,1])
  ntraps<- as.vector(rem[,2])
  
  tmp<- caps$cam.mon
  ncams<- unlist(lapply(tmp,function(x) dim(x)[1]))
  tmp<- lapply(tmp,function(x) apply(x,1,sum))
  cams<- unlist(lapply(tmp,function(x) length(x[x>0])))
  
  tmp<- caps$sign.mon
  nsign<- unlist(lapply(tmp,function(x) dim(x)[1]))
  tmp<- lapply(tmp,function(x) apply(x,1,sum))
  sign<- unlist(lapply(tmp,function(x) length(x[x>0])))
  
  
  #data = list(K=nrow(rem),n=n,eff=eff,ncams=ncams,cams=cams,nsign=nsign,sign=sign)
  data = list(K=K,n=n,ntraps=ntraps,ncams=ncams,cams=cams,nsign=nsign,sign=sign,PopMax=PopMax)
  
  # arguments for bugs()
  
  params = c('N','Nresid','a1','a2','a3','b1')
  inits = function() {
    list(a1=0, a2=0, a3=0, b1=0, u=5)
    #list(a=rep(0,5), b=rep(0,5),u=10,prior.p=0.5)
  }
  
  
  # call to bugs()
  library(R2OpenBUGS)
  fit = bugs(data, inits, params, model.file=modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1,summary.only=T,clearWD=T)
  
  
  est<- fit$stats[,c(1:5)]
  est
  
}
#==================================================

Fit.CamMon<- function(caps, PopMax) {
  # model specification in WinBUGS
  # model specification in WinBUGS
  modelFilename = 'composite-monitor.txt'
  cat('
      model {
      # declining catchability model
      for(i in 1:K) {
      
      # removal part
      n[i] ~ dbin(p1[i], N[i]) # removal trapping 
      cams[i] ~ dbin(p2[i],ncams[i]) # cameras
      N[i+1] <- N[i] - n[i]
      p2[i] <- 1-pow(1-r2[i],N[i])
      
      logit(p1[i])<- a1 + b1*ntraps[i]
      logit(r2[i])<- a2 
      
      }
      
      # prior on initial population size N0
      
      N[1] <- round(exp(u))
      nsum<- log(sum(n[]))
      u ~ dunif(nsum,PopMax)
      a1 ~ dunif(-20,10)
      a2 ~ dunif(-20,10) 
      b1 ~ dunif(-10,10)
      
      Nresid<- N[1] - sum(n[])
      }
      
      ', fill=TRUE, file=modelFilename)
  
  # import data
  rem<- caps$removed
  
  # simulated monitoring data: Sign and cameras
  K<- nrow(rem)
  n<- as.vector(rem[,1])
  ntraps<- as.vector(rem[,2])
  
  tmp<- caps$cam.mon
  ncams<- unlist(lapply(tmp,function(x) dim(x)[1]))
  tmp<- lapply(tmp,function(x) apply(x,1,sum))
  cams<- unlist(lapply(tmp,function(x) length(x[x>0])))
  
  #data = list(K=nrow(rem),n=n,eff=eff,ncams=ncams,cams=cams,nsign=nsign,sign=sign)
  data = list(K=K,n=n,ntraps=ntraps,ncams=ncams,cams=cams,PopMax=PopMax)
  
  # arguments for bugs()
  
  params = c('N','Nresid','a1','a2','b1')
  inits = function() {
    list(a1=0, a2=0, b1=0, u=5)
    #list(a=rep(0,5), b=rep(0,5),u=10,prior.p=0.5)
  }
  
  
  # call to bugs()
  library(R2OpenBUGS)
  fit = bugs(data, inits, params, model.file=modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1,summary.only=T,clearWD=T)
  
  
  est<- fit$stats[,c(1:5)]
  est
  
}
#==================================================
Fit.RemMon<- function(caps, PopMax) {
  # model specification in WinBUGS
  # model specification in WinBUGS
  modelFilename = 'removal-monitor.txt'
  cat('
      model {
      # declining catchability model
      for(i in 1:K) {
      
      # removal part
      n[i] ~ dbin(p1[i], N[i]) # removal trapping 
      N[i+1] <- N[i] - n[i]
      
      logit(p1[i])<- a1 + b1*ntraps[i]
      
      
      }
      
      # prior on initial population size N0
      
      N[1] <- round(exp(u))
      nsum<- log(sum(n[]))
      u ~ dunif(nsum,PopMax)
      a1 ~ dunif(-20,10)
      b1 ~ dunif(-10,10)
      
      Nresid<- N[1] - sum(n[])
      }
      
      ', fill=TRUE, file=modelFilename)
  
  # import data
  rem<- caps$removed
  
  # simulated monitoring data: Sign and cameras
  K<- nrow(rem)
  n<- as.vector(rem[,1])
  ntraps<- as.vector(rem[,2])  
  
  
  #data = list(K=nrow(rem),n=n,eff=eff,ncams=ncams,cams=cams,nsign=nsign,sign=sign)
  data = list(K=K,n=n,ntraps=ntraps,PopMax=PopMax)
  
  # arguments for bugs()
  
  params = c('N','Nresid','a1','b1')
  inits = function() {
    list(a1=0, b1=0, u=5)
    #list(a=rep(0,5), b=rep(0,5),u=10,prior.p=0.5)
  }
  
  
  # call to bugs()
  library(R2OpenBUGS)
  fit = bugs(data, inits, params, model.file=modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1,summary.only=T,clearWD=T)
  
  est<- fit$stats[,c(1:5)]
  est
  
}
#------------------------------------------------------------
caps.summary<- function(caps) {
  rem<- caps$removed
  sess<- 1:nrow(rem)
  nc<- NA
  ns<- NA
  if(!is.null(caps$cam.mon)) {
    tmp1<- lapply(caps$cam.mon,function(x) apply(x,1,sum))
    nc<- unlist(lapply(tmp1, function(x) length(x[x>0])))
  }
  if(!is.null(caps$sign.mon)) {
    tmp2<- lapply(caps$sign.mon,function(x) apply(x,1,sum))
    ns<- unlist(lapply(tmp2, function(x) length(x[x>0])))        
  }            
  data.frame(Sess=sess,Removed=rem[,1],Effort=rem[,2],Marked=rem[,3],Cams=nc,Sign=ns)
}
#-----------------------------------------------------------

Update.prior <- function(par1,par2=NULL,prior.par,ncams,nsign=NULL,roads=NULL,region=NULL,rad=NULL,n=10000){
  
  # Update posterior given prior and parameters
  
  prior<- rbeta(n,prior.par[1],prior.par[2])  # uniform
  
  # determine monitoring coverage
  ctraps<- runifpoint(ncams,win=region)
  straps<- runifpointOnLines(nsign,roads)
  cdist<- distmap(ctraps)
  sdist<- distmap(straps)
  cdist$v[cdist$v>rad]<- NA
  sdist$v[sdist$v>rad]<- NA
  ccov<- area.owin(cdist)/area.owin(region)
  scov<- area.owin(sdist)/area.owin(region)
  
  if(is.null(par2)) {
    lp<- rnorm(n,par1[1],par1[2])
    p<- plogis(lp)
    sens<- 1-(1-p)^ncams
    sens<- sens*ccov
    posterior<- (1-sens)*prior/((1-sens)*prior + (1-prior))
  }
  else
    if(!is.null(par2)) {
      lp1<- rnorm(n,par1[1],par1[2])
      lp2<- rnorm(n,par2[1],par2[2])
      p1<- plogis(lp1)
      p2<- plogis(lp2)
      sens1<- 1-(1-p1)^ncams
      sens2<- 1-(1-p2)^nsign
      sens1<- sens1*ccov
      sens2<- sens2*scov
      sens<- 1-(1-sens1)*(1-sens2)
      posterior<- (1-sens)*prior/((1-sens)*prior + (1-prior))    
    }
  
  #   Calculate parameters of posterior distribution [which is Beta-Binomial]
  varhat <- var(posterior)
  sd.hat<- sd(posterior)
  mu  <- mean(posterior)
  alpha.post <- mu*((mu*(1-mu))/varhat - 1)
  beta.post <- (1-mu)*((mu*(1-mu))/varhat - 1)
  quan<- quantile(posterior,p=c(0.025,0.975))
  
  c(mu,sd.hat,quan,alpha.post,beta.post)
}

#-------------------------------------------
calc.trap.sens <- function(par1,par2,ntraps,cover,n.resid,n=10000){
  
  n.units<- length(ntraps)
  sens.mat<- matrix(NA,nrow=n,ncol=n.units)
  rem.mat<- matrix(NA,nrow=n,ncol=n.units)
  for(i in 1:n.units) {
    
    a<- rnorm(n,par1[1],par1[2])
    b<- rnorm(n,par2[1],par2[2])
    lp<- a + b*ntraps[i]
    p<- plogis(lp)
    sens<- p*cover[i]
    sens.mat[,i]<- sens
    rem.mat[,i]<- rbinom(n,n.resid,sens)
    
  }
  
  #   Calculate parameters of posterior distribution [which is Beta-Binomial]
  mu  <- apply(sens.mat,2,mean)
  quan<- apply(sens.mat, 2, quantile,c(0.025,0.975))
  mu.rem<- apply(rem.mat,2,mean)
  quan.rem<- apply(rem.mat, 2, quantile,c(0.025,0.975))
  data.frame(mu=mu,lcl=quan[1,],ucl=quan[2,],n.rem=mu.rem,r.lcl=quan.rem[1,],r.ucl=quan.rem[2,])
}
#-------------------------------------------
calc.mon.sens <- function(par1,ntraps,cover,n=10000){
  
  n.units<- length(ntraps)
  sens.mat<- matrix(NA,nrow=n,ncol=n.units)
  for(i in 1:n.units) {
    
    a<- rnorm(n,par1[1],par1[2])
    p<- plogis(a)
    sens<- 1-(1-p)^ntraps[i]
    sens<- sens*cover[i]
    sens.mat[,i]<- sens
    
  }
  
  #   Calculate parameters of posterior distribution [which is Beta-Binomial]
  mu  <- apply(sens.mat,2,mean)
  quan<- apply(sens.mat, 2, quantile,c(0.025,0.975))
  
  data.frame(mu=mu,lcl=quan[1,],ucl=quan[2,])
}

#-----------------------------------------------
calc.coverage.cams<- function(ncams,rad=10,region) {
  # determine monitoring coverage
  ctraps<- runifpoint(ncams,win=region)
  cdist<- distmap(ctraps)
  cdist$v[cdist$v>rad]<- NA
  ccov<- area.owin(cdist)/area.owin(region) 
  ccov[ccov>1]<- 1
  ccov[ccov<0]<- 0
  ccov
}
#-----------------------------------------------
calc.coverage.sign<- function(nsign,rad=10,roads,region) {
  # determine monitoring coverage
  straps<- runifpointOnLines(nsign,roads)
  sdist<- distmap(straps)
  sdist$v[sdist$v>rad]<- NA
  scov<- area.owin(sdist)/area.owin(region)
  scov[scov>1]<- 1
  scov[scov<0]<- 0
  scov
}
#-----------------------------------------------
