#===================================================
#
# Functions used for the eradications tools DSS
#
#===================================================
sim.trap<- function(traps, animals, timehorizon=1, detector.type="traps") {
  # spatial trapping/detection simulator.  Assume half-normal detection
  HN<- function(x, g0, sig) {
    g0 * exp(-x/2/sig^2) 
  }
  
  pars<- animals$marks
  nt<- traps$n
  na<- animals$n
  if(detector.type=="traps") g0<- pars$g0 else g0<- pars$g01
  minprob<- 1e-5 # trap for trivial p
  if(na > 1) {
    dx<- outer(traps$x, animals$x, FUN="-")
    dy<- outer(traps$y, animals$y, FUN="-")
    d2<- dx^2 + dy^2
    p<- t(apply(d2,1,HN, g0=g0,sig=pars$sig))
    p<- -log(1-p)
    U<-  matrix(runif(nt*na),nrow=nt,ncol=na)
    ok<- (p > minprob) & (U > 0) 
    ttime<- -log(U)/p # Random time to next event
    ttime[!ok]<- 1e6 # practically infinite
    trapanimal<- which(ttime <= timehorizon, arr.ind=T)
  }
  else
  { # case of 1 animal left
    dx<- traps$x - animals$x
    dy<- traps$y - animals$y
    d2<- dx^2 + dy^2
    p<- HN(d2, g0=g0, sig=pars$sig)
    p<- -log(1-p)
    U<-  runif(nt*na)
    ok<- (p > minprob) & (U > 0) 
    ttime<- -log(U)/p 
    ttime[!ok]<- 1e6 # practically infinite
    trapID<- which(ttime <= timehorizon, arr.ind=T)
    trapanimal<- cbind(trapID,rep(1,length(trapID)))
  }
  
  trapanimal<- data.frame(trapanimal)
  names(trapanimal)<- c("trap","ID")
  trapanimal
}
#------------------------------
sim.remove<- function(nights, months, traps, animals, move.traps=T, trap.r=500) {
  library(spatstat)
  removed<- matrix(NA, months,2)
  for(i in 1: months) {
    cat(paste("Simulating removal for month ",i,sep=""),"\n")
    count.removed<- 0 # reset each month
    if(move.traps) traps<- rSSI(trap.r, traps, win=region)
    for(j in 1:nights) {
      if(animals$n==0) next
      cap<- sim.trap(traps=traps, animals=animals)
      id.to.remove<- unique(cap$ID)
      count.removed<- count.removed + length(id.to.remove)
      if(nrow(cap) > 0) animals<- animals[-id.to.remove,]
    }
    removed[i,1]<- count.removed
    removed[i,2]<- traps$n*nights
  }
  list(pop=animals,removed=removed)
}
#------------------------------
sim.count<- function(nights, traps, animals) {
  counts<- matrix(0, traps$n, nights) 
  for(j in 1:nights) {
    if(animals$n==0) next
    cap<- sim.trap(traps=traps, animals=animals, detector.type = "cameras")    
    traps.seen<- table(cap$trap)
    counts[as.numeric(names(traps.seen)),j]<- traps.seen    
  }
  counts
}
#------------------------------
sim.remove.monitor<- function(nights, months, animals, traps, monitor=c("none","cams","sign","both"),
                              cams=NULL, sign=NULL, roads=NULL, min.r=1, move.devices=FALSE) {
  library(spatstat)
  removed<- matrix(NA, months,3)
  pop.list<- list()
  monitor<- match.arg(monitor)
  
  if(monitor=="cams") count1<- list()
  if(monitor=="sign") count2<- list()
  if(monitor=="both") {
        count1<- list()
        count2<- list()
        }

  for(i in 1:months) {
    cat(paste("Simulating removal for month ",i,sep=""),"\n")
    count.removed<- 0 # reset each month
    marks.removed<- 0 # reset each month
    if(move.devices)
      traps<- rSSI(min.r, traps$n, win=region)
    for(j in 1:nights) {
      if(animals$n==0) next
      cap<- sim.trap(traps=traps, animals=animals)
      id.to.remove<- unique(cap$ID)
      count.removed<- count.removed + length(id.to.remove)
      if(nrow(cap) > 0) {
            marks.caught<- length(which(animals[id.to.remove]$marks$marked>0))
            marks.removed<- marks.removed+marks.caught
            animals<- animals[-id.to.remove,]
            }
    }
    
    if(monitor=="cams") {
      if(move.devices)    
        cams<- rSSI(min.r, cams$n, win=region)
        count1[[i]]<- sim.count(nights,cams,animals)
    }
    if(monitor=="sign") {
      if(move.devices)    
        sign<- runifpointOnLines(sign$n, roads) 
        count2[[i]]<- sim.count(nights,sign,animals)
    }
    if(monitor=="both") {
      if(move.devices) {
        cams<- rSSI(min.r, cams$n, win=region)
        sign<- runifpointOnLines(sign$n, roads) 
      }
      count1[[i]]<- sim.count(nights,cams,animals)
      count2[[i]]<- sim.count(nights,sign,animals)
    }
    removed[i,1]<- count.removed
    removed[i,2]<- traps$n
    removed[i,3]<- marks.removed
    pop.list[[i]]<- animals
  }
  removed<- data.frame(removed)
  names(removed)<- c("n.removed","n.traps","n.marks")
  if(monitor=="none") res<- list(pop=pop.list,removed=removed)
  if(monitor=="cams") res<- list(pop=pop.list,removed=removed,cam.mon=count1)
  if(monitor=="sign") res<- list(pop=pop.list,removed=removed,sign.mon=count2)
  if(monitor=="both") res<- list(pop=pop.list,removed=removed,cam.mon=count1,sign.mon=count2)
  res
}
#--------------------------------------

gen.pars<- function(N, ming0, maxg0, g0mu, g0shape, minsig, maxsig, sigmu, sigshape, camfct=1) {
# Generate random detection parameters (scaled)  
  g0<- rbeta(N, g0mu/g0shape, (1-g0mu)/g0shape)
  g0<- ming0 + (maxg0 - ming0) * g0
  g01<- g0 * camfct
  
  sig<- rbeta(N, sigmu/sigshape, (1-sigmu)/sigshape)
  sig<- minsig + (maxsig - minsig) * sig
  
  data.frame(g0=g0, g01=g01, sig=sig)
  
}
#---------------------------------------
Fit.NmixModel<- function(Z, nq, n.iter=2000, n.burn=1000, n.chains=3, n.thin=1,msgs=FALSE) {
  # Simple N-mix model assuming independence of sampling locations
  modelFilename = 'NMixModel.txt'
  cat('
      model{
      
      for(i in 1:M){
        N[i] ~ dpois(lambda[i]) 
        log(lambda[i])<- b
        logit(p[i]) <- a       
        for(j in 1:J) {
          Z[i,j] ~ dbin(p[i], N[i])
        }
      }
      a ~ dnorm(0, 0.1) 
      b ~ dnorm(0, 0.1) 
      
      TotalN<- sum(N)
      PopEst<- exp(b)*nq
      }
      ', fill=TRUE, file=modelFilename)
  
  
  M = dim(Z)[1]
  J = dim(Z)[2]
  
  data <- list(M = M, J=J, Z=Z, nq=nq)
  
  # arguments for bugs()
  
  params = c("a","b","TotalN","PopEst")
  inits = function() {
    list(a=rnorm(1), b=rnorm(1), N=apply(Z,1,max)+10)
  }
  
  
  # call to bugs()
  library(jagsUI)
  fit = jags(data, inits, params, modelFilename,
             n.chains=n.chains, n.iter=n.iter, n.burnin=n.burn, n.thin=n.thin,verbose=msgs)
  
  fit$summary
}

#--------------------------------------

Fit.RNModel<- function(Z, nq, n.iter=2000, n.burn=1000, n.chains=3, n.thin=1, msgs=FALSE) {
  # Royle Nichols Occupancy model to estimate N from presence/absence data  
  modelFilename = 'RNmodel1.txt'
  cat('
      model{
      
      for(i in 1:M){
        N[i] ~ dpois(lambda[i]) 
        log(lambda[i])<- b
        mu[i] <- 1-pow(1-r[i],N[i])
        logit(r[i])<- a
        Y[i] ~ dbin(mu[i], J)
      }
      a ~ dnorm(0, 0.1) 
      b ~ dnorm(0, 0.1) 
      
      TotalN<- sum(N)
      PopEst<- exp(b)*nq
      }
      ', fill=TRUE, file=modelFilename)
  
  # Create presence/absence data: Sign and cameras
  Y<- Z
  Y[Y>0]<- 1 
  Y<- apply(Y, 1, sum)
  M = dim(Z)[1]
  J = dim(Z)[2]
  
  data <- list(M = M, J=J, Y=Y, nq=nq)
  
  
  # arguments for bugs()
  
  params = c("a","b","TotalN","PopEst")
  inits = function() {
    list(a=rnorm(1), b=rnorm(1), N=rep(1,M))
  }
  
  
  library(jagsUI)
  fit = jags(data, inits, params, modelFilename,
             n.chains=n.chains, n.iter=n.iter, n.burnin=n.burn, n.thin=n.thin,verbose=msgs)
  
  fit$summary
}
#==================================================
Fit.RemMon<- function(caps, msgs=FALSE) {
  # removal monitoring model
  
  modelFilename = 'removal-monitor.txt'
  cat('
  data {
    # Data statement to calculate cumulative removals just prior to current
    for(j in 1:K) {
      cumy[j]<- sum(y[1:j]) - y[j]
    }
  }
  model {
  for (i in 1:K) {
    y[i] ~ dbin(p1[i], N[i]) # removal trapping
    N[i] <- N0 - cumy[i]
    cloglog(p1[i]) <- a1 + log(ntraps[i])
  }

      N0 ~ dpois(lambda)
      log(lambda)<- u 
      u ~ dnorm(0, 0.1)
      a1 ~ dnorm(0, 0.1)
    
      Nresid<- N0 - sum(y[1:K])
}
      
      ', fill=TRUE, file=modelFilename)
  
  # import data
  rem<- caps$removed
  
  # simulated monitoring data: Sign and cameras
  K<- nrow(rem)
  y<- as.vector(rem[,1])
  ntraps<- as.vector(rem[,2])
  
  data = list(K=K,y=y,ntraps=ntraps)
  
  # arguments for bugs()
  
  params = c('N0','Nresid','a1','N')
  inits = function() {
    list(a1=-3, u=5)
  }
  
  
  # call to bugs()
  library(jagsUI)
  fit = jags(data, inits, params, modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1, verbose=msgs)
  
  fit  
  
}

#==================================================
Fit.OneMon<- function(caps, type=c("cams","sign"), msgs=FALSE) {
  # Combined removal & camera monitoring model
  
  modelFilename = 'removal-monitor.txt'
  cat('
  data {
    # Data statement to calculate cumulative removals just prior to current
    for(j in 1:K) {
      cumy[j]<- sum(y[1:j]) - y[j]
    }
  }
  model {
  for (i in 1:K) {
    y[i] ~ dbin(p1[i], N[i]) # removal trapping
    detected[i] ~ dbin(p2[i], ndevices[i]) # monitoring data
    N[i] <- N0 - cumy[i]
    p2[i] <- 1 - pow(1 - r2[i], N[i])
  
    cloglog(p1[i]) <- a1 + log(ntraps[i])
    logit(r2[i]) <- a2
  }

      N0 ~ dpois(lambda)
      log(lambda)<- u 
      u ~ dnorm(0, 0.1)
      a1 ~ dnorm(0, 0.1)
      a2 ~ dnorm(0, 0.1) 
      
      Nresid<- N0 - sum(y[1:K])
}
      
      ', fill=TRUE, file=modelFilename)
  
  # import data
  rem<- caps$removed
  
  # simulated monitoring data: Sign and one monitor
  K<- nrow(rem)
  y<- as.vector(rem[,1])
  ntraps<- as.vector(rem[,2])
  
  type<- match.arg(type)
  if(type=="cams") tmp<- caps$cam.mon
  if(type=="sign") tmp<- caps$sign.mon
  
  nd<- unlist(lapply(tmp,function(x) dim(x)[1]))
  tmp<- lapply(tmp,function(x) apply(x,1,sum))
  detected<- unlist(lapply(tmp,function(x) length(x[x>0])))
  
  data = list(K=K,y=y,ntraps=ntraps,ndevices=nd,detected=detected)
  
  # arguments for bugs()
  
  params = c('N0','Nresid','a1','a2','N')
  inits = function() {
    list(a1=-3, a2=-2, u=5)
  }
  
  
  # call to bugs()
  library(jagsUI)
  fit = jags(data, inits, params, modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1, verbose=msgs)
  
  fit  
  
}

#=================================================
Fit.BothMon<- function(caps, msgs=FALSE) {
# Combined removal, camera + sign monitoring model
  
  modelFilename = 'removal-monitor.txt'
  cat('
  data {
    # Data statement to calculate cumulative removals just prior to current
    for(j in 1:K) {
      cumy[j]<- sum(y[1:j]) - y[j]
    }
  }
  model {
  for (i in 1:K) {
    y[i] ~ dbin(p1[i], N[i]) # removal trapping
    cams[i] ~ dbin(p2[i], ncams[i]) # cameras
    sign[i] ~ dbin(p3[i], nsign[i]) # sign
    N[i] <- N0 - cumy[i]
    p2[i] <- 1 - pow(1 - r2[i], N[i])
    p3[i] <- 1 - pow(1 - r3[i], N[i])
  
    cloglog(p1[i]) <- a1 + log(ntraps[i])
    logit(r2[i]) <- a2
    logit(r3[i]) <- a3
  }

      N0 ~ dpois(lambda)
      log(lambda)<- u 
      u ~ dnorm(0, 0.1)
      a1 ~ dnorm(0, 0.1)
      a2 ~ dnorm(0, 0.1) 
      a3 ~ dnorm(0, 0.1)
      
      Nresid<- N0 - sum(y[1:K])
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
  data = list(K=K,y=n,ntraps=ntraps,ncams=ncams,cams=cams,nsign=nsign,sign=sign)
  
  # arguments for bugs()
  
  params = c('N0','Nresid','a1','a2','a3','N')
  inits = function() {
    list(a1=-3, a2=-2, a3=-2, u=5)
  }
  
  
  # call to bugs()
  library(jagsUI)
  fit = jags(data, inits, params, modelFilename,
             n.chains=1, n.iter=60000, n.burnin=50000, n.thin=1, verbose=msgs)
  
  #fit$summary  
  fit
}

#===============================================================================

#-----------------------------------------------
# Subset capture data to a set no. of months
caps.subset<- function(caps, n.months) {
 
  caps$removed<- caps$removed[1:n.months,]
  caps$cam.mon<- caps$cam.mon[1:n.months]
  caps$sign.mon<- caps$sign.mon[1:n.months]
  caps
}
#------------------------------------------------
# Summarise simulations
caps.summary<- function(caps) {
  rem<- caps$removed
  sess<- 1:nrow(rem)
  nc<- ns<- NULL
  if(!is.null(caps$cam.mon)) {
    tmp1<- lapply(caps$cam.mon,function(x) apply(x,1,sum))
    nc<- unlist(lapply(tmp1, function(x) length(x[x>0])))
  }
  if(!is.null(caps$sign.mon)) {
    tmp2<- lapply(caps$sign.mon,function(x) apply(x,1,sum))
    ns<- unlist(lapply(tmp2, function(x) length(x[x>0])))           
  } 
  if(is.null(ns) & !is.null(nc))
    data.frame(Sess=sess,Removed=rem[,1],Effort=rem[,2],Marked=rem[,3],Cams=nc)
  else if(is.null(nc) & !is.null(ns))
    data.frame(Sess=sess,Removed=rem[,1],Effort=rem[,2],Marked=rem[,3],Sign=ns)
  else if(!is.null(ns) & !is.null(nc))
    data.frame(Sess=sess,Removed=rem[,1],Effort=rem[,2],Marked=rem[,3],Cams=nc,Sign=ns)
  else
    data.frame(Sess=sess,Removed=rem[,1],Effort=rem[,2],Marked=rem[,3])
}

###################################################
# 
# erad Functions
#
#--------------------------------------------


cell.detect<- function(ndevice, ndays, params) {
  # calculate the probability of detection for a single cell for a given trapping effort
  # alpha is logit daily detection probability per device  
  alpha<- rnorm(ndevice, params$dpar$mean, params$dpar$sd)
  p<- plogis(alpha)
  device.p<- 1-(1-p)^ndays 
  cell.dp<- 1-prod(1-device.p)
  cell.dp
}

#---------------------------------------------

calc.erad<- function(map, nsims, prior.params, params, incRR=FALSE) {
  # calculate posterior probability of eradication, given monitoring detects zero individuals
  # Uses methods given in Anderson, Ramsey et al (2013) Epidemiology and Infection 141:1509-1521 
  ncells<- params$ncells
  ntraps<- params$ntraps
  ndays<- params$ndays
  habcells<- Which(map,cells=T) # ids of cells with habitat
  N<- length(habcells) # total available habitat cells
  if(ncells > N) stop("requested sample size exceeds available habitat cells")
  Pu<- 1/N # Minimum proportion of cells expected to contain individuals given presence
  samp.cells<- sampleRandom(map, size=ncells, cells=T)[,1]
  if(incRR) {
    # map contains relative habitat suitability values
    habk<- map[habcells]
    rri<- habk/min(habk)
    ari<- N * rri/sum(rri)
    ids<- which(!is.na(match(habcells, samp.cells)))
    ar.sum<- sum(ari[ids])
    EPIavg<- Pu * ar.sum/ncells
    expon<- EPIavg * N
  }
  else expon<- Pu * N
  Prior.p<- rbeta(nsims, prior.params$a, prior.params$b)
  cell.sess.dp<- matrix(NA, nsims, ncells)
  for(i in 1:nsims) {
      cell.sess.dp[i, ]<- replicate(ncells, cell.detect(ntraps,ndays,params))
  }
  cell.dp.ave<- apply(cell.sess.dp, 1, mean)
  sims.dp.ave<- apply(cell.sess.dp, 2, mean)
  map[habcells]<- 0 # non sampled cells have 0 sensitivity
  map[samp.cells]<- sims.dp.ave
  SSe<- 1-(1-cell.dp.ave*ncells/N)^expon
  Post.erad<- Prior.p/(1-SSe*(1-Prior.p))
  list(map=map,P.erad=Post.erad,SSe=SSe)
}

#------------------------------------------
mymode<- function(x) {
  z<- density(x,adjust=3,from=0,to=1)
  z$x[which.max(z$y)]
}
#-------------------------------------------
calc.beta.parms<- function(x=NULL, mux=0.5, varx=1){
  #calculate parameters of a beta distribution from a posterior using MOM.
  if(!is.null(x)) {mu<- mean(x);varx<- var(x)}
  a<- mu*((mu*(1-mu))/varx - 1)
  b<- (1-mu)*((mu*(1-mu))/varx - 1)
  xq<- seq(0,1,0.001)
  maxind<- which.max(dbeta(xq, a, b))
  list(a=a,b=b,mode=xq[maxind])
}

#-----------------------------------------------
raster_data <- function(x, maxpixels = 50000)  {
  # extract cell coords and data as tibble
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- as.data.frame(raster::getValues(x)) 
  names(dat) <- c('value')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

#--------------------------------------------------
draw.scale<- function(centre, lth, lab, offset=0.15, ...)
{
  # Draw scale of chosen length at centre co-ords
  lendx <- centre[1.] - lth/2.
  rendx <- centre[1.] + lth/2.
  segments(lendx, centre[2.], rendx, centre[2.])
  segments(lendx, centre[2.], lendx, centre[2.] - lth/10.)
  segments(rendx, centre[2.], rendx, centre[2.] - lth/10.)
  text(centre[1.], centre[2.] - lth * offset, lab, ...)
}
