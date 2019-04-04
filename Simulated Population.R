#===================================================================================
#
# Script to simulate the eradication of a simulated population from a defined area
#
#====================================================================================

library(spatstat)
library(tidyverse)
library(sf)
library(raster)
library(maptools)
#source utility functions
source("Framework_functions.R")

shape<- read_sf("Stuff/.","San_Nicolas_projected")
shape<- filter(shape, SP_ID==0) # ignore outlying islands

bb<- st_bbox(shape)
res<- 1000 # cell size in meters
rast<- raster(xmn=bb$xmin,xmx=bb$xmax,ymn=bb$ymin,ymx=bb$ymax,resolution=res)
rast<- rasterize(shape, rast, field=1)
rast_df<- raster_data(rast)

shape %>% ggplot() +
  geom_tile(aes(x, y, fill=value), color="grey50", data=rast_df, show.legend = F) +
  scale_fill_gradientn(colors="lightblue",na.value="white") +
  geom_sf(fill=NA) +
  labs(x="Easting",y="Northing",title="San Nicolas Island")

area<- units::set_units(st_area(shape), km^2) 

region<- as.owin.SpatialPolygons(as(shape,"Spatial")) # convert to owin object
#----------------------------------------------------------------
# generate cat population (HR centre locations)
pop<- runifpoint(150, win=region) # Cat population
plot(pop, pch=16)

#-----------------------------------------------------------------
# Set some parameters for trap capture

catHR<- list(min=5e5,max=3e6) # min and max cat home range size in m2
sigmin<- sqrt(catHR$min /pi)/2.45   # min and max sigma
sigmax<- sqrt(catHR$max /pi)/2.45

sigmu<- 0.4             
sigshape=0.1

g0min<- 0.01
g0max<- 0.2

g0mu<- 0.2
g0shape=0.1

cam.fct<- 2 # Increase in g0 for cameras cw traps..

pars<- gen.pars(pop$n, g0min, g0max, g0mu, g0shape, sigmin, sigmax, sigmu, sigshape, cam.fct)
marks(pop)<- pars
          
#----------------------------
# put some extra marks on a random sample of 20 individuals
marked<- rbinom(npoints(pop), 1, 0.2)
pars$marked<- marked
marks(pop)<- pars
#----------------------------------------------
# Plot distributions of parameters and other stuff


layout(matrix(1:4, nrow=2, byrow=TRUE))
  xx<- seq(0,1,0.001)
  dg0<- dbeta(xx,g0mu/g0shape,(1-g0mu)/g0shape)
  xx<- g0min + (g0max - g0min) * xx
  plot(xx,dg0,type='l',xlim=c(0,g0max),xlab="g0",ylab="Probability density")
  
  xx<- seq(0,1,0.001)
  dsig<- dbeta(xx,sigmu/sigshape,(1-sigmu)/sigshape)
  xx<- sigmin + (sigmax - sigmin) * xx
  xx<- pi*(xx*2.45)^2/1e6 # 95% HR in km2
  plot(xx,dsig,type='l',xlim=c(0 ,3),xlab=expression(paste("Home range size (",km^2,")")),ylab="Probability   density")
  
  plot(region, main="")
  symbols(pop$x,pop$y,circles=pop$marks$sig*2.45,fg="gray50",add=T,inches=F)
  points(pop$x,pop$y,pch=19,cex=0.5)

  plot(pop, use.marks=F,main="")
  draw.scale(c(1854000,471000),1000,"1 km",cex=0.8,offset=0.3)

#=======================================================================================
#
# # Pre - implementation monitoring
#
#=======================================================================================

# Camera Trap Locations

ncams<- 40  # Set number of cameras
  
traps<- as_tibble(sampleRandom(rast, size=ncams, xy=TRUE)) %>% dplyr::select(x, y)  # one camera per 1km2

layout(matrix(1))
shape %>% ggplot() +
  geom_tile(aes(x, y, fill=value), color="grey50", data=rast_df, show.legend = F) +
  scale_fill_gradientn(colors="lightblue",na.value="white") +
  geom_point(aes(x, y), data=traps) +
  geom_point(aes(x, y), data=as.data.frame(pop),shape=1, col="red") +
  geom_sf(fill=NA) +
  labs(x="Easting",y="Northing",title="San Nicolas Island")

straps<- as.ppp(traps, region)

nights<- 30
counts<- sim.count(nights, straps, pop)

#===================================
# Save simulated data
saveRDS(counts, "Stuff/pcounts.rds")

#===================================
#  do some estimates

nq<- sum(rast_df$value, na.rm=TRUE)  # total number of 1km2 quadrats

est1<- Fit.NmixModel(counts, nq, n.iter=20000, n.burn=10000, n.chains=3, msgs=T)

est2<- Fit.RNModel(counts, nq, n.iter=20000, n.burn=10000, n.chains=3, msgs=T)

# Estimates 
cat(paste("Nmix Population size = ",round(est1["PopEst","mean"])," with 95% credible interval of ",
          round(est1["PopEst","2.5%"]),"-",round(est1["PopEst","97.5%"]),sep=""),"\n")

cat(paste("RN Population size = ",round(est2["PopEst","mean"])," with 95% credible interval of ",
          round(est2["PopEst","2.5%"]),"-",round(est2["PopEst","97.5%"]),sep=""),"\n")


cat(paste("Mean detection probability Nmix model = ", round(plogis(est1["a","mean"]),3)),"\n")
cat(paste("Mean detection probability RN model = ", round(plogis(est2["a","mean"]),3)),"\n")

# Compare to theoretical mean detection probability
# Assume effective radius is 2.45 * sigma (conservative)

HN<- function(x, g0, sig, w) {
  px<- 2 * x / w^2
  e<- g0 * exp(-x^2/(2*sig^2))
  px * e
}

sig<- mean(marks(pop)$sig)
g01<-mean(marks(pop)$g01) # mean camera g0
w<- 2.45*sig
integrate(HN, lower=0, upper=w, g0=g01, sig=sig, w=w)

#=======================================================================================
#
# Implementation - removal trapping + extra camera + sign monitoring
#
#=======================================================================================

# Example monitoring using leg-hold traps (removal) and two additional monitoring methods
# Cameras and sand plots on roads

# Simulate some monitoring locations

# removal traps (150 traps with min 250m spacing)
rtraps<- rSSI(r=250, n=150, win=region)

# 30 random camera traps at min 1000m spacing
cam.traps<- rSSI(r=1000, n=30, win=region)

# Random roads for sign surveys
roads<- rpoisline(2e-4, win=region) 

# random plots on roads
sign.traps<- runifpointOnLines(40, roads)

# Plot these
layout(matrix(1))
plot(roads,main="",col="red")
plot(sign.traps,add=T)
plot(rtraps,add=T,chars=15)
plot(cam.traps, add=T, chars=16, col="blue")
draw.scale(c(1854000,471000),1000,"1 km",cex=0.8,offset=0.3)

# Eradication program parameters----------------

sessions<- 30  # unique sessions
nights<- 10  #consecutive nights per session

caps<- sim.remove.monitor(nights,sessions,animals=pop, traps=rtraps,monitor="both",
                          cams=cam.traps,sign=sign.traps,roads=roads,min.r=250,move.devices = TRUE)

caps.summary(caps)

# Estimation
#rem1<- Fit.RemMon(caps)
#rem1<- Fit.OneMon(caps, type="cams", msgs = T)
#rem1<- Fit.BothMon(caps)
#rem1

saveRDS(caps, file="Stuff/removal_data.rds")

#------------------------------------------
# Subsample sessions if needed
nm<-10  # only 10 sessions
caps.sub<- caps.subset(caps, nm)

saveRDS(caps.sub, file="Stuff/removal_data_sub.rds")

#=======================================================================================
#
# Decision support - to be applied once monitoring detects no individuals
#
#=======================================================================================

par1<- data.frame(mean=-5.0, sd=1) # detection probability device per-day (logit scale)
nsims<- 10000 
n.sessions<-20  # Montioring sessions
ndays<- 20 # days per monitoring session
ndevices<- 2 # devices per cell
ncells<-  nq * 0.8  # Number of cells monitored as proportion of total 


mat<- matrix(NA, nrow=n.sessions, ncol=4)
colnames(mat)<- c("Session","P.erad","LCL","UCL")
params<- list(ntraps=ndevices, ndays=ndays, ncells=ncells, dpar=par1)

prior.params<- list(a=1,b=1) # initial uninformative prior
res= calc.erad(rast, nsims, prior.params, params)
updated.priors<- calc.beta.parms(res$P.erad)
mat[1,]<- c(1,qbeta(0.5,1.5,1.5),qbeta(0.025,1.5,1.5),qbeta(0.975,1.5,1.5))

# Update for {n} monitoring sessions

for(i in 2:n.sessions) {
  cat("Updating session ",i,"\n")
  res= calc.erad(rast, nsims, updated.priors, params)
  updated.priors<- calc.beta.parms(res$P.erad)
  mat[i,]<- c(i,median(res$P.erad),quantile(res$P.erad,c(0.025,0.975)))
}

mat<- as_tibble(mat)
#------------------------------------------------------------------------------------------
# Various plots


mat %>% ggplot(aes(Session)) +
  geom_ribbon(aes(ymin=LCL,ymax=UCL),fill="grey70") +
  geom_line(aes(y=P.erad), size=1.5) +
  labs(x="Session",y="Probability of Eradication") +
  geom_hline(yintercept=0.95, linetype=2, colour="red", size=1.5)+
  theme_bw()

tmp<- data.frame(SSe=res$SSe)
tmp %>% ggplot(aes(SSe)) +
  geom_histogram(binwidth=0.01, fill="lightblue", colour="black") +
  xlim(0, 1)+
  theme_bw()

tmp<- raster_data(res$map)
shape %>% ggplot() +
  geom_raster(aes(x, y, fill=value), data=tmp, interpolate=F) +
  scale_fill_distiller(palette="Spectral",na.value="white", limits=c(0, 1))  +
  geom_sf(fill=NA) +
  labs(x="Easting",y="Northing",title="Map of system sensitivity (SSe)")+
  theme_bw()



