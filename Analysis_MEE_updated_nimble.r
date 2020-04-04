#To run this example code, set the working directory below and place the example dataset in that directory
#Also, be sure to have installed the 'rjags' library from CRAN
rm(list=ls())
gc()

setwd("C:\\Users\\conve\\Documents\\Projects - Old\\Hierarchical Capture-Recapture Projects\\Andy Royle HSCR") 
library(nimble)

load("Example_Data_MEE.RData")
peromyscus.data <- peromyscus.data
group.mem <- group.mem
S.st <- S.st          

################################################################################
#BUGS CODE 

replicated_scr_model <- nimbleCode({

#PRIORS
#abundance model
for(i in 1:M){
  group.mem[i] ~ dcat(gprobs[1:n.groups])
  z[i] ~ dbern(psi)
}
psi ~ dunif(0,1)
for(i in 1:n.sites){
  b.site[i] ~ dnorm(int.lam,tau.site)
}
int.lam <- 0
sigma.site ~ dunif(0,10)
tau.site<-1/(sigma.site*sigma.site)
for(i in 1:2){
  b.season[i] <- b.seas[i]
  b.seas[i] ~ dunif(-10,10)
}
b.season[3] <- -1*(b.season[1]+b.season[2])

b.fire ~ dunif(-10,10)
b.thin ~ dunif(-10,10)

#observation model
for(i in 1:M){
  cent[i,1] ~ dunif(Xl,Xu)
  cent[i,2] ~ dunif(Yl,Yu)
}

for(s in 1:24){
  bgroup.p[s] ~ dnorm(int.p,tau.p)
}

int.p ~ dunif(-10,10)
sigma.p ~ dunif(0,10)
tau.p <- 1/(sigma.p*sigma.p)

bcap.p ~ dunif(-10,10)

bcap.am ~ dunif(-10,10)

for(i in 1:24){
  b.dist[i] ~ dnorm(int.dist,tau.dist)
}
int.dist ~ dunif(-10,10)
sigma.dist ~ dunif(0,10)
tau.dist <- 1/(sigma.dist*sigma.dist)

#LIKELIHOOD
#abundance model
for(j in 1:n.groups){
  log(lam[j]) <- b.site[site[j]] + b.season[season[j]] + b.thin*thin[j] + b.fire*fire[j]
  gprobs[j] <- lam[j]/sum(lam[1:n.groups])
}

#observation model for observed individuals
for(i in 1:nind){
 # assumes traplocs are the SAME for all groups 
 for(k in 1:last.cap[i]){
  cp[i,k,1:(n.traps+1)] <- GetMnlDetProb(cent=cent[i,1:2], bgroup.p=bgroup.p[1:n.groups], group.mem=group.mem[i], am.effect=am.effect[k], bcap.am=bcap.am,
                                         bcap.p=bcap.p, reencounter=reencounter[i,k], b.dist=b.dist[1:n.groups],  
                                         traps.avail = traps.avail[1:n.traps,1:n.groups], z=z[i],
                                         n.groups=n.groups, n.traps=n.traps, trap.locs=trap.locs[1:n.traps,1:2])
  Ycat[i,k] ~ dcat(cp[i,k,1:(n.traps+1)])
  }  
}   

#observation model for unobserved individuals (chandler trick)
for(i in (nind+1):M){
 Y0[i] ~ dbern(pStar[i])
 # 1 - probablity missed on every occasion 
 pStar[i] <-  GetMnlPstar(cent=cent[i,1:2], bgroup.p=bgroup.p[1:n.groups], group.mem=group.mem[i], am.effect=am.effect[1:last.cap[i]], bcap.am=bcap.am,
                          bcap.p=bcap.p, reencounter=reencounter[i,1:last.cap[i]], b.dist=b.dist[1:n.groups],  
                          traps.avail = traps.avail[1:n.traps,1:n.groups], z=z[i],
                          n.groups=n.groups, n.traps=n.traps, trap.locs=trap.locs[1:n.traps,1:2], 
                          last.cap=last.cap[i])
}  

#DERIVED PARAMETERS
#The G.N are the total population sizes by group
N.tot <- sum(z[1:M]) 
G.N[1:n.groups] <- getGroupAbundance(z = z[1:M], group.mem=group.mem[1:M], n.groups=n.groups, M=M)

}) # model


## function to derive occasion-specific multinomial detection probability 
GetMnlDetProb<- nimbleFunction(
  run = function(cent=double(1), bgroup.p=double(1), group.mem=double(0), am.effect=double(0), bcap.am=double(0),
                                         bcap.p=double(0), reencounter=double(0), b.dist=double(1), 
                                         traps.avail = double(2), z=double(0),
                                         n.groups=double(0), n.traps=double(0), trap.locs=double(2)){ 
   returnType(double(1))
   cp <- rep(0, n.traps+1)

   if(z==0) return(c( rep(0, n.traps), 1) ) # cannot be detected

   if(z==1) {
     d <- pow(pow(cent[1]-trap.locs[1:n.traps,1],2) + pow(cent[2]-trap.locs[1:n.traps,2],2),.5)
     lp <- (exp(bgroup.p[group.mem] + am.effect*bcap.am + bcap.p*reencounter + b.dist[group.mem]*d[1:n.traps])*traps.avail[1:n.traps,group.mem])          
     cp[1:n.traps] <- lp[1:n.traps]/(1+sum(lp[1:n.traps]))
     cp[n.traps+1] <- 1-sum(cp[1:n.traps])  # last cell = not captured
     return(cp)
   }#if(z==1)
  }
)

## function to derive total detection probabiltiy for augmented individuals
GetMnlPstar<- nimbleFunction(
 run = function(cent=double(1), bgroup.p=double(1), group.mem=double(0), am.effect=double(1), bcap.am=double(0),
                                         bcap.p=double(0), reencounter=double(1), b.dist=double(1), 
                                         traps.avail = double(2), z=double(0),
                                         n.groups=double(0), n.traps=double(0), trap.locs=double(2),last.cap=double(0) ){
   returnType(double(0))
   if(z==0)  pstar <- 0       #detection prob = 0 if z=0

   if(z==1) {
    cp <- matrix(0, nr=n.traps+1, nc=last.cap)
    # distance doesn't change
    d <- pow(pow(cent[1]-trap.locs[1:n.traps,1],2) + pow(cent[2]-trap.locs[1:n.traps,2],2),.5)
    for(k in 1:last.cap){
     lp <- (exp(bgroup.p[group.mem] + am.effect[k]*bcap.am + bcap.p*reencounter[k] + b.dist[group.mem]*d[1:n.traps])*traps.avail[1:n.traps,group.mem])          
     cp[1:n.traps,k] <- lp[1:n.traps]/(1+sum(lp[1:n.traps]))
     cp[n.traps+1,k] <- 1-sum(cp[1:n.traps,k])  # last cell = not captured
    }
    pstar<- 1-prod(cp[ (n.traps+1) ,1:last.cap]) # 1-probablity missed on every occasion
   }#if(z==1)

  return(pstar)
  }
)


## function to derive group specific abundance
getGroupAbundance <- nimbleFunction(
  run = function(z = double(1), group.mem=double(1), n.groups=double(0), M=double(0)){
   returnType(double(1))
   G.N <- rep(0, n.groups)
   group.out <- group.mem[1:M]*z[1:M]
   for(j in 1:n.groups){
    G.N[j] <- sum( (group.out[1:M]==j)*1 )
   } 
  return(G.N)
  }
)


## data
Y0 <- rep(0, peromyscus.data$M)
nind <- sum(!is.na(group.mem))
nz <- sum(is.na(group.mem))
peromyscus.data$z <- c( rep(1, nind), rep(NA, nz) )
zst<-c( rep(NA, sum(!is.na(group.mem))), rep(0, sum(is.na(group.mem))) )
gst <- group.mem
gst[is.na(group.mem)]<- sample(1:24,sum(is.na(group.mem)),replace=TRUE)
gst[!is.na(group.mem)]<-NA
inits <- list (z=zst,group.mem=gst,psi=runif(1,.9,1),b.site=runif(8),sigma.site=runif(1),b.seas=runif(2),b.thin=runif(1),b.fire=runif(1),cent=S.st,bgroup.p=runif(24),int.p=runif(1),sigma.p=runif(1),bcap.p=runif(1),bcap.am=runif(1),b.dist=runif(24),int.dist=runif(1),sigma.dist=runif(1))              

# Bundle data
jags.data <- peromyscus.data
nim.data <- list(Ycat = peromyscus.data$Ycat,group.mem = peromyscus.data$group.mem, traps.avail= peromyscus.data$traps.avail, z=peromyscus.data$z, Y0=Y0) 

constants <- list(nind=nind, n.groups = peromyscus.data$n.groups, n.traps= peromyscus.data$n.traps, 
                 M= peromyscus.data$M, last.cap= peromyscus.data$last.cap, n.sites= peromyscus.data$n.sites, 
season=peromyscus.data$season,site=peromyscus.data$site, thin=peromyscus.data$thin,fire=peromyscus.data$fire, 
Yl = peromyscus.data$Yl, Yu = peromyscus.data$Yu,Xl = peromyscus.data$Xl,Xu = peromyscus.data$Xu,trap.locs = peromyscus.data$trap.locs,
reencounter = peromyscus.data$reencounter, trap.space= peromyscus.data$trap.space,
am.effect= peromyscus.data$am.effect)

parameters <- c("psi","sigma.site","b.season","b.fire","b.thin","int.p","sigma.p","bcap.p","bcap.am","int.dist","sigma.dist","N.tot","G.N")

ni <- 50
nt <- 1
nb <- 20
nc <- 1
adaptInterval <- 1000

#out.1 <- jags(jags.data, inits, parameters, "replicated_scr_model.txt", n.chains = nc, n.burnin = nb, n.thin = nt, n.iter = ni, n.adapt=10)
start.time<-Sys.time()

Rmodel <- nimbleModel(code=replicated_scr_model, data=nim.data, constants=constants, calculate=F, check=F, inits=inits)
conf <- configureMCMC(Rmodel,monitors=parameters, control = list(adaptInterval = adaptInterval ), thin=nt, useConjugacy = FALSE)  
 # categorical samplers can sometimesbe slow. could try slice as an alternative. This is not the ultimate speed issue though.
 #conf$printSamplers("group.mem")
 #conf$removeSampler("group.mem")
 #for(i in (sum(!is.na(group.mem))+1):peromyscus.data$M){
 # conf$addSampler(target = paste("group.mem[",i,"]", sep=""), type = "slice",control = list(adaptInterval = adaptInterval)) 
 #}

Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)#, showCompilerOutput = TRUE)
samp.time<-Sys.time()
out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-samp.time
end.time-start.time


plot(out[,"N.tot"])
