# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module
library(runjags)
library(rjags)
load.module("wiener")
# setwd("/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Clases/Clases_UDD/Rafael/scripts")
setwd("/Volumes/GoogleDrive-112808863907079649330/Mi\ unidad/Working\ papers/Martinez_theta/OLDs/DATA/control_noTMS")
source("../HDIofMCMC.r") 

# Draw random samples with JAGS
data <- read.table("COR.txt", header=TRUE,sep="\t")

names(data)
dime =dim(data)
# [1]  "rt"    "est"   "resp"  "laten" "good"  "SU" "seq"     

data$seq[1] = 1
for (t in 2:dime[1]) {
  if ( data$est[t]==1  && data$SU[t-1]==data$SU[t] ){
  data$seq[t] = data$seq[t-1] +1 }else{
    data$seq[t] = 0 
  }
}

# creara relevan variables 

data$Uno= as.numeric(data$seq==1)
data$error = (data$rt<0 & data$est==1) + (data$rt>0 & data$est==2) 
data$perror = c(0, data$error[1:(dime[1]-1)])
su =levels(data$SU)
data$Dsu <- 0 
n=0;
for (e  in su) {
  n=n+1;
  data$Dsu =  data$Dsu + (data$SU==e)*n
}
# correlative id for subjects 
sort(unique(data$Dsu))
(Nsubj=length(unique(data$Dsu)))



# only correct, go 
dataG = data[data$rt<1000 & data$rt>40 & data$est==1,]
dataG$Dseq_n = dataG$seq
Q=0.25
dataG$seq = 1-((1-Q)^(dataG$seq-1)) # espectative 
dime =dim(dataG)
dataG$pC = c(0, dataG$Dseq_n[1:(dime[1]-1)]>=dataG$Dseq_n[2:(dime[1])] );

# first trial for each subject 
dime =dim(data)
data$T1 = 1-c(0, data$Dsu[1:(dime[1]-1)]==data$Dsu[2:(dime[1])])

# only correct, go 
Trt = which(data$rt<1500 & data$rt>40 & data$est==1 & !data$T1)
pTrt = (data$rt<1500 & data$rt>40 & data$est==1 & !data$T1) 
pTrt = c(0,pTrt[1:(dime[1]-1)])

#
NGO = as.numeric(data$est==2)
Tng=which( NGO==1 & pTrt==1 )
Tg = which(data$est==1 & pTrt==1)





Ntotal = length(dataG$rt)

# data for JASG
dat <- dump.format(list(Drt=dataG$rt, nT=Ntotal, Nsubj=Nsubj, idSub = dataG$Dsu, Dseq=dataG$Dseq_n, Dseq_n=dataG$Dseq_n, Uno=dataG$Uno, pError=dataG$perror
                        ))

Ntotal = length(data$rt)

# data for JASG
datall <- dump.format(list(Drt=data$rt, nT=Ntotal, Nsubj=Nsubj, idSub = data$Dsu, Dseq=data$seq, Dseq_n=data$seq, Uno=data$Uno, pError=data$perror, Trt=Trt,T1=data$T1
))


# Visualize the data
hist(dataG$rt, breaks=50)

# prepare the data for JAGS


# Initialize chains
inits1 <- dump.format(list(alph=200,ta=15,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=10,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-1,mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,mubeta10=-1,mubeta11=0,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=300,ta=20,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=10,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list(alph=400,ta=10,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=10,mubeta1=-1,mubeta2=-1,mubeta3=1,mubeta4=1,mubeta5=-1,mubeta6=1,mubeta7=-1,mubeta8=-1,mubeta9=1,mubeta10=-1,mubeta11=-1,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8','mubeta9', 'mubeta10','mubeta11',"alph","ta","bet","delt","mu.Q","deviance")

# Run the function that fits the models using JAGS
results <- run.jags(model="model_mixed.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

results

chains = rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# 79791.98

# Run the function that fits the models using JAGS
results_mLR <- run.jags(model="model_mixed_LR.txt",
                    monitor=monitor, data=datall, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=10000, sample=1000, thin=10,
                    modules=c("wiener"), method=c("parallel"))

results_mLR 

chains = rbind(results_mLR $mcmc[[1]], results_mLR $mcmc[[2]], results_mLR $mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# 79457.38


# Run the function that fits the models using JAGS
results_DDM <- run.jags(model="model_mixed_DDM.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

results_DDM

chains = rbind(results_DDM$mcmc[[1]], results_DDM$mcmc[[2]], results_DDM$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# Boundering 
# [1] 79601

# Run the function that fits the models using JAGS
results_DDM_LR <- run.jags(model="model_DDM_LR.txt",
                        monitor=monitor, data=datall, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=10000, sample=1000, thin=10,
                        modules=c("wiener"), method=c("parallel"))

results_DDM_LR 

chains = rbind(results_DDM_LR $mcmc[[1]], results_DDM_LR $mcmc[[2]], results_DDM_LR $mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# 79285.67









# Run the function that fits the models using JAGS
results_DDMdr <- run.jags(model="model_mixed_DDM_dr.txt",
                        monitor=monitor, data=dat, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

results_DDMdr

chains = rbind(results_DDMdr$mcmc[[1]], results_DDMdr$mcmc[[2]], results_DDMdr$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# Drif
# [1] Con UNO 
# [1] 79458.49 sin UNO con mas ajuste de cadenas 

# Run the function that fits the models using JAGS
results_DDMmx <- run.jags(model="model_mixed_DDM_mix.txt",
                          monitor=monitor, data=dat, n.chains=3, 
                          inits=c(inits1, inits2, inits3), plots = TRUE,
                          burnin=2000, sample=1000, thin=5,
                          modules=c("wiener"), method=c("parallel"))

results_DDMmx

chains = rbind(results_DDMmx$mcmc[[1]], results_DDMmx$mcmc[[2]], results_DDMmx$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# Drif:Seq  //  baoudering:Eror
# [1]  79013.58 aka Con UNO 
# [1] [1] 79099.59 sin UNO

# Run the function that fits the models using JAGS
results_DDMmx2 <- run.jags(model="model_mixed_DDM_mix2.txt",
                          monitor=monitor, data=dat, n.chains=3, 
                          inits=c(inits1, inits2, inits3), plots = TRUE,
                          burnin=2000, sample=1000, thin=5,
                          modules=c("wiener"), method=c("parallel"))

results_DDMmx2

chains = rbind(results_DDMmx2$mcmc[[1]], results_DDMmx2$mcmc[[2]], results_DDMmx2$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# Drif:error  //  baoudering:Seq
# [1] Con UNO 
# [1] 79082.69 sin UNO

save.image("models_.RData")
