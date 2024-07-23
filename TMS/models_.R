# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module
library(runjags)
library(rjags)
load.module("wiener")


# setwd("/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Clases/Clases_UDD/Rafael/scripts")
setwd("/Volumes/GoogleDrive-112808863907079649330/Mi\ unidad/Working\ papers/Martinez_theta/OLDs/DATA")
setwd("/Volumes/GoogleDrive/Mi\ unidad/Working\ papers/Martinez_theta/OLDs/DATA")

source("HDIofMCMC.r") 

# Draw random samples with JAGS
data <- read.table("COR.txt", header=TRUE,sep="\t")
names(data)
dime =dim(data)
# [1] "rt"       "est"      "resp"     "laten"    "good"     "Sub"      "TMS"      "TMS_SITE" "nnb"      "seq"     

# creara relevan variables 
data$Dtms51 = as.numeric(data$TMS ==51);
data$Dtms52 = as.numeric(data$TMS==52);
data$Dtms50 = as.numeric(data$TMS==50);
data$Dtms=data$Dtms51+data$Dtms52;
# levels(data$TMS_SITE)
data$DtmsB = as.numeric(data$TMS_SITE==" B") # ojo CON ESE ESPACIO !!!
data$DtmsA = as.numeric(data$TMS_SITE==" A") # ojo CON ESE ESPACIO !!!
data$Uno= as.numeric(data$seq==1)

data$error = (data$rt<0 & data$est==10) + (data$rt>0 & data$est==21) 
data$perror = c(0, data$error[1:(dime[1]-1)])


su =levels(factor(data$Sub))
data$Dsu <- 0 
n=0;
for (e  in su) {
    n=n+1;
    data$Dsu =  data$Dsu + (data$Sub==e)*n
}
# correlative id for subjects 
sort(unique(data$Dsu))
(Nsubj=length(unique(data$Dsu)))

# only correct, go 
dataG = data[data$rt<1000 & data$rt>10 & data$est==10,]
dataG$Dseq_n = dataG$seq

Q=0.25
dataG$seq = 1-((1-Q)^(dataG$seq-1)) # espectative 
dime =dim(dataG)
dataG$pC = c(0, dataG$Dseq_n[1:(dime[1]-1)]>=dataG$Dseq_n[2:(dime[1])] );
  
  
Ntotal = length(dataG$rt)

# data for JASG
dat <- dump.format(list(Drt=dataG$rt, nT=Ntotal, Nsubj=Nsubj, idSub = dataG$Dsu, Dseq=dataG$seq, Dseq_n=dataG$Dseq_n, Uno=dataG$Uno, pError=dataG$perror, Dtms=dataG$Dtms, Dtms51=dataG$Dtms51, DtmsA=dataG$DtmsA))


# Visualize the data
hist(dataG$rt, breaks=50)

# prepare the data for JAGS


# Initialize chains
inits1 <- dump.format(list(alph=20,ta=8,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=9,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-1,mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,mubeta10=-1,mubeta11=0,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=150,ta=5,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=9,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list(alph=180,ta=6,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=9,mubeta1=-1,mubeta2=-1,mubeta3=1,mubeta4=1,mubeta5=-1,mubeta6=1,mubeta7=-1,mubeta8=-1,mubeta9=1,mubeta10=-1,mubeta11=-1,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8','mubeta9', 'mubeta10','mubeta11',"alph","ta","bet","delt","mu.Q","deviance")

# Run the function that fits the models using JAGS
results <- run.jags(model="model_mixed.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5, #modules=c("wiener"), 
                    method=c("parallel"))

results

chains = rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC




# Run the function that fits the models using JAGS
results_robust <- run.jags(model="model_mixed_robust.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

results_robust

chains = rbind(results_robust$mcmc[[1]], results_robust$mcmc[[2]], results_robust$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

save.image(file = "modelo_mixt_robust_sn.RData") 

# Run the function that fits the models using JAGS
results_DDM <- run.jags(model="model_mixed_DDM.txt",
                           monitor=monitor, data=dat, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=1000, sample=1000, thin=5,
                           modules=c("wiener"), method=c("parallel"))

results_DDM

chains = rbind(results_DDM$mcmc[[1]], results_DDM$mcmc[[2]], results_DDM$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

save.image(file = "modelo_mixt_robust_DDM_sn.RData") 

# Run the function that fits the models using JAGS
results_DDMmx <- run.jags(model="model_mixed_DDM_mx.txt",
                        monitor=monitor, data=dat, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

results_DDMmx

chains = rbind(results_DDMmx$mcmc[[1]], results_DDMmx$mcmc[[2]], results_DDMmx$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC

# Run the function that fits the models using JAGS
results_DDMmx2 <- run.jags(model="model_mixed_DDM_mx2.txt",
                        monitor=monitor, data=dat, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

results_DDMmx2

chains = rbind(results_DDMmx2$mcmc[[1]], results_DDMmx2$mcmc[[2]], results_DDMmx2$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC




save.image(file = "modelo_mixt_robust_DDMmx_sn.RData") 


#chains = rbind(results_robus$mcmc[[1]], results_robus$mcmc[[2]], results_robus$mcmc[[3]])
chains = rbind(results_DDM$mcmc[[1]], results_DDM$mcmc[[2]], results_DDM$mcmc[[3]])


(pval = min(mean(chains[,"mubeta5"]<0)*2, mean(chains[,"mubeta5"]>0)*2))
(pval = min(mean(chains[,"mubeta6"]<0)*2, mean(chains[,"mubeta6"]>0)*2))
(pval = min(mean(chains[,"mubeta7"]<0)*2, mean(chains[,"mubeta7"]>0)*2))
(pval = min(mean(chains[,"mubeta8"]<0)*2, mean(chains[,"mubeta8"]>0)*2))




cadena =  c(chains[,"mubeta8"],chains[,"mubeta7"],chains[,"mubeta6"],chains[,"mubeta5"])
plotData= data.frame(cadena)
plotData$beta = c(rep("8TMS:Th:A",3000),rep("7TMS:A",3000),rep("6TMS:Th",3000),rep("5TMS",3000))


### plot 
library(latex2exp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)


data_summary <- function(x) {
  m <- median(x)
  hdi = HDIofMCMC(x , credMass=0.95)
  ymin <- hdi[1]
  ymax <-  hdi[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#tikzDevice::tikz(file = "./taus_TMS.tex", width = 3, height = 3)

ggplot(plotData, aes(y=cadena, x = beta, color=beta, fill = beta)) + 
  geom_violin( )+
  #theme_minimal()+ 
  geom_hline(yintercept=0,col="red") +
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00", "#42C01A"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00", "#42C01A"))+
  ggtitle("Reaction time slowing by spectative of nogo stimuli")+  #  
  ylab("beta estimated (ms)")+
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  # geom_boxplot(width=0.1) +
  theme(legend.position = "none") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")

#dev.off()
