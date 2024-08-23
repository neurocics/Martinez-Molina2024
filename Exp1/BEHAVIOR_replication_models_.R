# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module
library(runjags)
library(rjags)
library(loo)
load.module("wiener")
library(lme4)
library(rstatix)
library(ggplot2)
# setwd("/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Clases/Clases_UDD/Rafael/scripts")
#setwd("/Volumes/GoogleDrive-112808863907079649330/Mi\ unidad/Working\ papers/Martinez_theta/OLDs/DATA/control_noTMS")
#setwd("/Volumes/GoogleDrive/Mi\ unidad/Working\ papers/Martinez_theta/OLDs/DATA/control_noTMS")
setwd("/Users/pablobilleke/Library/CloudStorage/OneDrive-udd.cl/Working\ papers/Martinez_theta/OLDs/DATA/control_noTMS")
source("../HDIofMCMC.r") 

pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST=""}}
    }  
}



# Draw random samples with JAGS
data <- read.table("replication.txt", header=TRUE,sep="\t")

names(data)
dime =dim(data)
# [1]  "rt"    "est"   "resp"  "laten" "good"  "SU" "seq"     

# for new scrot with BIDS compativility
data$est=(data$value==10)*1 + (data$value==21)*2



# create relevan variable 

data$Uno= as.numeric(data$seq==1)
data$error = (data$rt<0 & data$est==1) + (data$rt>0 & data$est==2) 
data$perror = c(0, data$error[1:(dime[1]-1)])

data$n2 = data$seq<3 & data$rt>40 & data$rt<1500
data$m2 = data$seq>=3 & data$rt>40 &data$rt<1500 & data$seq<6
data$faster = data$rt <  median(data$rt [data$rt>40 &data$rt<1500])


Dml = aggregate(faster~ n2, FUN=mean,data=data[(data$n2==1)|(data$m2==1) & data$rt>40 & data$rt<1500,] )


sum(data$n2) 
sum(data$n2) / (sum(data$m2) +sum(data$n2))

Dm = aggregate(cbind(rt,error)~ m2+SU, FUN=mean,data=data[(data$n2==1)|(data$m2==1) & data$rt>40,] )
wilcox.test(error~m2,data=Dm,paired=T)
t.test(error~m2,data=Dm,paired=T)
aggregate(error~m2,data=Dm, FUN=mean)
wilcox_effsize(Dm, error~m2,paired=T )

Dm = aggregate(rt~ m2+SU, FUN=median,data=data[(data$n2==1)|(data$m2==1) & data$rt>40,] )
Rwil = wilcox.test(rt~m2,data=Dm,paired=T)
t.test(rt~m2,data=Dm,paired=T)
wilcox_effsize(Dm, rt~m2,paired=T )
aggregate(rt~m2,data=Dm, FUN=mean)



MP_p <- ggplot(Dm, aes(x=m2, y=rt, col=m2)) + 
  geom_boxplot(linewidth = 0.5, width=0.4,outlier.shape = NA) + theme_classic()+
  geom_line(aes(x= (m2 - (m2==1)*0.2  + (m2==0)*0.2)+ 1 , y=rt, group = SU),color="grey",linetype=2,size=0.25)+
  geom_jitter( position=position_jitter(0.1),size = 1,alpha=0.7) +
  ylab("Reaction Time (ms)") +
  xlab("Go Stimui") +
  scale_x_discrete(labels=c("TRUE" = "late", "FALSE" = "early"))+ 
  
  theme(legend.position = "none") +  
  geom_segment(aes(x = 1.3, y = 525, xend = 1.7, yend = 525), color = "black") +
  annotate(geom="text", x=1.5, y=535, label=pv2str(Rwil $p.value),
           color="black") 











su =levels(as.factor(data$SU))
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
dataG = data[data$rt<1500 & data$rt>40 & data$est==1,]
dataG$Dseq_n = dataG$seq
Q=0.25 # a prior probability of occurrence of a ngo stimuli 
dataG$seq = 1-((1-Q)^(dataG$seq-1)) # expectation [[ EQUATION exp  ]]
dime =dim(dataG)
dataG$pC = c(0, dataG$Dseq_n[1:(dime[1]-1)]>=dataG$Dseq_n[2:(dime[1])] );

# first trial for each subject 
dime =dim(data)
data$T1 = 1-c(0, data$Dsu[1:(dime[1]-1)]==data$Dsu[2:(dime[1])])

# only correct, go trial for RTs model 
Trt = which(data$rt<1500 & data$rt>40 & data$est==1 & !data$T1)

#
NGO = as.numeric(data$est==2)



Ntotal = length(dataG$rt)

# data for JASG, only correct go trials 
dat <- dump.format(list(Drt=dataG$rt, nT=Ntotal, Nsubj=Nsubj, idSub = dataG$Dsu, Dseq=dataG$Dseq_n, Dseq_n=dataG$Dseq_n, Uno=dataG$Uno, pError=dataG$perror
                        ))

Ntotal = length(data$rt)

# data for JASG, all trials, ans index for correct go trials (Trt) 
datall <- dump.format(list(Drt=data$rt, nT=Ntotal, Nsubj=Nsubj, idSub = data$Dsu, Dseq=data$seq, Dseq_n=data$seq, Uno=data$Uno, pError=data$perror, Trt=Trt,T1=data$T1
))


# Visualize the data
hist(dataG$rt, breaks=50)




# Initialize chains
inits1 <- dump.format(list(alph=200,ta=15,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=10,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-1,mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,mubeta10=-1,mubeta11=0,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=300,ta=20,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=10,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list(alph=250,ta=17,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=10,mubeta1=0.1,mubeta2=-.1,mubeta3=.1,mubeta4=.1,mubeta5=-.1,mubeta6=1,mubeta7=-1,mubeta8=-1,mubeta9=1,mubeta10=-1,mubeta11=-1,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','nor',"ta","bet","delt","mu.LR","deviance")


resultsl0 <- run.jags(model="model_null_l.txt",
                                              monitor=monitor, data=datall, n.chains=3, 
                                              inits=c(inits1, inits2, inits3), plots = TRUE,
                                              burnin=2000, sample=1000, thin=5,
                                              modules=c("wiener"), method=c("parallel"))

resultsl0
chains = rbind(resultsl0$mcmc[[1]], resultsl0$mcmc[[2]], resultsl0$mcmc[[3]])
DIC.M0lgn = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M0lgn

para_loo <- extend.jags(resultsl0,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M0lgn = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M0lgn <- loo(para_loo_M0lgn)
waic.M0lgn <- waic(para_loo_M0lgn)


# Run the function that fits the models using JAGS
# model_mixed.txt  --> a mixed linear model with lognormal function a a link function for RT
#                  --> Sequence of go trials adjusted as expectation  that this secuence occurer , that if     
resultsl <- run.jags(model="model_mixed_l.txt",
                    monitor=monitor, data=datall, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

summary(resultsl)

chains = rbind(resultsl$mcmc[[1]], resultsl$mcmc[[2]], resultsl$mcmc[[3]])
DIC.M1lgn = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M1lgn

# 79572 datall - 79958.02 (dat)

para_loo <- extend.jags(resultsl,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M1lgn = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M1lgn <- loo(para_loo_M1lgn)
waic.M1lgn <- waic(para_loo_M1lgn)


# Run the function that fits the models using JAGS
# model_mixed.txt  --> a mixed linear model with lognormal function a a link function for RT
#                  --> Sequence of go trials adjusted as expectation  that this secuence occurer , that if     
results <- run.jags(model="model_mixed.txt",
                    monitor=monitor, data=datall, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

results

chains = rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
DIC.M2lgn = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M2lgn

#  79923.24 (dat)
# [1] 79528.88 (datall)

para_loo <- extend.jags(results,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M2lgn = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M2lgn <- loo(para_loo_M2lgn)
waic.M2lgn <- waic(para_loo_M2lgn)


# Run the function that fits the models using JAGS
results_mLR <- run.jags(model="model_mixed_LR.txt",
                    monitor=monitor, data=datall, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=10000, sample=1000, thin=10,
                    modules=c("wiener"), method=c("parallel"))

results_mLR 

chains = rbind(results_mLR $mcmc[[1]], results_mLR$mcmc[[2]], results_mLR$mcmc[[3]])
DIC.M3lgn = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3lgn
# 79446.85

para_loo <- extend.jags(results_mLR,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M3lgn = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M3lgn <- loo(para_loo_M3lgn)
waic.M3lgn <- waic(para_loo_M3lgn)

# Run the function that fits the models using JAGS
results0_DDM <- run.jags(model="model_null_DDM.txt",
                        monitor=monitor, data=datall, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

results0_DDM

chains = rbind(results0_DDM$mcmc[[1]], results0_DDM$mcmc[[2]], results0_DDM$mcmc[[3]])
DIC.M0ddm = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M0ddm

para_loo <- extend.jags(results0_DDM,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M0ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M0ddm <- loo(para_loo_M0ddm)
waic.M0ddm <- waic(para_loo_M0ddm)




# Run the function that fits the models using JAGS
results_DDM <- run.jags(model="model_mixed_DDM.txt",
                    monitor=monitor, data=datall, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=2000, sample=1000, thin=5,
                    modules=c("wiener"), method=c("parallel"))

summary(results_DDM)

chains = rbind(results_DDM$mcmc[[1]], results_DDM$mcmc[[2]], results_DDM$mcmc[[3]])
DIC.M1ddm = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M1ddm
#Drt[514]

# 
#  79387.42 (datall)

para_loo <- extend.jags(results_DDM,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M1ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M1ddm <- loo(para_loo_M1ddm)
waic.M1ddm <- waic(para_loo_M1ddm)


# Run the function that fits the models using JAGS
results_DDMex <- run.jags(model="model_mixed_DDM_ex.txt",
                        monitor=monitor, data=datall, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

summary(results_DDMex)

chains = rbind(results_DDMex$mcmc[[1]], results_DDMex$mcmc[[2]], results_DDMex$mcmc[[3]])
DIC.M2ddm = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M2ddm

# 
# [1] 79581 79357.07[datall]

para_loo <- extend.jags(results_DDMex,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M2ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M2ddm <- loo(para_loo_M2ddm)
waic.M2ddm <- waic(para_loo_M2ddm)


# Run the function that fits the models using JAGS
results_DDM_LR <- run.jags(model="model_DDM_LR.txt",
                        monitor=monitor, data=datall, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=10000, sample=1000, thin=10,
                        modules=c("wiener"), method=c("parallel"))

results_DDM_LR 
results_DDM_LR  <- extend.jags(results_DDM_LR,#drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, 
                               sample = 1000, adapt=0,
                        burnin=1000, #summarise=FALSE, 
                        method="parallel")



chains = rbind(results_DDM_LR $mcmc[[1]], results_DDM_LR $mcmc[[2]], results_DDM_LR $mcmc[[3]])
DIC.M3ddm  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3ddm 


save.image("data_repli.Rdata")
#
rlps = extend.jags(results_DDM_LR,drop.monitor = results_DDM_LR$monitor,add.monitor = c("LR"),burnin=0, sample=200,adapt = 0)



rlpss =summary(rlps)
for (s in 1:max(data$Dsu)){
  data$LRmedian[data$Dsu==s] = rlpss[s,2]
  data$LRmean[data$Dsu==s] = max(0.0015,rlpss[s,4])
  data$LRmode[data$Dsu==s] = rlpss[s,6]
}


hist(data$LRmean)
# [1] 3510.939 (sin Q con Exp)  
#  3699 (exp sinQ sin pC)

# data for JASG, all trials, ans index for correct go trials (Trt) 
datall <- dump.format(list(Drt=data$rt, nT=Ntotal, Nsubj=Nsubj,
                           idSub = data$Dsu, Dseq=data$seq, Dseq_n=data$seq, Uno=data$Uno, 
                           pError=data$perror, Trt=Trt,T1=data$T1,
                           LRi=data$LRmean
))


# Run the function that fits the models using JAGS
results_DDM_LRfix <- run.jags(model="model_DDM_LRfix.txt",
                           monitor=monitor, data=datall, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=10000, sample=1000, thin=10,
                           modules=c("wiener"), method=c("parallel"))

results_DDM_LRfix 

chains = rbind(results_DDM_LRfix $mcmc[[1]], results_DDM_LRfix $mcmc[[2]], results_DDM_LRfix $mcmc[[3]])
DIC.M3ddm  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3ddm 




# 79285.28

para_loo <- extend.jags(results_DDM_LRfix,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M3ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M3ddm <- loo(para_loo_M3ddm)
waic.M3ddm <- waic(para_loo_M3ddm)


# Run the function that fits the models using JAGS
results_DDM_LR_c <- run.jags(model="model_DDM_LR_c.txt",
                           monitor=monitor, data=datall, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=10000, sample=1000, thin=10,
                           modules=c("wiener"), method=c("parallel"))

results_DDM_LR_c 

chains = rbind(results_DDM_LR_c$mcmc[[1]], results_DDM_LR_c$mcmc[[2]], results_DDM_LR_c$mcmc[[3]])
DIC.M3ddm_c = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3ddm_c
#  79330.65

para_loo <- extend.jags(results_DDM_LR_c,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")



para_loo <- extend.jags(results_DDM_LR_c,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M3ddm_c = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M3ddm_c <- loo(para_loo_M3ddm_c)
waic.M3ddm_c <- waic(para_loo_M3ddm_c)


setwd("~/Documents/GitHub/Martinez-Molina2024/Exp1")
# Run the function that fits the models using JAGS
results_DDM_Q_LR <- run.jags(model="model_Q_LR_DDM.txt",
                             monitor=monitor, data=datall, n.chains=3, 
                             inits=c(inits1, inits2, inits3), plots = TRUE,
                             burnin=10000, sample=1000, thin=10,
                             modules=c("wiener"), method=c("parallel"))

results_DDM_Q_LR 
results_DDM_Q_LR2 <- extend.jags(results_DDM_Q_LR,drop.monitor = results_DDM_Q_LR$monitor,
                             add.monitor = results_DDM_Q_LR$monitor,
                             burnin=10000, sample=1000, thin=10)

results_DDM_Q_LR2 



chains = rbind(results_DDM_Q_LR2$mcmc[[1]], results_DDM_Q_LR2$mcmc[[2]], results_DDM_Q_LR2$mcmc[[3]])
DIC.M4ddm  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M4ddm 

# 79350.75

para_loo <- extend.jags(results_DDM_Q_LR2,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M4ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M4ddm <- loo(para_loo_M4ddm)
waic.M4ddm <- waic(para_loo_M4ddm)




# Compare the two models
loo_compare(loo.M1ddm,loo.M1lgn,loo.M2ddm,loo.M2lgn)

loo_compare(loo.M3ddm,
            loo.M0ddm,
            loo.M1ddm,
            loo.M2ddm,
            loo.M4ddm)*2


c(DIC.M3ddm,DIC.M3lgn,DIC.M0lgn,DIC.M0ddm,DIC.M1lgn,DIC.M1ddm,DIC.M2lgn,DIC.M2ddm)-DIC.M3ddm


save.image(file="data_repli.22.08.2024.RData")

##3  Figure 
chains = rbind(results_DDM_LRfix $mcmc[[1]], results_DDM_LRfix $mcmc[[2]], results_DDM_LRfix $mcmc[[3]])
cadena =  c(chains[,"mubeta1"]/1,
            chains[,"mubeta2"]/1,
            chains[,"mubeta3"]/1
            #chains[,"nor"],
            #chains[,"ta"]/200,
            #chains[,"mu.LR"]
)
plotData= data.frame(cadena)
plotData$beta = c(rep("Ex",3000),
                  rep("pCE",3000),
                  rep("pOE",3000)
                  #rep("Driff",3000),
                  #rep("tau",3000),
                  #rep("LR",3000)
)

mean(plotData$cadena[plotData$beta=="Ex"])*1000
HDIofMCMC(plotData$cadena[plotData$beta=="Ex"])*1000
ppoe = min(mean(plotData$cadena[plotData$beta=="pOE"]>0),mean(plotData$cadena[plotData$beta=="pOE"]<0))*2
ppce = min(mean(plotData$cadena[plotData$beta=="pCE"]>0),mean(plotData$cadena[plotData$beta=="pCE"]<0))*2
pexp = min(mean(plotData$cadena[plotData$beta=="Ex"]>0),mean(plotData$cadena[plotData$beta=="Ex"]<0))*2

### plot 
library(latex2exp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)

pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST=""}}
    }  
}


data_summary <- function(x) {
  m <- median(x)
  hdi = HDIofMCMC(x , credMass=0.95)
  ymin <- hdi[1]
  ymax <-  hdi[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#tikzDevice::tikz(file = "./taus_TMS.tex", width = 3, height = 3)

PD = ggplot(plotData, aes(y=cadena, x = beta, color=beta, fill = beta)) + 
  geom_violin( )+
  #theme_minimal()+ 
  geom_hline(yintercept=0,col="red") +
  scale_fill_manual(values=c("#77ac30ff","#808080ff","#808080ff"))+
  scale_color_manual(values=c("#77ac30ff","#808080ff","#808080ff"))+
   #ggtitle("Reaction time slowing")+  #  
  ylab("Posterior distribution (a.u.)")+
  xlab("Model parameters")+
  theme(legend.position = "none") +
  coord_flip()+
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=10, x=3.3, label=pv2str(ppoe),
           color="black") +
  annotate(geom="text", y=25, x=2.3, label=pv2str(ppce),
           color="black") +
  annotate(geom="text", y=65, x=1.3, label=pv2str(pexp),
           color="black") 


# Model comparison
DICS=c(loo.M0ddm$estimates[3,1],#DIC.M0lgn, 
       loo.M1ddm$estimates[3,1],#
       loo.M2ddm$estimates[3,1],#
       loo.M3ddm$estimates[3,1],#
       loo.M4ddm$estimates[3,1],#
       DIC.M0ddm, 
       DIC.M1ddm,
       DIC.M2ddm,
       DIC.M3ddm,
       DIC.M4ddm)
#DICS=c(loo.M0lgn$estimates[3], loo.M1lgn$estimates[3],loo.M2lgn$estimates[3],loo.M3lgn$estimates[3] ,
       #  loo.M0ddm$estimates[3], loo.M1ddm$estimates[3],loo.M2ddm$estimates[3],loo.M3ddm$estimates[3])

Link=factor(c("LOOIC","LOOIC","LOOIC","LOOIC","LOOIC",
              "DIC","DIC","DIC","DIC","DIC"))
df <- data.frame(modelos = c("M0", "M1", "M2","M3(seq)","M3(Q)","M0", "M1", "M2","M3(seq)","M3(Q)"),
                 Link = Link,
                 dics = c(DICS[1:5]- min(DICS[1:5])+1  ,DICS[6:10]- min(DICS[6:10])+1  ) )
DIC = ggplot(df) +
  geom_col(aes(x=modelos, y=dics,fill=Link),position = position_dodge())+
  theme(legend.position=c(0.7, 0.9)) +
  labs(#title="Model Comparison by Metric", 
    x="Models", 
    y="Adjustment Difference", 
    fill="Metric")
library(latex2exp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)
library(patchwork)

layout <- "
AB
CB
"
#p1 + p2 + p3 + p4 + 
#  plot_layout(design = layout)

MP_p + PD+ DIC+plot_layout(design = layout) +plot_annotation(tag_levels = list(c('A','C','B')))

#save.image(file="all_behave_noTMS_07062023.RData")
