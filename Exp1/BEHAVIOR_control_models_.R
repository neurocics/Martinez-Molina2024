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
library(coda)

setwd("~/Documents/GitHub/Martinez-Molina2024/Exp1/")
source("../HDIofMCMC.r") 

pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST=""}}
    }  
}



# Draw random samples with JAGS
# in COR.txt all de trails of behavioral experiments 
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
  geom_boxplot(size = 0.5, width=0.4) + theme_classic()+
  geom_line(aes(x= (m2 - (m2==1)*0.2  + (m2==0)*0.2)+ 1 , y=rt, group = SU),color="grey",linetype=2,size=0.25)+
  geom_jitter( position=position_jitter(0.1),size = 1,alpha=0.7) +
  ylab("Reaction Time") +
  xlab("Go Stimui") +
  scale_x_discrete(labels=c("TRUE" = "late", "FALSE" = "early"))+ 
  
  theme(legend.position = "none") +  
  geom_segment(aes(x = 1.3, y = 425, xend = 1.7, yend = 425), color = "black") +
  annotate(geom="text", x=1.5, y=435, label=pv2str(Rwil $p.value),
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
datall <- dump.format(list(Drt=data$rt, nT=Ntotal, Nsubj=Nsubj, 
                           idSub = data$Dsu, Dseq=data$seq, Dseq_n=data$seq,
                           Uno=data$Uno, pError=data$perror, 
                           Trt=Trt,T1=data$T1,ISI=data$ISI
))


# Visualize the data
hist(dataG$rt, breaks=50)




# Initialize chains
inits1 <- dump.format(list(alph=200,ta=15,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=10,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-1,
                           mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,
                           mubeta10=-1,mubeta11=0,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=300,ta=20,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=10,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list(alph=400,ta=10,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=10,mubeta1=-1,mubeta2=-1,mubeta3=1,mubeta4=1,mubeta5=-1,mubeta6=1,mubeta7=-1,mubeta8=-1,mubeta9=1,mubeta10=-1,mubeta11=-1,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','nor',
            "ta","bet","delt","mu.LR","deviance",
            "th","muthr")


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
para_stat <- extend.jags(results_DDMex,drop.monitor = monitor,add.monitor = monitor, sample = 1000, adapt=0,
                        burnin=0,  method="parallel")

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

chains = rbind(results_DDM_LR $mcmc[[1]], results_DDM_LR $mcmc[[2]], results_DDM_LR $mcmc[[3]])
DIC.M3ddm  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3ddm 

# 79285.28

para_loo <- extend.jags(results_DDM_LR,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M3ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M3ddm <- loo(para_loo_M3ddm)
waic.M3ddm <- waic(para_loo_M3ddm)


para_loo <- extend.jags(results_DDM_LR_c,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M3ddm_c = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M3ddm_c <- loo(para_loo_M3ddm_c)
waic.M3ddm_c <- waic(para_loo_M3ddm_c)



# M4 is the  paper's model M3(Q)
# Run the function that fits the models using JAGS
results_DDM_Q_LR <- run.jags(model="model_Q_LR_DDM.txt",
                             monitor=monitor, data=datall, n.chains=3, 
                             inits=c(inits1, inits2, inits3), plots = TRUE,
                             burnin=10000, sample=1000, thin=10,
                             modules=c("wiener"), method=c("parallel"))

results_DDM_Q_LR 

chains = rbind(results_DDM_Q_LR$mcmc[[1]], results_DDM_Q_LR$mcmc[[2]], results_DDM_Q_LR$mcmc[[3]])
DIC.M4ddm  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M4ddm 

# 79350.75

para_loo <- extend.jags(results_DDM_Q_LR,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M4ddm = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M4ddm <- loo(para_loo_M4ddm)
waic.M4ddm <- waic(para_loo_M4ddm)


#
# Run the function that fits the models using JAGS
results_lgn_Q_LR <- run.jags(model="model_Q_LR_lognorm.txt",
                             monitor=monitor, data=datall, n.chains=3, 
                             inits=c(inits1, inits2, inits3), plots = TRUE,
                             burnin=10000, sample=1000, thin=10,
                             modules=c("wiener"), method=c("parallel"))

results_lgn_Q_LR ## ensure adequate/good convergence of the chains
results_lgn_Q_LR2 <- extend.jags(results_lgn_Q_LR,check.stochastic = FALSE, sample = 1000, adapt=0,
                                 burnin=1000, summarise=FALSE, method="parallel")
results_lgn_Q_LR3 <- extend.jags(results_lgn_Q_LR2,check.stochastic = FALSE, sample = 1000, adapt=0,
                                 drop.monitor = results_lgn_Q_LR2$monitor, add.monitor = monitor,
                                 burnin=1000, summarise=FALSE, method="parallel")
results_lgn_Q_LR4 <- extend.jags(results_lgn_Q_LR3,check.stochastic = FALSE, sample = 1000, adapt=0,thin=10,
                                 drop.monitor = results_lgn_Q_LR2$monitor, add.monitor = monitor,
                                 burnin=10000, summarise=FALSE, method="parallel")
results_lgn_Q_LR5 <- extend.jags(results_lgn_Q_LR4,check.stochastic = FALSE, sample = 1000, adapt=0,thin=20,
                                 drop.monitor = results_lgn_Q_LR4$monitor, add.monitor = monitor,
                                 burnin=0, summarise=FALSE, method="parallel")
summary(results_lgn_Q_LR5)

chains = rbind(results_lgn_Q_LR5$mcmc[[1]], results_lgn_Q_LR5$mcmc[[2]], results_lgn_Q_LR5$mcmc[[3]])
DIC.M4lgn  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M4lgn 

# 

para_loo <- extend.jags(results_lgn_Q_LR5,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M4lgn = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M4lgn <- loo(para_loo_M4lgn)
waic.M4lgn <- waic(para_loo_M4lgn)

#save.image(file="all_behave_noTMS_15072024.RData")






#save.image(file="all_behave_noTMS_10072024.RData")

##3  Figure 
chains = rbind(results_DDM_LR $mcmc[[1]], results_DDM_LR $mcmc[[2]], results_DDM_LR $mcmc[[3]])
cadena =  c(chains[,"mubeta1"]/1000,
            chains[,"mubeta2"]/1000,
            chains[,"mubeta3"]/1000
            #chains[,"nor"],
            #chains[,"ta"]/200,
            #chains[,"mu.LR"]
)
plotData= data.frame(cadena)
plotData$beta = c(rep("Exp",3000),
                  rep("pCE",3000),
                  rep("pOE",3000)
                  #rep("Driff",3000),
                  #rep("tau",3000),
                  #rep("LR",3000)
)

mean(plotData$cadena[plotData$beta=="pOE"])*1000
HDIofMCMC(plotData$cadena[plotData$beta=="pOE"])*1000
ppoe = min(mean(plotData$cadena[plotData$beta=="pOE"]>0),mean(plotData$cadena[plotData$beta=="pOE"]<0))*2
ppce = min(mean(plotData$cadena[plotData$beta=="pCE"]>0),mean(plotData$cadena[plotData$beta=="pCE"]<0))*2
pexp = min(mean(plotData$cadena[plotData$beta=="Exp"]>0),mean(plotData$cadena[plotData$beta=="Exp"]<0))*2

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
  ylab("Posterior distribution")+
  xlab("Model parameters")+
  theme(legend.position = "none") +
  coord_flip()+
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=0.01, x=3.3, label=pv2str(ppoe),
           color="black") +
  annotate(geom="text", y=0.03, x=2.3, label=pv2str(ppce),
           color="black") +
  annotate(geom="text", y=0.06, x=1.3, label=pv2str(pexp),
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


layout <- "
AB
CB
"
#p1 + p2 + p3 + p4 + 
#  plot_layout(design = layout)

MP_p + PD+ DIC+plot_layout(design = layout) +plot_annotation(tag_levels = list(c('C','E','D')))




####.  control paper v2 
### # M5 is the  paper's cotrol model M4

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4',
            'nor',"ta","bet","delt","mu.LR",
            'th',
            "deviance")


# Run the function that fits the models using JAGS
results_DDM_LR_theta <- run.jags(model="model_DDM_LR_theta.txt",
                             monitor=monitor, data=datall, n.chains=3, 
                             inits=c(inits1, inits2, inits3), plots = TRUE,
                             burnin=10000, sample=1000, thin=10,
                             modules=c("wiener"), method=c("parallel"))

results_DDM_LR_theta

chains = rbind(results_DDM_LR_theta$mcmc[[1]], results_DDM_LR_theta$mcmc[[2]], results_DDM_LR_theta$mcmc[[3]])
DIC.M5  = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M5

mean(chains[,"th"])
HDIofMCMC(chains[,"th"])
mean(chains[,"th"]<0.5)

para_loo <- extend.jags(results_DDM_LR_theta,drop.monitor = monitor,add.monitor = c("Drtlog"),check.stochastic = FALSE, sample = 200, adapt=0,
                        burnin=0, summarise=FALSE, method="parallel")
para_loo_M5 = rbind(para_loo$mcmc[[1]], para_loo$mcmc[[2]],para_loo$mcmc[[3]])
loo.M5 <- loo(para_loo_M5)
waic.M5 <- waic(para_loo_M5)


##

#save.image(file="Model0_1_2_3_4_5_.RData")

###. valores predicho del modle M3ddm


para_muestars <- extend.jags(results_DDM_LR,
                             drop.monitor = results_DDM_LR$monitor,
                             add.monitor = c("y","ta","nor"),
                             check.stochastic = FALSE, sample = 100, adapt=0,
                             burnin=0, summarise=FALSE, method="parallel")


muestras <- as.mcmc.list(para_muestars)
muestras_matriz <- as.matrix(muestras)
# 300 7234

muestras_matriz_par <- muestras_matriz[,7233:7234]
muestras_matriz_par_m =colMeans(muestras_matriz_par) 
y_new = colMeans(muestras_matriz[,1:7232])  
y_new_all = muestras_matriz[,1:7232] 
library(brms)

#RTpre = numeric()
for (i in length(RTpre):7232){
d <- rwiener(n = 1, alpha = y_newall[i], tau = muestras_matriz_par_m[1] ,
        beta = 0.5, delta = muestras_matriz_par_m[2])
RTpre[i] = as.numeric(d[1])
}



RTpre2 = numeric()
for (i in 1:7232){
  # i = 1
  np <- sample(1:300, 1, replace = TRUE)
  d <- rwiener(n = 1, alpha = y_new_all[np,i], tau = muestras_matriz_par[np,1] ,
               beta = 0.5, delta = muestras_matriz_par[np,2])
  RTpre2[i] = as.numeric(d[1])
}










RTreal = data$rt[Trt]

data2 <- data.frame(
  Value = c(RTpre, RTreal),
  Data = rep(c("Predicted", "Observed"), each = length(RTpre2))
)

GA <- ggplot(data2, aes(x = Value, fill = Data)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 80) +
  labs(title = "RT Comparison",
       x = "Reaction Time (ms)",
       y = "Frequency") +
  theme_minimal() +
 # scale_fill_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  theme(legend.position = c(0.7, 0.85),  # Ajustar la posición de la leyenda (x, y) dentro del área del gráfico
        #legend.background = element_rect(fill = alpha('white', 0.5)),
        legend.background = element_blank(),  # Eliminar el fondo de la leyenda
        legend.key = element_blank())  # Hacer el fondo de la leyenda semi-transparente



# Calcular cuartiles
observed_quantiles <- quantile(RTreal, probs = seq(0.1,0.9,by=0.1))
predicted_quantiles <- quantile(RTpre, probs = seq(0.1,0.9,by=0.1))

# Crear un data frame para los cuartiles
quartile_data <- data.frame(
  Deciles = factor(c("p1", "p2","p3", "p4", "M","p6", "p7","p8", "p9"), 
                    levels = c("p1", "p2","p3", "p4", "M","p6", "p7","p8", "p9")),
  Observed = observed_quantiles,
  Predicted = predicted_quantiles
)



# Crear el gráfico
GC= ggplot(quartile_data, aes(x = Deciles)) +
  geom_point(aes(y = Observed, color = "Observed"), size = 4) +
  geom_point(aes(y = Predicted, color = "Predicted"), size = 4) +
  geom_line(aes(y = Observed, group = 1, color = "Observed")) +
  geom_line(aes(y = Predicted, group = 1, color = "Predicted")) +
  labs(title = "Quantile Plot",
       y = "Reaction Time (ms)",
       color = "Data") +
  theme_minimal() +
  theme(legend.position = c(0.4, 0.85),  # Ajustar la posición de la leyenda (x, y) dentro del área del gráfico
        #legend.background = element_rect(fill = alpha('white', 0.5)),
        legend.background = element_blank(),  # Eliminar el fondo de la leyenda
        legend.key = element_blank()) 
  #theme_minimal()

GB = ggplot(qqplot_data, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlim(100, 900) +  # Limitar el eje x
  ylim(100, 900) +  # Limitar el eje y
  labs(title = "QQ Plot",
       x = "Observed Quantiles",
       y = "Predicted Quantiles") +
  theme_minimal()
library(latex2exp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)
library(patchwork)

layout <- "
AABB
"
#p1 + p2 + p3 + p4 + 
#  plot_layout(design = layout)

GA + GC +  plot_layout(design = layout) +plot_annotation(tag_levels = list(c('A','B','C')))

##


sequencia = c(1,2,3,4,5,6,7,0,1,2,3,4,0,1,2,3,4,5,6,7)
#sequencia = c(1,2,3,4,5,6,7)
y1=sequencia * 1.931 # escalado al beta del modelo ajustado
Q2=0.25
y2=(1-(1-Q2)^(sequencia-1))*14.45 # escalado al beta del modelo ajustado
Q3=0.25;
for (n in 1:(length(sequencia)-1)){Q3[n+1] = Q3[n]+0.26*((1-(sequencia[n]>0))-Q3[n])} # alpha del modelo ajustado
y3=(1-(1-Q3)^(sequencia-1) )*42.093

Q4=0.25;
for (n in 1:(length(sequencia)-1)){Q4[n+1] = Q4[n]+0.01*((1-(sequencia[n]>0))-Q4[n])} # alpha del modelo ajustado
y4=(-Q4*250)+65 # escalado al beta del modelo ajustado



# Crear un data frame con los datos
data <- data.frame(
  x = rep(1:(length(sequencia)-2), 4),
  y = c(y1[which(sequencia>0)], y2[which(sequencia>0)], y3[which(sequencia>0)],y4[which(sequencia>0)]),
  #y = c(y1,y2,y3,y4),
  group = rep(c("M1", "M2", "M3(seq)","M3(q)"), each = length(sequencia)-2)
)

# Crear el gráfico con ggplot2
ggplot(data, aes(x = x, y = y, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Model-Based Expectation of Conflict", x = "Sequences of GO", y = "Expectation [A.U]") +
  theme_minimal() +
  geom_vline(xintercept = c(7.5, 11.5, 18.5), linetype = "dashed", color = "black") +
  annotate("text", x = c(18.5), y = 3, label = "Nogo", angle = 90, vjust = -0.5, color = "black") +
  scale_x_continuous(breaks = 1:18, labels = c(1:7,1:4,1:7)) +
  #scale_color_manual(values = c("blue", "red", "green")) +
  theme(legend.title = element_blank())


data$mn_new = as.factor(((data$seq<3)*-1 + (data$seq>5)*1 ) +2)
data$mn_new = as.factor(data$seq)
rt_final = aggregate(rt ~ mn_new, FUN=mean, data=data[data$rt>100 & data$est==1 & data$seq>0 & data$perror==0,]) 


library(dplyr)

# Calcular las medias y errores estándar
rt_summary <- data %>%
  filter(rt > 100, est == 1, seq > 0, perror == 0) %>%
  group_by(mn_new) %>%
  summarise(mean_rt = mean(rt),
            se_rt = sd(rt) / sqrt(n()))

# Crear el gráfico de barras con barras de error
ggplot(rt_summary, aes(x = mn_new, y = mean_rt)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = mean_rt - se_rt, ymax = mean_rt + se_rt), width = 0.2) +
  labs(x = "Mn New", y = "Reaction Time (Mean ± SE)", title = "Mean Reaction Time with Error Bars") +
  coord_cartesian(ylim = c(310, 350)) +  # Recortar el eje y entre 310 y 400
  theme_minimal()
