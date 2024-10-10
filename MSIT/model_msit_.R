# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module
library(runjags)
library(rjags)
library(loo)
library(lme4)
library(ggplot2)
library(patchwork)
library(rstatix)
load.module("wiener")

#setwd("/Users/pablobilleke/Library/CloudStorage/OneDrive-udd.cl/Working papers/Martinez_theta/OLDs/DATA/MSIT")
setwd("~/Documents/GitHub/Martinez-Molina2024/MSIT")
source("../HDIofMCMC.r") 



DATA <- read.table("Behavioral_MSIT.csv", header=TRUE,sep=",")
names(DATA)
#  "rt"    "est"   "resp"  "laten" "good"  "acu"   "seq"  "subN" 

hist(DATA$rt)

#


Trt =which(DATA$rt>100)

Data_f = DATA[Trt,]

early = Data_f$seq<3
late = Data_f$seq>2 & Data_f$seq<6

Data_f$p_conf_mayorQ = early * (Data_f$est==11) + late * (Data_f$est==12)
Data_f$p_conf_menorQ = early * (Data_f$est==12) + late * (Data_f$est==11)

mean(Data_f$p_conf_mayorQ)
mean(Data_f$p_conf_menorQ)



forTest = aggregate(acu ~ p_conf_mayorQ*p_conf_menorQ*subN, FUN=mean, data=Data_f[(Data_f$p_conf_mayorQ | Data_f$p_conf_menorQ) & Data_f$est==12 ,  ] )
aggregate(acu ~ p_conf_mayorQ, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ,  ], FUN=mean)
wilcox.test(acu ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])
wilcox_effsize(acu ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])
t.test(acu ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])

forTest = aggregate(rt ~ p_conf_mayorQ*p_conf_menorQ*subN, FUN=median, data=Data_f[Data_f$p_conf_mayorQ | Data_f$p_conf_menorQ,  ] )
aggregate(rt ~ p_conf_mayorQ, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ,  ], FUN=mean)
wilcox.test(rt ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])
wilcox_effsize(rt ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])
t.test(rt ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])

##

pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST=""}}
    }  
}



w_a0 = wilcox.test(rt ~ p_conf_mayorQ, paired=T, data=forTest[forTest$p_conf_mayorQ | forTest$p_conf_menorQ ,  ])




MP_p <- ggplot(forTest, aes(x=(p_conf_mayorQ==1), y=rt, col=as.factor(p_conf_mayorQ))) + 
  geom_boxplot(outlier.shape = NA,size = 0.5, width=0.4) + theme_classic()+
  geom_line(aes(x= ((p_conf_mayorQ==1) - (p_conf_mayorQ==1)*0.22  + ((p_conf_mayorQ==0)*0.22) )+ 1 , y=rt, group = subN),color="grey",linetype=2,size=0.25)+
  #geom_line(aes(x= (p_conf_mayorQ==1), y=rt, group = subN),color="grey",linetype=2,size = 0.25,)+
  
  #geom_hline(yintercept=0.5,color="grey")+ 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ylab("RT (ms)") +
  xlab("Trials divided by p(cnf)") +
  scale_x_discrete(labels=c("TRUE" = "greater", "FALSE" = "lower"))+ 
  geom_jitter( position=position_jitter(0.15),size = 2,alpha=0.7) +
  
  theme(legend.position = "none") +
  #annotate(geom="text", x=1, y=0.05, label=pv2str(w_a0$p.value),
  #        color="black") +
  geom_segment(aes(x = 1.3, y = max(rt)+10, xend = 1.7, yend = max(rt)+10), color = "black") +
  annotate(geom="text", x=1.5, y=max(forTest$rt)+20, label=pv2str(w_a0$p.value),
           color="black") 
MP_p

d <- data.frame(Qmenor = forTest$rt[forTest$p_conf_menorQ==1], Qmayor = forTest$rt[forTest$p_conf_mayorQ==1])
#ggpaired(d, cond1 = "Qmenor", cond2 = "Qmayor",
#         fill = "condition", palette = "jco")





##


Wrt = DATA$rt * (DATA$acu + (DATA$acu-1))

(Nsubj=length(unique(DATA$subN)))

hist(Wrt[Trt],breaks = 100)

dime =dim(DATA)


#M1 = lmer(rt ~ I(est==12)*seq +  pError + pCNF+   (1 + I(est==12)*seq+  pError + pCNF| subN), data=DATA[DATA$rt>100,])
#summary(M1)



#
DATA$CNF = as.numeric(DATA$est==12)
Ntotal = length(DATA$rt)
DATA$pError = c(0,DATA$acu[1:dime[1]-1])
DATA$pCNF = c(0,DATA$CNF[1:dime[1]-1])
DATA$T1 = 1-c(0, DATA$subN[1:(dime[1]-1)]==DATA$subN[2:(dime[1])])

DATA$ISI = c(0,diff(DATA$laten))
# data for JAS
dat <- dump.format(list(Wrt=Wrt/1000, nT=Ntotal, Nsubj=Nsubj, idSub = DATA$subN, 
                        Dseq=DATA$seq, pError=DATA$pError,CNF=DATA$CNF,pCNF=DATA$pCNF,
                        T1=DATA$T1,Trt=Trt,ISI=DATA$ISI
))



# Initialize chains
inits3 <- dump.format(list(alph=.020,ta=.008,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=9,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-.1,mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,mubeta10=-1,mubeta11=0,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=.050,ta=.005,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=9,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=.1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits1 <- dump.format(list(alph=.030,ta=.006,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=9,mubeta1=.1,mubeta2=.1,mubeta3=-.1,mubeta4=-.1,mubeta5=-.1,mubeta6=0,mubeta7=0,mubeta8=1,mubeta9=1,mubeta10=0,mubeta11=1,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=646 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3','mubeta4', 'mu.LR',"alph","ta",
            "bet","delt","mu.Q","deviance",
            "th","muthr")



# Run the function that fits the models using JAGS
resultsM1 <- run.jags(model="DDM_ln.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=1000, sample=1000, thin=5, modules=c("wiener"), 
                    method=c("parallel"))

resultsM1
chains = rbind(resultsM1$mcmc[[1]], resultsM1$mcmc[[2]], resultsM1$mcmc[[3]])
DIC.M1 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M1 #[1]  3431.99





# Run the function that fits the models using JAGS
resultsM2 <- run.jags(model="DDM_exp.txt",
                      monitor=monitor, data=dat, n.chains=3, 
                      inits=c(inits1, inits2, inits3), plots = TRUE,
                      burnin=2000, sample=1000, thin=5, modules=c("wiener"), 
                      method=c("parallel"))

resultsM2
chains = rbind(resultsM2$mcmc[[1]], resultsM2$mcmc[[2]], resultsM2$mcmc[[3]])
DIC.M2 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M2#[1] 3506.518


# Run the function that fits the models using JAGS
resultsM3expRL <- run.jags(model="DDM_LR.txt",
                    monitor=monitor, data=dat, n.chains=3, 
                    inits=c(inits1, inits2, inits3), plots = TRUE,
                    burnin=1000, sample=1000, thin=5, modules=c("wiener"), 
                    method=c("parallel"))

resultsM3expRL2 <- extend.jags(resultsM3expRL, sample=1000, thin=5)

# Run the function that fits the models using JAGS
resultsM4 <- run.jags(model="DDM_LR_M4.txt",
                           monitor=monitor, data=dat, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=1000, sample=1000, thin=5, modules=c("wiener"), 
                           method=c("parallel"))


chains = rbind(resultsM3expRL2$mcmc[[1]], resultsM3expRL2$mcmc[[2]], resultsM3expRL2$mcmc[[3]])

rlps = extend.jags(resultsM3expRL,drop.monitor = resultsM3expRL$monitor,add.monitor = c("LR"),burnin=0, sample=200,adapt = 0)


rlpss =summary(rlps)
for (s in 1:26){
  DATA$LRmedian[DATA$subN==s] = rlpss[s,2]
  DATA$LRmean[DATA$subN==s] = max(0.015,rlpss[s,4])
  DATA$LRmode[DATA$subN==s] = rlpss[s,6]
}

hist(DATA$LRmean)
# [1] 3510.939 (sin Q con Exp)  
#  3699 (exp sinQ sin pC)

#save.image("msit_13092022.RData")

# data for JAS
dat <- dump.format(list(Wrt=Wrt/1000, nT=Ntotal, Nsubj=Nsubj, idSub = DATA$subN, 
                        Dseq=DATA$seq, pError=DATA$pError,CNF=DATA$CNF,pCNF=DATA$pCNF,
                        T1=DATA$T1,Trt=Trt, LRi=DATA$LRmean,Uno=as.numeric(DATA$seq==1)
))



# Run the function that fits the models using JAGS
resultsM3fixLR <- run.jags(model="DDM_LR_lr.txt",
                           monitor=monitor, data=dat, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=20000, sample=1000, thin=20, modules=c("wiener"), 
                           method=c("parallel"))
resultsM3fixLR
chains = rbind(resultsM3fixLR$mcmc[[1]], resultsM3fixLR$mcmc[[2]], resultsM3fixLR$mcmc[[3]])
DIC.M3fixLR = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3fixLR # [1] 3458.836

#save.image("msit_04102022.RData")


##3  Figure 
chains = rbind(resultsM3fixLR$mcmc[[1]], resultsM3fixLR$mcmc[[2]], resultsM3fixLR$mcmc[[3]])
cadena =  c(chains[,"mubeta1"]/2,
            chains[,"mubeta2"]/1,
            chains[,"mubeta3"]/1
            #chains[,"nor"],
            #chains[,"ta"]/200,
            #chains[,"mu.LR"]
)
plotData= data.frame(cadena)
plotData$beta = c(rep("Exp",3000),
                  rep("CNF",3000),
                  rep("pEr",3000)
                  #rep("Driff",3000),
                  #rep("tau",3000),
                  #rep("LR",3000)
)


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

pb1=mean(chains[,"mubeta1"]<0)*2
pb2=mean(chains[,"mubeta2"]<0)*2
pb3=mean(chains[,"mubeta3"]<0)*2


PP = ggplot(plotData, aes(y=cadena, x = beta, color=beta, fill = beta)) + 
  geom_violin( )+
  #theme_minimal()+ 
  geom_hline(yintercept=0,col="red") +
  scale_fill_manual(values=c("#00bfc4ff","#77ac30ff","#808080ff"))+
  scale_color_manual(values=c("#00bfc4ff","#77ac30ff","#808080ff"))+
  #scale_fill_manual(values=c("#999899","#999999", "#56B4E9", "#E69F00", "#42C01A"))+
  #scale_color_manual(values=c("#999899","#999999", "#56B4E9", "#E69F00", "#42C01A"))+
  ggtitle("Reaction time slowing MSIT task")+  #  
  ylab("Posterior probability (a.u.)")+
  xlab("Model parameters")+
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  # geom_boxplot(width=0.1) +
  theme(legend.position = "none") +
  #scale_y_log10() +
  coord_flip()+
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=0.2, x=3.3, label=pv2str(pb3),
           color="black") +
  annotate(geom="text", y=1, x=2.3, label=pv2str(pb2),
           color="black") +
  annotate(geom="text", y=0.45, x=1.3, label=pv2str(pb1),
           color="black") 

PP





layout <- "
AABBB
"

MP_p + PP+plot_layout(design = layout) +plot_annotation(tag_levels = list(c('C','D')))


################  modelo controles 

# Run the function that fits the models using JAGS
resultsM3real <- run.jags(model="DDM_LR_exp.txt",
                           monitor=monitor, data=dat, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=5000, sample=1000, thin=5, modules=c("wiener"), 
                           method=c("parallel"))

chains = rbind(resultsM3real$mcmc[[1]], resultsM3real$mcmc[[2]], resultsM3real$mcmc[[3]])
DIC.M3real = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M3real # [1] 





# Run the function that fits the models using JAGS
resultsM5<- run.jags(model="../Exp1/model_DDM_LR_theta.txt",
                          monitor=monitor, data=dat, n.chains=3, 
                          inits=c(inits1, inits2, inits3), plots = TRUE,
                          burnin=5000, sample=1000, thin=5, modules=c("wiener"), 
                          method=c("parallel"))

chains = rbind(resultsM5$mcmc[[1]], resultsM5$mcmc[[2]], resultsM5$mcmc[[3]])
DIC.M5 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC.M5 # [1] 






