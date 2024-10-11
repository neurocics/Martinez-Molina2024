# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module
library(runjags)
library(rjags)
load.module("wiener")


setwd("~/Documents/GitHub/Martinez-Molina2024/TMS")
source("../HDIofMCMC.r") 

# Draw random samples with JAGS
#data <- read.table("COR.txt", header=TRUE,sep="\t")
#data$Sub
#write.csv(data[ , c(1:5, 7,8,9,10,11,13)], file="Behavioral_TMS.csv")
data <- read.table("Behavioral_TMS.csv", header=TRUE,sep=",")

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




# correlative id for sessions
sort(unique(data$Dsu))
(Nsubj=length(unique(data$Dsu)))

# correlative id for subjects 
sort(unique(data$Dsu2))
(Nsubj=length(unique(data$Dsu2)))


# first trial for each subject 
dime =dim(data)
data$T1 = 1-c(0, data$Dsu[1:(dime[1]-1)]==data$Dsu[2:(dime[1])])

# n mini bloque
nmb=1
for (i in 2:dime[1]){
  if (data$nnb[i]>data$nnb[i-1]){
    nmb[i] = nmb[i-1]
  } else if (data$Dsu[i] != data$Dsu[i-1]  | nmb[i-1]==5 ){
    nmb[i] =1
  }   else { nmb[i]=nmb[i-1]+1 }
}
data$nmb=nmb


#

data$n2 = data$seq<3 & data$rt>40 & data$rt<1500
data$m2 = data$seq>=3 & data$rt>40 &data$rt<1500
data$faster = data$rt <  median(data$rt [data$rt>40 &data$rt<1500])


sum(data$n2) 
sum(data$n2) / (sum(data$n2) +sum(data$m2))

Dm = aggregate(rt~ m2+Dsu2+Dtms+Dtms51+DtmsA, FUN=median,data=data[(data$n2==1)|(data$m2==1) & data$rt>40,] )

wilcox.test(rt~m2,data=Dm[Dm$Dtms==0 ,],paired=T)
t.test(rt~m2,data=Dm[Dm$Dtms==0 ,],paired=T)
aggregate(rt~m2,data=Dm, FUN=median,)
wilcox_effsize(Dm[Dm$Dtms==0 ,], rt~m2,paired=T)

wilcox.test(rt~m2,data=Dm,paired=T)
t.test(rt~m2,data=Dm,paired=T)
aggregate(rt~m2,data=Dm, FUN=median,)
wilcox_effsize(Dm, rt~m2,paired=T)




E= Dm$rt[Dm$Dtms==0 &  Dm$m2==0] 
L= Dm$rt[Dm$Dtms==0 & Dm$m2==1] 

A = Dm$rt[Dm$Dtms==1 & Dm$Dtms51==0 & Dm$DtmsA==1 & Dm$m2==1] - Dm$rt[Dm$Dtms==1 & Dm$Dtms51==0 & Dm$DtmsA==1 & Dm$m2==0]
B = Dm$rt[Dm$Dtms==1 & Dm$Dtms51==1 & Dm$DtmsA==1 & Dm$m2==1] - Dm$rt[Dm$Dtms==1 & Dm$Dtms51==1 & Dm$DtmsA==1 & Dm$m2==0]

#Bp = Dm$rt[Dm$Dtms==0 & Dm$Dtms51==0 & Dm$DtmsA==1 & Dm$m2==1] - Dm$rt[Dm$Dtms==0 & Dm$Dtms51==0 & Dm$DtmsA==1 & Dm$m2==0]

#wilcox.test(B)
wilcox.test(A,B,paired=T)
t.test(A,B,paired=T)
wilcox_effsize(data.frame(D=c(A,B), G=c(rep(1,22), rep(2,22))), D~G,paired=T )
mean(A-B)

C = Dm$rt[Dm$Dtms==1 & Dm$Dtms51==0 & Dm$DtmsA==0 & Dm$m2==1] - Dm$rt[Dm$Dtms==1 & Dm$Dtms51==0 & Dm$DtmsA==0 & Dm$m2==0]
D = Dm$rt[Dm$Dtms==1 & Dm$Dtms51==1 & Dm$DtmsA==0 & Dm$m2==1] - Dm$rt[Dm$Dtms==1 & Dm$Dtms51==1 & Dm$DtmsA==0 & Dm$m2==0]
wilcox.test(C,D,paired=T)
t.test(C,D,paired=T)
wilcox_effsize(data.frame(D=c(C,D), G=c(rep(1,22), rep(2,22))), D~G,paired=T )

wilcox.test((A-B)-(C-D))
t.test((A-B)-(C-D))
wilcox_effsize(data.frame(D=c(A-B,C-D), G=c(rep(1,22), rep(2,22))), D~G,paired=T )




pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST=""}}
    }  
}



w_a0 = wilcox.test(E,L,paired=T)

forPlot = data.frame(rt=c(E,L) , late=c(rep(FALSE,44), rep(TRUE,44)) , subN = c(rep(1:44,2)) )


MP_A <- ggplot(forPlot , aes(x=late, y=rt, col=as.factor(late))) + 
  geom_boxplot(outlier.shape = NA,size = 0.5, width=0.4) + theme_classic()+
  geom_line(aes(x= ((late==1) - (late==1)*0.22  + ((late==0)*0.22) )+ 1 , y=rt, group = subN),color="grey",linetype=2,size=0.25)+
  #geom_line(aes(x= (p_conf_mayorQ==1), y=rt, group = subN),color="grey",linetype=2,size = 0.25,)+
  #geom_hline(yintercept=0.5,color="grey")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ylab("RT") +
  xlab("Trials in a sequence") +
  scale_x_discrete(labels=c("TRUE" = "late", "FALSE" = "early"))+ 
  geom_jitter( position=position_jitter(0.15),size = 2,alpha=0.7) +
  theme(legend.position = "none") +
  annotate(geom="text", x=1.5, y=max(forPlot$rt)+10, label=pv2str(w_a0$p.value),
          color="black") +
  geom_segment(aes(x = 1.3, y = max(rt), xend = 1.7, yend = max(rt)), color = "black") 
  #annotate(geom="text", x=1.5, y=max(rt)+20, label=pv2str(w_a0$p.value),
  #         color="black")

w_a0 = wilcox.test(A-B)
w_a01 = wilcox.test(A-B-C+D)

forPlot = data.frame(rt=c(B-A, D-C) , cite=c(rep(TRUE,22), rep(FALSE,22)) , subN = c(rep(1:22,2)) )

MP_B <- ggplot(forPlot , aes(x=cite, y=rt, col=as.numeric(cite*1))) + 
  geom_boxplot(outlier.shape = NA,size = 0.5, width=0.4) + theme_classic()+
  geom_line(aes(x= ((cite==1) - (cite==1)*0.22  + ((cite==0)*0.22) )+ 1 , y=rt, group = subN),color="grey",linetype=2,size=0.25)+
  #geom_line(aes(x= (p_conf_mayorQ==1), y=rt, group = subN),color="grey",linetype=2,size = 0.25,)+
  geom_hline(yintercept=0,color="grey")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  ylab("RT diff (late - early)") +
  xlab("diff (Theta - no Theta)") +
  scale_x_discrete(labels=c("TRUE" = "TMS SFG", "FALSE" = "TMS IFG"))+ 
  geom_jitter( position=position_jitter(0.15),size = 2,alpha=0.7) +
  theme(legend.position = "none") +
  annotate(geom="text", x=1.5, y=max(forPlot$rt)+5, label=pv2str(w_a01$p.value),
           color="black") +
  geom_segment(aes(x = 1.3, y = max(forPlot$rt)+2.5, xend = 1.7, yend = max(forPlot$rt)+2.5), color = "black") +
  annotate(geom="text", x=2, y=max(forPlot$rt)+2.5, label=pv2str(w_a0$p.value),
         color="black")



layout <- "
AB
"
#p1 + p2 + p3 + p4 + 
#  plot_layout(design = layout)

MP_A + MP_B +plot_layout(design = layout)





# only correct, go trial for RTs model 
Trt = which(data$rt<1500 & data$rt>40 & data$est==10 & !data$T1)
NGO = as.numeric(data$est==21)
Tng = intersect(round(Trt+1),round(which(data$est==21)))
data$error = as.numeric( (data$est==21 & data$rt>-99) | (data$est==10 & data$rt<40))

Ntotal = length(data$rt)

# data for JASG
#dat <- dump.format(list(Drt=data$rt, nT=Ntotal, Trt=Trt , Nsubj=Nsubj, idSub = data$Dsu,T1=data$T1, NGO=NGO,
#·                        Dseq=data$seq, Uno=data$Uno, pError=data$perror, Dtms=data$Dtms, Dtms51=data$Dtms51, DtmsA=data$DtmsA))


# Visualize the data
hist(dataG$rt, breaks=50)

# prepare the data for JAGS


# Initialize chains
inits1 <- dump.format(list(alph=20,ta=8,#bet=0.5,delt=0,mu.Q=0.2, 
                           mubeta0=9,mubeta1=1,mubeta2=0,mubeta3=0,mubeta4=-1,mubeta5=0,mubeta6=0,mubeta7=1,mubeta8=0,mubeta9=0,mubeta10=-1,mubeta11=0,mu.LR=0.2,
                           .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list(alph=150,ta=5,#bet=0.5,delt=0,mu.Q=0.3,  
                           mubeta0=9,mubeta1=0,mubeta2=1,mubeta3=-1,mubeta4=0,mubeta5=1,mubeta6=-1,mubeta7=0,mubeta8=1,mubeta9=-1,mubeta10=0,mubeta11=1,mu.LR=0.4,
                           .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list(alph=180,ta=6,#bet=0.5,delt=0, mu.Q=0.4, 
                           mubeta0=9,mubeta1=-1,mubeta2=-1,mubeta3=1,mubeta4=1,mubeta5=-1,mubeta6=1,mubeta7=-1,mubeta8=-1,mubeta9=1,mubeta10=-1,mubeta11=-1,mu.LR=0.6,
                           .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))

# Tell JAGS which latent variables to monitor
monitor = c('mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8',
            'mubeta9', 'mubeta10','mubeta11',"muWi0","ta","muWi1","delt","deviance","mu.LR")#
# data for JASG



dataall=data
data = data[data$nmb>=5,]
Trt = which(data$rt<1500 & data$rt>40 & data$est==10 & !data$T1)

NGO = as.numeric(data$est==21)
Ntotal = length(data$rt)

dat <- dump.format(list(Drt=data$rt, nT=Ntotal, Trt=Trt , Nsubj=Nsubj, idSub = data$Dsu2,
                        T1=data$T1, NGO=NGO,
                        Dseq=data$seq, Uno=data$Uno, pError=data$perror, Dtms=data$Dtms, 
                        Dtms51=data$Dtms51, DtmsA=data$DtmsA# , LR=0.37
                        ))

# partial data to find Learnig  Rate 
# Run the function that fits the models using JAGS
results_M3_DDM <- run.jags(model="model_M3_DDM_rl.txt",
                        monitor=monitor, data=dat, n.chains=3, 
                        inits=c(inits1, inits2, inits3), plots = TRUE,
                        burnin=2000, sample=1000, thin=5,
                        modules=c("wiener"), method=c("parallel"))

M3LR_ps =summary(results_M3_DDM)
for (s in 1:37){
  dataall$LRmedian[dataall$Dsu==s] = M3LR_ps[13+s,2]
  dataall$LRmean[dataall$Dsu==s] = M3LR_ps[13+s,4]
  dataall$LRmode[dataall$Dsu==s] = M3LR_ps[13+s,6]
}

# LR per subjet 
Trt = which(dataall$rt<1500 & dataall$rt>40 & dataall$est==10 & !dataall$T1)
Tng = intersect(round(Trt+1),round(which(dataall$est==21)))
NGO = as.numeric(dataall$est==21)
Ntotal = length(dataall$rt)

# AL darta to find all other parameters 
date <- dump.format(list(Drt=dataall$rt, nT=Ntotal, Trt=Trt , Nsubj=Nsubj, idSub = dataall$Dsu2,
                         T1=dataall$T1, NGO=NGO,Tng=Tng,error=dataall$error,
                         Dseq=dataall$seq, Uno=dataall$Uno, pError=dataall$perror, Dtms=dataall$Dtms, 
                         Dtms51=dataall$Dtms51, DtmsA=dataall$DtmsA, LR=as.numeric(dataall$LRmode)# , LR=0.37
))


# Run the function that fits the models using JAGS
results_M3_DDMrlps <- run.jags(model="model_M3_DDM_2.txt",
                           monitor=monitor, data=date, n.chains=3, 
                           inits=c(inits1, inits2, inits3), plots = TRUE,
                           burnin=2000, sample=1000, thin=5,
                           modules=c("wiener"), method=c("parallel"))


summary(results_M3_DDMrlps)

results_M3_DDMrlpsex = extend.jags(results_M3_DDMrlps,drop.monitor = c("LR"), add.monitor = monitor, burnin=0, sample=1000,adapt = 0)
summary(results_M3_DDMrlpsex)



chains = rbind(results_M3_DDMrlpsex$mcmc[[1]], results_M3_DDMrlpsex$mcmc[[2]], results_M3_DDMrlpsex$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC # .29   -->  462531.4   ind --> 462431.3
(pval = min(mean(chains[,"mubeta1"]<0)*2, mean(chains[,"mubeta1"]>0)*2))
(pval = min(mean(chains[,"mubeta5"]<0)*2, mean(chains[,"mubeta5"]>0)*2))
(pval = min(mean(chains[,"mubeta6"]<0)*2, mean(chains[,"mubeta6"]>0)*2))
(pval = min(mean(chains[,"mubeta7"]<0)*2, mean(chains[,"mubeta7"]>0)*2))#.29 [1] 0.02666667
                                                                        # ind    0.002666667
           #0.0006666667
(pval = min(mean(chains[,"mubeta4"]<0)*2, mean(chains[,"mubeta4"]>0)*2))


(pval = min(mean(chains[,"mubeta8"]<0)*2, mean(chains[,"mubeta8"]>0)*2))

#save.image("TMS_04092022.RData")


results_M3_DDM_cohe= extend.jags(results_M3_DDMrlps,drop.monitor = results_M3_DDMrlps$monitor, add.monitor = c("mubeta7","sigmabeta7"), burnin=0, sample=200,adapt = 0)
summary(results_M3_DDM_cohe)
chains = rbind(results_M3_DDM_cohe$mcmc[[1]], results_M3_DDM_cohe$mcmc[[2]], results_M3_DDM_cohe$mcmc[[3]])

D = (chains[,"mubeta7"]/chains[,"sigmabeta7"])/2
mean(D)
# -------------------------------------
# 
# -------------------------------------


# Run the function that fits the models using JAGS
results_M3_DDMlog <- run.jags(model="model_M3_DDM_log.txt",
                               monitor=monitor, data=date, n.chains=3, 
                               inits=c(inits1, inits2, inits3), plots = TRUE,
                               burnin=2000, sample=1000, thin=5,
                               modules=c("wiener"), method=c("parallel"))


summary(results_M3_DDMlog)

data$Trt=0
data$Trt[Trt]=1
data$Tng=0
data$Tng[Tng]=1

chains = rbind(results_M3_DDMlog$mcmc[[1]], results_M3_DDMlog$mcmc[[2]], results_M3_DDMlog$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC # .29   -->  462531.4   ind --> 462431.3  --> 471518.3 log + DDM

(pval = min(mean(chains[,"mubeta1"]<0)*2, mean(chains[,"mubeta1"]>0)*2))
(pval = min(mean(chains[,"mubeta5"]<0)*2, mean(chains[,"mubeta5"]>0)*2))
(pval = min(mean(chains[,"mubeta6"]<0)*2, mean(chains[,"mubeta6"]>0)*2))
(pval = min(mean(chains[,"mubeta7"]<0)*2, mean(chains[,"mubeta7"]>0)*2))#.29 [1] 0.02666667
# ind    0.002666667
(pval = min(mean(chains[,"mubeta4"]<0)*2, mean(chains[,"mubeta4"]>0)*2))

(pval = min(mean(chains[,"muWi1"]<0)*2, mean(chains[,"muWi1"]>0)*2))

(pval = min(mean((chains[,"muWi1"]*chains[,"mubeta7"])<0)*2, mean((chains[,"muWi1"]*chains[,"mubeta7"])>0)*2))


#

results_M3_log <- run.jags(model="model_M3_DDM_log.txt",
                              monitor=monitor, data=date, n.chains=3, 
                              inits=c(inits1, inits2, inits3), plots = TRUE,
                              burnin=2000, sample=1000, thin=5,
                              modules=c("wiener"), method=c("parallel"))


summary(results_M3_log)
chains = rbind(results_M3_log$mcmc[[1]], results_M3_log$mcmc[[2]], results_M3_log$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC # .29   -->  462531.4   ind --> 462431.3  --> 471518.3 log + DDM

(pval = min(mean(chains[,"mubeta1"]<0)*2, mean(chains[,"mubeta1"]>0)*2))
(pval = min(mean(chains[,"mubeta5"]<0)*2, mean(chains[,"mubeta5"]>0)*2))
(pval = min(mean(chains[,"mubeta6"]<0)*2, mean(chains[,"mubeta6"]>0)*2))
(pval = min(mean(chains[,"mubeta7"]<0)*2, mean(chains[,"mubeta7"]>0)*2))#.29 [1] 0.02666667
# ind    0.002666667
(pval = min(mean(chains[,"mubeta4"]<0)*2, mean(chains[,"mubeta4"]>0)*2))

(pval = min(mean(chains[,"muWi1"]<0)*2, mean(chains[,"muWi1"]>0)*2))

(pval = min(mean((chains[,"muWi1"]*chains[,"mubeta7"])<0)*2, mean((chains[,"muWi1"]*chains[,"mubeta7"])>0)*2))








cadena =  c(chains[,"mubeta1"],chains[,"mubeta4"],chains[,"mubeta5"],chains[,"mubeta6"],chains[,"mubeta7"])
plotData= data.frame(cadena)
plotData$beta = c(rep("1Exp",3000),rep("4TMS",3000),rep("5:TMS:Th",3000),rep("6TMS:A",3000),rep("7TMS:Th:A",3000))


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

p7=mean(chains[,"mubeta7"]<0)*2
p1=mean(chains[,"mubeta1"]<0)*2
mm=median(chains[,"mubeta1"])


PD_plot = ggplot(plotData, aes(y=cadena, x = beta, color=beta, fill = beta)) + 
  
  theme_minimal()+ 
  #geom_hline(yintercept=0,col="red") +
  # rep("5:TMS:Th",3000),rep("6TMS:A",3000),rep("7TMS:Th:A",3000))
  scale_x_discrete(labels=c("1Exp" = "Exp", 
                            "4TMS" = "TMS",
                            "5:TMS:Th" = "TMS:Th",
                            "6TMS:A" = "TMS:SFG",
                            "7TMS:Th:A" = "TMS:Th:SFG"))+ 
  annotate('rect', xmin=1.5, xmax=5.6, ymin=-12, ymax=22, alpha=.2, fill='#999999') +
  scale_fill_manual(values=c("#77ac30ff","#999999", "#999999", "#999999", "#77ac30ff"))+
  scale_color_manual(values=c("#77ac30ff","#999999", "#999999", "#999999", "#77ac30ff"))+
  geom_violin( )+
  geom_segment(aes(x = 1.6, y = 0, xend =5.6 , yend = 0), color = "red") +
  geom_segment(aes(x = 1.6, y = 0, xend =1.4 , yend = mm), color = "red",linetype="dotted") +
  geom_segment(aes(x = 1.4, y =mm, xend =0.5 , yend = mm), color = "red") +
  
  ggtitle("Parameter Distribution")+  #  
  ylab("Posterior distribution (a.u.)")+
  xlab("Parameters")+
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  # geom_boxplot(width=0.1) +
  coord_flip()+
  annotate(geom="text", y=15, x=5.3, label=pv2str(p7),
           color="black") +
  annotate(geom="text", y=19, x=1.3, label=pv2str(p1),
           color="black") +
  annotate(geom="text", y=-5.5, x=5.4, label="Exp interaction",
           color="black") +
  theme(legend.position = "none") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")

#dev.off()


layout <- "
AACC
BBCC
"
#p1 + p2 + p3 + p4 + 
#  plot_layout(design = layout)

MP_A + MP_B + PD_plot +plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c('B','C','D')))







