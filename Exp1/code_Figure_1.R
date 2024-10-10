# code Figure 1

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module

library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)


setwd("~/Documents/GitHub/Martinez-Molina2024/Exp1/")
source("../HDIofMCMC.r") 
### plot 


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



#save(chains,loo.M0ddm,
#     loo.M1ddm,#
#     loo.M2ddm,#
#     loo.M3ddm,#
#     loo.M4ddm,#
#     DIC.M0ddm, 
#     DIC.M1ddm,
#     DIC.M2ddm,
#     DIC.M3ddm,
#     DIC.M4ddm,
#     Dm,
#     file="data_Figure_1.RData")

load("data_Figure_1.RData")

Rwil = wilcox.test(rt~m2,data=Dm,paired=T)
##3  Figure 

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







#chains = rbind(results_DDM_LR $mcmc[[1]], results_DDM_LR $mcmc[[2]], results_DDM_LR $mcmc[[3]])
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



#tikzDevice::tikz(file = "./taus_TMS.tex", width = 3, height = 3)

PD = ggplot(plotData, aes(y=cadena, x = beta, color=beta, fill = beta)) + 
  geom_violin( )+
  #theme_minimal()+ 
  geom_hline(yintercept=0,col="red") +
  scale_fill_manual(values=c("#77ac30ff","#808080ff","#808080ff"))+
  scale_color_manual(values=c("#77ac30ff","#808080ff","#808080ff"))+
  #ggtitle("Reaction time slowing")+  #  
  ylab("Posterior distribution (arb. units)")+
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



