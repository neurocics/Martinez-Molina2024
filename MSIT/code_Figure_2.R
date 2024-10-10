# code Figure 2

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

# Load rjags library and the wiener module

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

setwd("~/Documents/GitHub/Martinez-Molina2024/MSIT/")
source("../HDIofMCMC.r") 



#save(chains,
#     forTest,
#     file="data_Figure_2.RData")

load("data_Figure_2.RData")



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





##3  Figure 
#chains = rbind(resultsM3fixLR$mcmc[[1]], resultsM3fixLR$mcmc[[2]], resultsM3fixLR$mcmc[[3]])
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

