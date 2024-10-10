# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
### plot 
library(latex2exp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggcharts)



setwd("~/Documents/GitHub/Martinez-Molina2024/TMS")
source("../HDIofMCMC.r") 

data_summary <- function(x) {
  m <- median(x)
  hdi = HDIofMCMC(x , credMass=0.95)
  ymin <- hdi[1]
  ymax <-  hdi[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}




#save(E,L,A,B,C,D,chains, file = "data_Figure_3.Rdata")
load("data_Figure_3.Rdata")

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


cadena =  c(chains[,"mubeta1"],chains[,"mubeta4"],chains[,"mubeta5"],chains[,"mubeta6"],chains[,"mubeta7"])
plotData= data.frame(cadena)
plotData$beta = c(rep("1Exp",3000),rep("4TMS",3000),rep("5:TMS:Th",3000),rep("6TMS:A",3000),rep("7TMS:Th:A",3000))


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

