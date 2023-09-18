library(Rsolnp)
library(survRM2)
library(nnet)
library(dplyr)
library(ggplot2)
library(twang)
library(pracma)
library(gridExtra)
library(tidyverse)

source("source_MRCT.R")



########### 
# Set-Ups
###########
seed=12321
S=1000
n=c(400,500,600)
lambda=c(0.5,0.3)
gamma=c(1,0.3)
tau=4
M=list(c(0.5,1),c(0.5,1,1/3,2))
a = c(0.3,0.5,-1,0.5,0.3,0.5,-1,-0.5,-0.6,0.3,-1,0.5)

# log-linear sampling setting
r1 <- matrix(c(-5,1,0.3,-5,0.8,0.35,-5,0.6,0.4),nrow = 3,ncol = 3)  # moderate overlap
r2 <- matrix(c(-5,3,0.8,-5,2.5,0.9,-5,2,1),nrow = 3,ncol = 3) # low overlap
# logistic sampling setting
r1_mis <- matrix(c(-3,1,-0.1,-3,0.8,-0.15,-3,0.6,-0.2),nrow = 3,ncol = 3) # moderate overlap
r2_mis <- matrix(c(-3,3,-0.7,-3,2.5,-0.8,-3,2,-0.9),nrow = 3,ncol = 3) # low overlap




########### 
# Standardized Mean Difference
###########
nl <- c(50000,50000,50000)
set.seed(12321)
df1 <- gen.dat1(nl,r1,a3,lambda,gamma)
df2 <- gen.dat1(nl,r2,a3,lambda,gamma)
df1_mis <- gen.dat2(nl,r1_mis,a3,lambda,gamma)
df2_mis <- gen.dat2(nl,r2_mis,a3,lambda,gamma)

#############SMD###########
SMD(df1$X1,df1$R,0.5,sqrt(1/12))
SMD(df2$X1,df2$R,0.5,sqrt(1/12))
SMD(df1_mis$X1,df1_mis$R,0.5,sqrt(1/12))
SMD(df2_mis$X1,df2_mis$R,0.5,sqrt(1/12))
SMD(df1$X2,df1$R,1,1)
SMD(df2$X2,df2$R,1,1)
SMD(df1_mis$X2,df1_mis$R,1,1)
SMD(df2_mis$X2,df2_mis$R,1,1)
#############################




########### 
# True Region-Specific RMST
###########
nR=3
true_rmst <- true.rmst.fix(tau,nR,a,lambda,gamma)


##############
# Simulation Results
##############
res_r1 <- sim.MRCT.Est(seed,S,n,a,r1,tau,lambda,gamma,M)
res_r2 <- sim.MRCT.Est(seed,S,n,a,r2,tau,lambda,gamma,M)
res_r1_mis <- sim.MRCT.Est(seed,S,n,a,r1_mis,tau,lambda,gamma,M)
res_r2_mis <- sim.MRCT.Est(seed,S,n,a,r2_mis,tau,lambda,gamma,M)




##############
# Figures
##############
plot_r1 <- plot.data(res_r1)
plot_r2 <- plot.data(res_r2)
plot_r1_mis <- plot.data(res_r1_mis)
plot_r2_mis <- plot.data(res_r2_mis)


plot_delta1 <- rbind(plot_r1$delta1, plot_r2$delta1,
                     plot_r1_mis$delta1, plot_r2_mis$delta1)
plot_delta2 <- rbind(plot_r1$delta2, plot_r2$delta2,
                     plot_r1_mis$delta2, plot_r2_mis$delta2)
plot_delta3 <- rbind(plot_r1$delta3, plot_r2$delta3,
                     plot_r1_mis$delta3, plot_r2_mis$delta3)
plot_delta1$Difference <- plot_delta2$Difference <- plot_delta3$Difference <- 
  rep(c("moderate","low","moderate","low"), each = 4000)
plot_delta1$Sampling <- plot_delta2$Sampling <- plot_delta3$Sampling <- 
  rep(c("log","logit"), each = 8000)


plot_mu1_1 <- rbind(plot_r1$mu1_1, plot_r2$mu1_1,
                     plot_r1_mis$mu1_1, plot_r2_mis$mu1_1)
plot_mu1_2 <- rbind(plot_r1$mu1_2, plot_r2$mu1_2,
                     plot_r1_mis$mu1_2, plot_r2_mis$mu1_2)
plot_mu1_3 <- rbind(plot_r1$mu1_3, plot_r2$mu1_3,
                     plot_r1_mis$mu1_3, plot_r2_mis$mu1_3)
plot_mu1_1$Difference <- plot_mu1_2$Difference <- plot_mu1_3$Difference <- 
  rep(c("moderate","low","moderate","low"), each = 4000)
plot_mu1_1$Sampling <- plot_mu1_2$Sampling <- plot_mu1_3$Sampling <- 
  rep(c("log","logit"), each = 8000)


plot_mu0_1 <- rbind(plot_r1$mu0_1, plot_r2$mu0_1,
                    plot_r1_mis$mu0_1, plot_r2_mis$mu0_1)
plot_mu0_2 <- rbind(plot_r1$mu0_2, plot_r2$mu0_2,
                    plot_r1_mis$mu0_2, plot_r2_mis$mu0_2)
plot_mu0_3 <- rbind(plot_r1$mu0_3, plot_r2$mu0_3,
                    plot_r1_mis$mu0_3, plot_r2_mis$mu0_3)
plot_mu0_1$Difference <- plot_mu0_2$Difference <- plot_mu0_3$Difference <- 
  rep(c("moderate","low","moderate","low"), each = 4000)
plot_mu0_1$Sampling <- plot_mu0_2$Sampling <- plot_mu0_3$Sampling <- 
  rep(c("log","logit"), each = 8000)





# make the box plot
bias.plot <- function(data,true_rmst,ylab,ylim,breaks){
  data$Difference <- factor(data$Difference, 
                          levels = c("moderate","low"),
                          labels = c("moderate overlap","low overlap"))
  data$Sampling <- factor(data$Sampling,
                          levels = c("log","logit"),
                          labels = c("log","logit"))
  my.color <- c("darkorange2","seagreen","darkslateblue","hotpink2")
  #my.color <- c("#E889BD","#67C2A3","#FC8A61","#8EA0C9")
  ggplot(data, aes(x = Method, y = RMST, fill = Method)) +
    geom_boxplot() +
    coord_cartesian(ylim = ylim) +
    scale_y_continuous(breaks = breaks) +
    geom_hline(yintercept = true_rmst, linetype = "dotted", color = "red4", size = 0.7) +
    theme_bw() +
    theme(text = element_text(size=20),
          axis.text = element_text(size=15),
          legend.position="bottom") +
    #scale_x_continuous("Scenario", labels = as.character(data$S), breaks = data$S) +
    labs(y = ylab, x = "Estimator") +
    scale_color_manual(values = my.color) +
    facet_grid(cols = vars(Difference), rows = vars(Sampling))
}



p_delta1 <- bias.plot(plot_delta1, true_rmst_a3$delta_1,"RMSTD in Region 1",c(0,3),seq(0,3,0.5))
p_delta2 <- bias.plot(plot_delta2, true_rmst_a3$delta_2,"RMSTD in Region 2",c(0,3),seq(0,3,0.5))
p_delta3 <- bias.plot(plot_delta3, true_rmst_a3$delta_3,"RMSTD in Region 3",c(0,3),seq(0,3,0.5))
p_mu1_1 <- bias.plot(plot_mu1_1, true_rmst_a3$mu1_1,"RMST1 in Region 1",c(1.5,4.5),seq(1.5,4.5,0.5))
p_mu1_2 <- bias.plot(plot_mu1_2, true_rmst_a3$mu1_2,"RMST1 in Region 2",c(1.5,4.5),seq(1.5,4.5,0.5))
p_mu1_3 <- bias.plot(plot_mu1_3, true_rmst_a3$mu1_3,"RMST1 in Region 3",c(1.5,4.5),seq(1.5,4.5,0.5))
p_mu0_1 <- bias.plot(plot_mu0_1, true_rmst_a3$mu0_1,"RMST0 in Region 1",c(0,3),seq(0,3,0.5))
p_mu0_2 <- bias.plot(plot_mu0_2, true_rmst_a3$mu0_2,"RMST0 in Region 2",c(0,3),seq(0,3,0.5))
p_mu0_3 <- bias.plot(plot_mu0_3, true_rmst_a3$mu0_3,"RMST0 in Region 3",c(0,3),seq(0,3,0.5))




get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


pdf("Figure2.pdf", height = 15, width = 8, onefile = T)
grid.arrange(p_delta1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta3+theme(legend.position="none"),
             #get_legend(p_delta1),
             ncol = 1, heights = c(5,5,5),widths=8)
dev.off()



pdf("Web Figure 1 (left panel).pdf", height = 15, width = 8, onefile = T)
grid.arrange(p_mu1_1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_mu1_2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_mu1_3+theme(legend.position="none"),
             #get_legend(p_mu1_1),
             ncol = 1, heights = c(5,5,5),widths=8)

dev.off()



pdf("Web Figure 1 (right panel).pdf", height = 15, width = 8, onefile = T)
grid.arrange(p_mu0_1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_mu0_2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_mu0_3+theme(legend.position="none"),
             #get_legend(p_mu0_1),
             ncol = 1, heights = c(5,5,5),widths=8)

dev.off()




