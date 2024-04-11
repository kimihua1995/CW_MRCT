library(Rsolnp)
library(survRM2)
library(nnet)
library(dplyr)
library(ggplot2)
library(twang)
library(pracma)
library(gridExtra)
library(tidyverse)
library(nleqslv)
library(geex)

source("source_estimator.R")
source("source_simulation.R")



########### 
# Set-Ups
###########
lambda <- c(0.5,0.3)
gamma <- c(1,0.3)
a = c(0.3,0.1,-1,0.5,0.3,0.7,-1,-0.5,-0.6,0.3,-0.7,0.3)
r1 <- matrix(c(-5,0.8,0.3,-5,0.7,0.27,-5,0.6,0.25),nrow = 3,ncol = 3)  # moderate
r2 <- matrix(c(-5,2.5,0.5,-5,2.3,0.55,-5,2,0.6),nrow = 3,ncol = 3) # large
r1_mis <- matrix(c(-3,0.6,-0.15,-3,0.5,-0.1,-3,0.4,-0.05),nrow = 3,ncol = 3) # moderate
r2_mis <- matrix(c(-2.3,3,-0.2,-2.3,2.5,-0.15,-2.3,2,-0.1),nrow = 3,ncol = 3) # large
n <- c(400,500,600)
tau <- 4
seed <- 12321
S <- 1000
g <- c(0.5,1,1/3,2)
true_rmst <- true.rmst.fix(tau,3,a,lambda,gamma)



res_r1 <- sim.MRCT.Est(seed,S,n,a,r1,tau,lambda,gamma,g,1)
res_r2 <- sim.MRCT.Est(seed,S,n,a,r2,tau,lambda,gamma,g,1)
res_r1_mis <- sim.MRCT.Est(seed,S,n,a,r1_mis,tau,lambda,gamma,g,2)
res_r2_mis <- sim.MRCT.Est(seed,S,n,a,r2_mis,tau,lambda,gamma,g,2)



########### 
# Standardized Mean Difference
###########
nl <- c(30000,30000,30000)
set.seed(12321)
df1 <- gen.dat1(nl,r1,a,lambda,gamma)
df2 <- gen.dat1(nl,r2,a,lambda,gamma)
df1_mis <- gen.dat2(nl,r1_mis,a,lambda,gamma)
df2_mis <- gen.dat2(nl,r2_mis,a,lambda,gamma)

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


### load saved results
load("res.RData")



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
# Check Results
##############
#res.check(res_r1$mu.IPW1, res_r1$sd.IPW1, true_rmst)
#res.check(res_r1$mu.IPW2, res_r1$sd.IPW2, true_rmst)
#res.check(res_r1$mu.EW, res_r1$sd.EW, true_rmst)


##############
# Figures
##############
plot_r1 <- plot.data(res_r1)
plot_r2 <- plot.data(res_r2)
plot_r1_mis <- plot.data(res_r1_mis)
plot_r2_mis <- plot.data(res_r2_mis)



plot_list_IPW1 <- plot_list_IPW2 <- plot_list_CW <- list()
for (i in 1:9){
  plot_list_IPW1[[i]] <- rbind(plot_r1$IPW1[[i]], plot_r1_mis$IPW1[[i]], 
                               plot_r2$IPW1[[i]], plot_r2_mis$IPW1[[i]])
  plot_list_IPW2[[i]] <- rbind(plot_r1$IPW2[[i]], plot_r1_mis$IPW2[[i]],
                               plot_r2$IPW2[[i]], plot_r2_mis$IPW2[[i]])
  plot_list_CW[[i]] <- rbind(plot_r1$CW[[i]], plot_r1_mis$CW[[i]],
                             plot_r2$CW[[i]], plot_r2_mis$CW[[i]])
  plot_list_IPW1[[i]]$Sampling <- plot_list_IPW2[[i]]$Sampling <- 
  plot_list_CW[[i]]$Sampling <- rep(c("log","logit","log","logit"), each = 700)
  plot_list_IPW1[[i]]$Difference <- plot_list_IPW2[[i]]$Difference <- 
  plot_list_CW[[i]]$Difference <-  rep(c("moderate","low"), each = 1400)
}




# make the plot
bias.plot <- function(data,true_rmst,ylab,ylim,breaks){
  data$Difference <- factor(data$Difference, 
                          levels = c("moderate","low"),
                          labels = c("Moderate SMD","Large SMD"))
  data$Sampling <- factor(data$Sampling,
                          levels = c("log","logit"),
                          labels = c("Log-Linear Sampling","Logistic Sampling"))
  #my.color <- c("darkorange2","seagreen","darkslateblue","hotpink2")
  #my.color <- c("#E889BD","#67C2A3","#FC8A61","#8EA0C9")
  data$Estimator <- data$Method
  #my.color <- c("gray73","darkorange3","darkorange1","darkgoldenrod3","darkgoldenrod1","brown2","orangered1",
  #              "dodgerblue2","slateblue2","cadetblue2","cornflowerblue","seagreen2","aquamarine2")
  my.color <- c("gray73","#8DD3C7","#FFED6F","#BC80BD","#FB8072","#80B1D3","#FDB462")
  ggplot(data, aes(x = Estimator, y = RMST, fill = Estimator)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun=mean, colour="black", geom="point", 
                 shape=18, size=2.5, show.legend=FALSE) +
    coord_cartesian(ylim = ylim) +
    scale_y_continuous(breaks = breaks) +
    geom_hline(yintercept = true_rmst, linetype = "dotted", color = "red2", size = 0.7) +
    theme_bw() +
    theme(text = element_text(size=15),
          axis.text = element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="right") +
    #scale_x_continuous("Scenario", labels = as.character(data$S), breaks = data$S) +
    labs(y = ylab, x = "Estimator") +
    scale_fill_manual(name = "Estimators", values = my.color) +
    facet_grid(cols = vars(Sampling), rows = vars(Difference))
}





pdf("RMSTD plots.pdf", height = 15, width = 10, onefile = T)
grid.arrange(p_delta1_IPW1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta1_IPW2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta1_CW+theme(legend.position="bottom"),
             #get_legend(p_delta1),
             ncol = 1, heights = c(5,5,5),widths=10)
grid.arrange(p_delta2_IPW1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta2_IPW2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta2_CW+theme(legend.position="bottom"),
             #get_legend(p_delta1),
             ncol = 1, heights = c(5,5,5),widths=10)
grid.arrange(p_delta3_IPW1+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta3_IPW2+theme(axis.title.x = element_blank(),legend.position="none"),
             p_delta3_CW+theme(legend.position="bottom"),
             #get_legend(p_delta1),
             ncol = 1, heights = c(5,5,5),widths=10)
dev.off()


p_delta1_IPW1 <- bias.plot(plot_list_IPW1[[1]], true_rmst$delta_1, 
                           expression("RMSTD in Region 1, IPSW, True"~rho[r](X)),
                           c(1.0,2.5),seq(1.0,2.5,0.5))
p_delta1_IPW2 <- bias.plot(plot_list_IPW2[[1]], true_rmst$delta_1,
                           expression("RMSTD in Region 1, IPSW, Est"~rho[r](X)),
                           c(1.0,2.5),seq(1.0,2.5,0.5))
p_delta1_CW <- bias.plot(plot_list_CW[[1]], true_rmst$delta_1,
                         "RMSTD in Region 1, CW",
                         c(1.0,2.5),seq(1.0,2.5,0.5))

p_delta2_IPW1 <- bias.plot(plot_list_IPW1[[2]], true_rmst$delta_2, 
                           expression("RMSTD in Region 2, IPSW, True"~rho[r](X)),
                           c(1.0,2.5),seq(1.0,2.5,0.5))
p_delta2_IPW2 <- bias.plot(plot_list_IPW2[[2]], true_rmst$delta_2,
                           expression("RMSTD in Region 2, IPSW, Est"~rho[r](X)),
                           c(1.0,2.5),seq(1.0,2.5,0.5))
p_delta2_CW <- bias.plot(plot_list_CW[[2]], true_rmst$delta_2,
                         "RMSTD in Region 2, CW",
                         c(1.0,2.5),seq(1.0,2.5,0.5))


p_delta3_IPW1 <- bias.plot(plot_list_IPW1[[3]], true_rmst$delta_3, 
                           expression("RMSTD in Region 3, IPSW, True"~rho[r](X)),
                           c(0.5,2.0),seq(0.5,2.0,0.5))
p_delta3_IPW2 <- bias.plot(plot_list_IPW2[[3]], true_rmst$delta_3,
                           expression("RMSTD in Region 3, IPSW, Est"~rho[r](X)),
                           c(0.5,2.0),seq(0.5,2.0,0.5))
p_delta3_CW <- bias.plot(plot_list_CW[[3]], true_rmst$delta_3,
                         "RMSTD in Region 3, CW",
                         c(0.5,2.0),seq(0.5,2.0,0.5))

