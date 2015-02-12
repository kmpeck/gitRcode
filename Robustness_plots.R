#-----------
# Figure 1: 
# Contour plot of environmental vs. genetic robustness for random individuals
#-----------

setwd("~/GitHub/gene_network/data/Manuscript data/Figure 1")
library("gplots")
dat <- read.table("10000_random_individuals.txt", header=T)
dat$g_rob <- round(dat$genetic_robustness,1)
dat$e_rob <- round(dat$environmental_robustness,1)

dat.tab2 <- as.matrix(table(dat$e_rob, dat$g_rob))
colors = c("white","grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10","black","black")

pdf("Contour_Plot.pdf", width=7, height=7)
par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)
filled.contour(dat.tab2, col=colors, nlevels=9, xlab="environmental robustness",
               ylab="genetic robustness")
dev.off()

out1 <- cor.test(dat$environmental_robustness, dat$genetic_robustness, method="spearman")
out1


#-----------
# Figure 2
#-----------

# Environmental and genetic robustness before and after evolution (in a control and perturbed environment)

data_summary <- function(treatment_dir){
  for(i in 1:length(treatment_dir)){
  setwd(treatment_dir[i])
  
  #read in directory names for each experiment
  experiment_dir <- dir(getwd(), pattern="founder+", full.names=T)
  message("Number of experiments: ", length(experiment_dir))
  if (length(experiment_dir) == 0) {
    stop("No experiment directories found! Exiting.")
  }
  
  #create data_summary file (will be stored in treatment folders)
  fileName = "data_summary.txt"
  fileHeader = c("Replicate", "Generation", "Mean_Fitness", "Mean_Robustness", "Median_Robustness", "Mean_Env_Robustness", "Median_Env_Robustness",
                 "path_length", "path_length_ancestral", "mean_weight_diagonal", "mean_weight_diagonal_nonzero", "prop_positive_on_diagonal", 
                 #"prop_positive_off_diagonal", "nonzero_diagonal", "mean_weight_all",
                 "nonzero_diagonal", "mean_weight_all",
                 "standard_deviation_all", "matrix_asymmetry", "genetic_variation", "genetic_variation_euclidian",
                 "genetic_robustness_ancestral", "environmental_robustness_ancestral")
  nCol = length(fileHeader)
  write(fileHeader, fileName, ncolumns=nCol)
  
  for(exp in 1:length(experiment_dir)){
    setwd(experiment_dir[exp])
    output_files <- list.files(getwd(),pattern = "data+")
    if (length(output_files) == 0) {
      stop("No output files found! Exiting.")
    }
    generations <- NULL
    for (e in 1:length(output_files)) {
      setwd(experiment_dir[exp])
      mydata = read.table(output_files[e], header=T)
      message(output_files[e])
      generations[e] <- gsub("([a-z_]+)([0-9]+)([a-z\\.]+)", "\\2", output_files[e])
      meanFitness <- mean(mydata$fitness)
      meanRobustness <- mean(mydata$genetic_robustness)
      medianRobustness <- median(mydata$genetic_robustness)
      meanEnvRobustness <- mean(mydata$environmental_robustness)
      medianEnvRobustness <- median(mydata$environmental_robustness)
      meanPL <- mean(mydata$path_length)
      meanPLA <- mean(mydata$path_length_ancestral)
      meanWD <- mean(mydata$mean_weight_diagonal)
      meanWDNZ <- mean(mydata$mean_weight_diagonal_nonzero)
      meanWDP <- mean(mydata$number_positive_on_diagonal)/10
      #meanWODP <- mean(mydata$number_positive_off_diagonal)/90
      meanND <- mean(mydata$nonzero_diagonal)
      meanMWA <- mean(mydata$mean_weight_all)
      meanSDA <- mean(mydata$standard_deviation_all)
      meanMA <- mean(mydata$matrix_asymmetry)
      meanGV <- mean(mydata$genetic_variation)
      meanGVE <- mean(mydata$genetic_variation_euclidian)
      medianGRA <- 0
      medianERA <- 0
      if (i == 2){
        medianGRA <- median(mydata$genetic_robustness_ancestral)
        medianERA <- median(mydata$environmental_robustness_ancestral)
      }
      
      results = c(exp, as.numeric(generations[e]), meanFitness, meanRobustness, medianRobustness, meanEnvRobustness, medianEnvRobustness,
                  meanPL,meanPLA,meanWD,meanWDNZ,meanWDP,
                  #meanWODP,
                  meanND,meanMWA,meanSDA,meanMA,meanGV, meanGVE,
                  medianGRA, medianERA)
      setwd(treatment_dir[i])
      write(results, fileName, ncolumns=nCol, append=T)
      
    }
  }
}
}

setwd("~/GitHub/gene_network/data/Manuscript data/Figure 2 new/default6")
main_dir = getwd()
treatment_dir <- NULL #asex and sex
treatment_dir[1] <- paste(main_dir,"/control_pop",sep="")
treatment_dir[2] <- paste(main_dir,"/perturb_pop",sep="")

data_summary(treatment_dir)

setwd(main_dir)
pdf("Robustness_Boxplots.pdf", width=10, height=7)
par(mfrow=c(1,2))
par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)

#get data needed to set up empty plot
setwd(treatment_dir[1])
mydat <- read.table("data_summary.txt", header=T)
env_start <- subset(mydat$Median_Env_Robustness, mydat$Generation == 0)
env_end <- subset(mydat$Median_Env_Robustness, mydat$Generation == max(mydat$Generation))
gen_start <- subset(mydat$Median_Robustness, mydat$Generation == 0)
gen_end <- subset(mydat$Median_Robustness, mydat$Generation == max(mydat$Generation))
setwd(treatment_dir[2])
mydat2 <- read.table("data_summary.txt", header=T)
env_end_pert <- subset(mydat2$Median_Env_Robustness, mydat2$Generation == max(mydat2$Generation))
gen_end_pert <- subset(mydat2$genetic_robustness_ancestral, mydat2$Generation == max(mydat2$Generation))
setwd(main_dir)

boxplot(env_start, env_end, env_end_pert, names=c("start", "end - control", "end"), main="Environmental Robustness", ylim=c(0,1), col="grey70")
boxplot(gen_start, gen_end, gen_end_pert, names=c("start", "end - control", "end"), main="Genetic Robustness", ylim=c(0,1), col="grey70")
dev.off()

#Wilcoxon sign-rank test
#Does not assume a normal distribution
library(MASS)
filestats="Population_stats.txt"
write("Wilcoxon rank tests", filestats)
wilcox.test(env_start, env_end_pert, paired=T) -> env_pert; write(c(env_pert$data.name, env_pert$statistic, signif(env_pert$p.value,3)), filestats, ncol=3, append=T)
wilcox.test(gen_start, gen_end_pert, paired=T) -> gen_pert; write(c(gen_pert$data.name, gen_pert$statistic, signif(gen_pert$p.value,3)), filestats, ncol=3, append=T)
wilcox.test(env_end, env_start, paired=T) -> env_cont; write(c(env_cont$data.name, env_cont$statistic, signif(env_cont$p.value,3)), filestats, ncol=3, append=T)
wilcox.test(gen_end, gen_start, paired=T) -> gen_cont; write(c(gen_cont$data.name, gen_cont$statistic, signif(gen_cont$p.value,3)), filestats, ncol=3, append=T)


#Scatterplot
pdf("Robustness_scatterplot.pdf", width=8, height=7)
par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)

setwd(treatment_dir[1])
mydat <- read.table("data_summary.txt", header=T)
start <- subset(mydat, mydat$Generation == 0)
end <- subset(mydat, mydat$Generation == max(mydat$Generation))
setwd(treatment_dir[2])
mydat2 <- read.table("data_summary.txt", header=T)
end_pert <- subset(mydat2, mydat2$Generation == max(mydat2$Generation))

setwd(main_dir)
plot(start$Median_Robustness~start$Median_Env_Robustness, pch=1,
     xlab="Environmental Robustness", ylab="Genetic Robustness", 
     ylim=c(0,1))
for(i in 1:100){
  segments(start$Median_Env_Robustness[i],start$Median_Robustness[i],end$Median_Env_Robustness[i],end$Median_Robustness[i], col=2)
  segments(start$Median_Env_Robustness[i],start$Median_Robustness[i],end_pert$Median_Env_Robustness[i],end_pert$Median_Robustness[i], col=4)
}
points(start$Median_Robustness~start$Median_Env_Robustness, pch=16, col="white")
points(end$Median_Robustness~end$Median_Env_Robustness, pch=16, col="white")
points(end_pert$Median_Robustness~end_pert$Median_Env_Robustness, pch=16, col="white")
points(start$Median_Robustness~start$Median_Env_Robustness, pch=1, col=1)
points(end$Median_Robustness~end$Median_Env_Robustness, pch=1, col=2)
points(end_pert$Median_Robustness~end_pert$Median_Env_Robustness, pch=1,col=4)

legend("topleft", c("Start", "End - control", "End - perturbed"),
       lty=c(NA,1,1), lwd=1.5, pt.cex=1,
       pch=21, col=c(1,2,4), pt.bg=c("white","white","white"),
       bty="n", cex=1)

dev.off()


#-----------
# Figure 3
#-----------

setwd(main_dir)

#function to plot error bars
error.bars <- function(x, y, err, barwidth, col){ 
  segments(x,y-err,x,y+err, col=col); segments(x-barwidth,y+err,x+barwidth,y+err, col=col); segments(x-barwidth,y-err,x+barwidth,y-err, col=col)
}

starting.plot <- function(x, ylab, y.adj=0.05, log=FALSE){
  xlim = c(0,500)
  if(log){
    xlim=c(0, log10(500))
  }
  plot(NA,xlab="Generation",ylab=ylab, ylim=c(min(x)-y.adj, 
                                              max(x)+y.adj), xlim=xlim, xaxt="n")
  if(log){
    axis(1, at=c(log10(1), log10(10), log10(100), log10(200), log10(500)), labels=c(1,10,100,200,500))
  }
  else{
    axis(1)
  }
}

plot.points <- function(param, eb=TRUE, legend="bottomright", log=FALSE){
  treatments <- c("control", "perturbed")
  colors <- c(1,"dodgerblue")
  for (i in 1:2) {
    setwd(treatment_dir[i]) 
    mydat = read.table("data_summary.txt", header=T)
    mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
    for (gen in unique(mydat$Generation)) {
      dataset <- subset(mydat,Generation == gen)
      if(log){points(log10(gen),mean(dataset[[param]]),pch=16,col=colors[i])}
      else{points(gen,mean(dataset[[param]]),pch=16,col=colors[i])}
      if(eb){
        if(log){
          error.bars(log10(gen), mean(dataset[[param]]), 1.96*sd(dataset[[param]])/sqrt(length(dataset[[param]])), .02, colors[i])}
        else{error.bars(gen, mean(dataset[[param]]), 1.96*sd(dataset[[param]])/sqrt(length(dataset[[param]])), 2, colors[i])}
      }
    }
  }
  legend(legend,treatments,col=colors,pch=16, bty='n')
}

plot.points.anc <- function(param, param2, eb=TRUE, legend="bottomright", log=FALSE){
  treatments <- c("control", "perturbed (perturbed env.)", "perturbed (ancestral env.)")
  colors <- c(1,"dodgerblue", "grey70")
  for (i in 1:2) {
    setwd(treatment_dir[i]) 
    mydat = read.table("data_summary.txt", header=T)
    mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
    for (gen in unique(mydat$Generation)) {
      dataset <- subset(mydat,Generation == gen)
      if(log){points(log10(gen),mean(dataset[[param]]),pch=16,col=colors[i])}#points(gen,mean(dataset[[param]]),pch=1)
      else{points(gen,mean(dataset[[param]]),pch=16,col=colors[i])}
      if(i==2){
        if(log){points(log10(gen),mean(dataset[[param2]]),pch=16,col=colors[i+1])}#points(gen,mean(dataset[[param2]]),pch=1)
        else{points(gen,mean(dataset[[param2]]),pch=16,col=colors[i+1])}
      }
      if(eb){
        if(log){error.bars(log10(gen), mean(dataset[[param]]), 1.96*sd(dataset[[param]])/sqrt(length(dataset[[param]])), .02, colors[i])}
        else{error.bars(gen, mean(dataset[[param]]), 1.96*sd(dataset[[param]])/sqrt(length(dataset[[param]])), 2, colors[i])}
        if(i==2){
          if(log){error.bars(log10(gen), mean(dataset[[param2]]), 1.96*sd(dataset[[param2]])/sqrt(length(dataset[[param2]])), .02, colors[i+1])}
          else{error.bars(gen, mean(dataset[[param2]]), 1.96*sd(dataset[[param2]])/sqrt(length(dataset[[param2]])), 2, colors[i+1])}
        }
      }
    }
  }
  legend(legend,treatments,col=colors,pch=16, bty='n')
}

plot.lines <- function(param, eb){
  colors <- c(1,"dodgerblue")
  setwd(treatment_dir[2]) 
  mydat = read.table("data_summary.txt", header=T)
  for (i in 1:length(param)){
    for (gen in unique(mydat$Generation)) {
      dataset <- subset(mydat,Generation == gen)
      lines(gen,mean(dataset[[param[i]]]),pch=16,col=colors[i])
      if(eb){
        error.bars(gen, mean(dataset[[param[i]]]), 1.96*sd(dataset[[param[i]]])/sqrt(length(dataset[[param[i]]])), 2, colors[i])}
    }
  }
  legend("bottomright",param,col=colors,pch=16, bty='n')
}

setwd(treatment_dir[2]) 
temp.dat <- read.table("data_summary.txt", header=T)
temp.dat$mean_weight_offdiagonal <- temp.dat$mean_weight_all - temp.dat$mean_weight_diagonal

setwd(main_dir)
pdf("Proximate_mechanisms_log.pdf", width=10, height=6)
par(mfrow=c(2,3))

#Fitness
starting.plot(c(0.35,.95), "Mean Fitness", log=TRUE)
plot.points("Mean_Fitness", log=TRUE)

#Genetic Robustness
starting.plot(c(0.35,.95), "Median Genetic Robustness",log=TRUE)
plot.points.anc("Median_Robustness", "genetic_robustness_ancestral", legend="bottomright", log=TRUE)

#Environmental Robustness
starting.plot(c(0.35,.95), "Median Environmental Robustness", log=TRUE)
plot.points("Median_Env_Robustness", log=TRUE)

#Path length
starting.plot(c(0,55), "Path Length", log=TRUE)
plot.points.anc("path_length","path_length_ancestral", legend="topright", log=TRUE)

#Matrix asymmetry
#starting.plot(temp.dat$matrix_asymmetry, "Matrix Asymmetry", .005)
#plot.points("matrix_asymmetry")

#Nonzero diagonal
#starting.plot(temp.dat$nonzero_diagonal, "Nonzero Diagonal")
#plot.points("nonzero_diagonal")

#Mean weight of the diagonals
#starting.plot(c(.2,.8), "Mean Weight of the Diagonal", .005, log=TRUE)
#plot.points("mean_weight_diagonal", legend="topright",log=TRUE)

starting.plot(c(.2,.5), "Strength of Autoregulation", .005, log=TRUE)
plot.points("mean_weight_diagonal_nonzero", legend="topright",log=TRUE)

#starting.plot(c(.2,.8), "Proportion positive on diagonal", .005, log=TRUE)
#plot.points("prop_positive_on_diagonal", legend="topright",log=TRUE)

#Mean weight of the off-diagonals
starting.plot(c(-.3,0), "Mean Weight of Off-Diagonal Interactions", .005, log=TRUE)
plot.points("mean_weight_offdiagonal", legend="topright", log=TRUE)

#Mean weight all
#starting.plot(temp.dat$mean_weight_all, "Matrix Weight All")
#plot.points("mean_weight_all")

#Standard deviation all
#starting.plot(temp.dat$standard_deviation_all, "Standard Deviation", 0.005)
#plot.points("standard_deviation_all")

dev.off()

#t-test on mean weight diagonal endpoints
setwd(treatment_dir[1]) 
mydat = read.table("data_summary.txt", header=T)
mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
end.on <- mydat$mean_weight_diagonal[mydat$Generation==max(mydat$Generation)]
end.nz.on <- mydat$mean_weight_diagonal_nonzero[mydat$Generation==max(mydat$Generation)]
end.np.on <- mydat$prop_positive_on_diagonal[mydat$Generation==max(mydat$Generation)]
end.off <- mydat$mean_weight_offdiagonal[mydat$Generation==max(mydat$Generation)]
PL.cont <- mydat$path_length[mydat$Generation==unique(mydat$Generation)[6]]
setwd(treatment_dir[2]) 
mydat = read.table("data_summary.txt", header=T)
mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
end_pert.on <- mydat$mean_weight_diagonal[mydat$Generation==max(mydat$Generation)]
end_pert.nz.on <- mydat$mean_weight_diagonal_nonzero[mydat$Generation==max(mydat$Generation)]
end_pert.np.on <- mydat$prop_positive_on_diagonal[mydat$Generation==max(mydat$Generation)]
end_pert.off <- mydat$mean_weight_offdiagonal[mydat$Generation==max(mydat$Generation)]
PLA.pert <-  mydat$path_length_ancestral[mydat$Generation==unique(mydat$Generation)[6]]

setwd(main_dir)
write("Weight on/off diagonal endpoint t-tests", filestats, append=T)
t.test(end.on, end_pert.on, paired=TRUE) -> on_diagonal; write(c("mean weight diagonal", on_diagonal$statistic, signif(on_diagonal$p.value,3)), filestats, ncol=3, append=T)
t.test(end.nz.on, end_pert.nz.on, paired=TRUE) -> on_diag_nz; write(c("mean weight diagonal nonzero", on_diag_nz$statistic, signif(on_diag_nz$p.value,3)), filestats, ncol=3, append=T)
t.test(end.np.on, end_pert.np.on, paired=TRUE) -> on_diag_np; write(c("proportion positive on diagonal", on_diag_np$statistic, signif(on_diag_np$p.value,3)), filestats, ncol=3, append=T)
t.test(end.off,end_pert.off,paired=TRUE) -> off_diagonal; write(c("mean weight off diagonal", off_diagonal$statistic, signif(off_diagonal$p.value,3)), filestats, ncol=3, append=T)
t.test(PL.cont,PLA.pert,paired=TRUE) -> path_length; write(c("path length", path_length$statistic, signif(path_length$p.value,3)), filestats, ncol=3, append=T)

#-----------
# Figure 4
#-----------

setwd(main_dir)

pdf("Genetic_Variation_Boxplots.pdf", width=5, height=7)
par(mfrow=c(1,1))
par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)

#get data needed to set up empty plot
setwd(treatment_dir[1])
mydat <- read.table("data_summary.txt", header=T)
control_end <- subset(mydat$genetic_variation_euclidian, mydat$Generation == max(mydat$Generation))
setwd(treatment_dir[2])
mydat2 <- read.table("data_summary.txt", header=T)
perturb_end <- subset(mydat2$genetic_variation_euclidian, mydat2$Generation == max(mydat2$Generation))
setwd(main_dir)

boxplot(control_end, perturb_end, names=c("control", "experimental"), ylim=c(0,.07), col="grey70", ylab="Pairwise genetic distance", xlab="Environments")
dev.off()

write("Genetic variation endpoint population wilcoxon rank test", filestats, append=T)
wilcox.test(control_end, perturb_end, paired=T) -> gen_var; write(c(gen_var$data.name, gen_var$statistic, signif(gen_var$p.value,3)), filestats, ncol=3, append=T) 



##### Random Individuals Plots

setwd("~/GitHub/gene_network/data/Manuscript data/Figure 3/boolean_1500")
mydata = read.table("1500_random_individuals.txt", header=T)
mydata$mean_weight_offdiagonal <- mydata$mean_weight_all - mydata$mean_weight_diagonal
model.g <- function(param){
  g.lm <- lm(asin(genetic_robustness)~param, data=mydata); return(g.lm)}
model.e <- function(param){
  e.lm <- lm(asin(environmental_robustness)~param, data=mydata)}

g.lm <- model.g(mydata$mean_weight_diagonal); e.lm <- model.e(mydata$mean_weight_diagonal); mean1 <- mean(mydata$mean_weight_diagonal); sd1 <- sd(mydata$mean_weight_diagonal)
g.lm2 <- model.g(mydata$path_length); e.lm2 <- model.e(mydata$path_length); mean2 <- mean(mydata$path_length); sd2 <- sd(mydata$path_length)
#g.lm3 <- model.g(mydata$nonzero_diagonal); e.lm3 <- model.e(mydata$nonzero_diagonal); mean3 <- mean(mydata$nonzero_diagonal); sd3 <- sd(mydata$nonzero_diagonal)
#g.lm4 <- model.g(mydata$mean_weight_all); e.lm4 <- model.e(mydata$mean_weight_all); mean4 <- mean(mydata$mean_weight_all); sd4 <- sd(mydata$mean_weight_all)
#g.lm5 <- model.g(mydata$matrix_asymmetry); e.lm5 <- model.e(mydata$matrix_asymmetry); mean5 <- mean(mydata$matrix_asymmetry); sd5 <- sd(mydata$matrix_asymmetry)
#g.lm6 <- model.g(mydata$standard_deviation_all); e.lm6 <- model.e(mydata$standard_deviation_all); mean6 <- mean(mydata$standard_deviation_all); sd6 <- sd(mydata$standard_deviation_all)
g.lm7 <- model.g(mydata$mean_weight_offdiagonal); e.lm7 <- model.e(mydata$mean_weight_offdiagonal); mean7 <- mean(mydata$mean_weight_offdiagonal); sd7 <- sd(mydata$mean_weight_offdiagonal)

model.y <- function(model,x){
  yg <- x*coef(model)[[2]] + coef(model)[[1]]
  return(yg)
}

z <- seq(-3,3,.1)

x.vals <- function(mu, sd){
  z <- seq(-3,3,.1)
  x <- z*sd + mu
  return(x)
}

#x <- seq(-1,1,.01)
yg <- model.y(g.lm,x.vals(mean1,sd1)); ye <- model.y(e.lm,x.vals(mean1,sd1))
yg2 <- model.y(g.lm2,x.vals(mean2,sd2)); ye2 <- model.y(e.lm2,x.vals(mean2,sd2))
#yg3 <- model.y(g.lm3,x.vals(mean3,sd3)); ye3 <- model.y(e.lm3,x.vals(mean3,sd3))
#yg4 <- model.y(g.lm4,x.vals(mean4,sd4)); ye4 <- model.y(e.lm4,x.vals(mean4,sd4))
#yg5 <- model.y(g.lm5,x.vals(mean5,sd5)); ye5 <- model.y(e.lm5,x.vals(mean5,sd5))
#yg6 <- model.y(g.lm6,x.vals(mean6,sd6)); ye6 <- model.y(e.lm6,x.vals(mean6,sd6))
yg7 <- model.y(g.lm7,x.vals(mean7,sd7)); ye7 <- model.y(e.lm7,x.vals(mean7,sd7))

pdf("Correlation_plots.pdf", width=10, height=6)
par(mfrow=c(1,2))
plot(NA, ylim=c(0.2,1), xlim=c(-3,3), xlab="Z-value", ylab="Genetic robustness")
lines(sin(yg)~z, col=2, lwd=2)
lines(sin(yg2)~z, col=4, lwd=3)
#lines(sin(yg4)~z, col=5, lwd=2)
lines(sin(yg7)~z, col=3, lwd=2)
#lines(sin(yg6)~z, col=6, lwd=2)
legend("bottomright", c("Strength of autoregulation", "Off-diagonal interactions", "Path length"), lty=1, col=c(2,3,4), lwd=2, bty='n')
plot(NA, ylim=c(0.2,1), xlim=c(-3,3), xlab="Z-value", ylab="Environmental robustness")
lines(sin(ye)~z, col=2, lwd=2)
lines(sin(ye2)~z, col=4, lwd=2)
#lines(sin(ye4)~z, col=5, lwd=2)
lines(sin(ye7)~z, col=3, lwd=2)
#lines(sin(ye6)~z, col=6, lwd=2)
legend("bottomright", c("Strength of autoregulation", "Off-diagonal interactions", "Path length"), lty=1, col=c(2,3,4), lwd=2, bty='n')
dev.off()