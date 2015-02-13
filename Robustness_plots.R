#-----------
# Set up functions for data summary and plotting
#-----------

library("gplots")
library("MASS")

contour.plot <- function(filename, wd=getwd()){
  dat <- read.table(filename, header=T)
  dat$g_rob <- round(dat$genetic_robustness,1)
  dat$e_rob <- round(dat$environmental_robustness,1)
  dat.tab2 <- as.matrix(table(dat$e_rob, dat$g_rob))
  colors = c("white","grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10","black","black")
  setwd(wd)
  pdf("Contour_Plot.pdf", width=7, height=7)
  par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)
  filled.contour(dat.tab2, col=colors, nlevels=9, xlab="environmental robustness",
               ylab="genetic robustness")
  dev.off()
  out1 <- cor.test(dat$environmental_robustness, dat$genetic_robustness, method="spearman")
  out1
}


data.summary <- function(treatment_dir){
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

robustness.boxplots <- function(treatment_dir, wd=getwd(), stats=FALSE){
  setwd(wd)
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

  if(stats){
    #Wilcoxon sign-rank test
    #Does not assume a normal distribution
    setwd(wd)
    filestats="Population_stats.txt"
    fileHeader = c("Statistical_text", "Comparison", "stat_value", "p_value")
    nCol = length(fileHeader)
    write(fileHeader, filestats, ncolumns=nCol)
    wilcox.test.write(env_start, env_end_pert, "Env. robustness start vs. end perturb", filestats, nCol)
    wilcox.test.write(env_start, env_end, "Env. robustness start vs. end control", filestats, nCol)
    wilcox.test.write(gen_start, gen_end_pert, "Gen. robustness start vs. end perturb", filestats, nCol)
    wilcox.test.write(gen_start, gen_end, "Gen. robustness start vs. end control", filestats, nCol)
    }
}

wilcox.test.write <- function(pop1, pop2, comparison, filename, ncol){
  wt <- wilcox.test(pop1,pop2, paired=T)
  write(c("Wilcox sign-rank test", comparison, wt$statistic, signif(wt$p.value,3)), filename, ncol=ncol, append=T)
}

robustness.scatterplot <- function(treatment_dir, wd=setwd(), lines=TRUE, legend="topleft")
  pdf("Robustness_scatterplot.pdf", width=8, height=7)
  par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)

  setwd(treatment_dir[1])
  mydat <- read.table("data_summary.txt", header=T)
  start <- subset(mydat, mydat$Generation == 0)
  end <- subset(mydat, mydat$Generation == max(mydat$Generation))
  setwd(treatment_dir[2])
  mydat2 <- read.table("data_summary.txt", header=T)
  end_pert <- subset(mydat2, mydat2$Generation == max(mydat2$Generation))

  setwd(wd)
  plot(start$Median_Robustness~start$Median_Env_Robustness, pch=1,
      xlab="Environmental Robustness", ylab="Genetic Robustness", 
      ylim=c(0,1))
  if(lines){
      segments(start$Median_Env_Robustness,start$Median_Robustness,end$Median_Env_Robustness,end$Median_Robustness, col=2)
      segments(start$Median_Env_Robustness,start$Median_Robustness,end_pert$Median_Env_Robustness,end_pert$Median_Robustness, col=4)
    }
  points(start$Median_Robustness~start$Median_Env_Robustness, pch=16, col="white")
  points(end$Median_Robustness~end$Median_Env_Robustness, pch=16, col="white")
  points(end_pert$Median_Robustness~end_pert$Median_Env_Robustness, pch=16, col="white")
  points(start$Median_Robustness~start$Median_Env_Robustness, pch=1, col=1)
  points(end$Median_Robustness~end$Median_Env_Robustness, pch=1, col=2)
  points(end_pert$Median_Robustness~end_pert$Median_Env_Robustness, pch=1,col=4)

  legend(legend, c("Start", "End - control", "End - perturbed"),
       lty=c(NA,1,1), lwd=1.5, pt.cex=1,
       pch=21, col=c(1,2,4), pt.bg=c("white","white","white"),
       bty="n", cex=1)

  dev.off()
}

#function to plot error bars
error.bars <- function(x, y, err, barwidth, col){ 
  segments(x,y-err,x,y+err, col=col); segments(x-barwidth,y+err,x+barwidth,y+err, col=col); segments(x-barwidth,y-err,x+barwidth,y-err, col=col)
}

starting.plot <- function(y, ylab, gens, y.adj=0, log=FALSE){
  if(log){
    xlim=c(0, log10(gens))
  }
  else{
    xlim=c(0,gens)
  }
  plot(NA,xlab="Generation",ylab=ylab, ylim=c(min(y)-y.adj, 
                                              max(y)+y.adj), xlim=xlim, xaxt="n")
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

proximate.mechanisms <- function(treatment_dir, wd=getwd(), log=TRUE){
  setwd(treatment_dir[2]) 
  temp.dat <- read.table("data_summary.txt", header=T)
  temp.dat$mean_weight_offdiagonal <- temp.dat$mean_weight_all - temp.dat$mean_weight_diagonal
  gens <- max(temp.dat$Generation)

  setwd(wd)
  pdf("Proximate_mechanisms_log.pdf", width=10, height=6)
  par(mfrow=c(2,3))

  #Fitness
  starting.plot(c(0.3,1), "Mean Fitness", gens, log=log)
  plot.points("Mean_Fitness", log=log)

  #Genetic Robustness
  starting.plot(c(0.3,1), "Median Genetic Robustness", gens, log=log)
  plot.points.anc("Median_Robustness", "genetic_robustness_ancestral", legend="bottomright", log=log)

  #Environmental Robustness
  starting.plot(c(0.3,1), "Median Environmental Robustness", gens, log=log)
  plot.points("Median_Env_Robustness", log=log)

  #Path length
  starting.plot(c(0,55), "Path Length", gens, log=log)
  plot.points.anc("path_length","path_length_ancestral", legend="topright", log=log)

  #Strength of autoregulation
  starting.plot(temp.dat$mean_weight_diagonal_nonzero, "Strength of Autoregulation", gens, y.adj=0.2, log=log)
  plot.points("mean_weight_diagonal_nonzero", legend="topright",log=log)

  #Mean weight of the off-diagonals
  starting.plot(temp.dat$mean_weight_offdiagonal, "Mean Weight of Off-Diagonal Interactions", gens, y.adj=0.2, log=log)
  plot.points("mean_weight_offdiagonal", legend="topright", log=log)

  dev.off()
}

proximate.mechanism.stats <- function(treatment_dir, wd=getwd()){
  #t-test on mean weight diagonal endpoints
  setwd(treatment_dir[1]) 
  mydat = read.table("data_summary.txt", header=T)
  mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
  pl.gen <- pl.difference(treatment_dir)
  
  end.on <- mydat$mean_weight_diagonal[mydat$Generation==max(mydat$Generation)]
  end.nz.on <- mydat$mean_weight_diagonal_nonzero[mydat$Generation==max(mydat$Generation)]
  end.np.on <- mydat$prop_positive_on_diagonal[mydat$Generation==max(mydat$Generation)]
  end.off <- mydat$mean_weight_offdiagonal[mydat$Generation==max(mydat$Generation)]
  PL.cont <- mydat$path_length[mydat$Generation==unique(mydat$Generation)[pl.gen]]
  setwd(treatment_dir[2]) 
  mydat = read.table("data_summary.txt", header=T)
  mydat$mean_weight_offdiagonal <- mydat$mean_weight_all - mydat$mean_weight_diagonal
  end_pert.on <- mydat$mean_weight_diagonal[mydat$Generation==max(mydat$Generation)]
  end_pert.nz.on <- mydat$mean_weight_diagonal_nonzero[mydat$Generation==max(mydat$Generation)]
  end_pert.np.on <- mydat$prop_positive_on_diagonal[mydat$Generation==max(mydat$Generation)]
  end_pert.off <- mydat$mean_weight_offdiagonal[mydat$Generation==max(mydat$Generation)]
  PLA.pert <-  mydat$path_length_ancestral[mydat$Generation==unique(mydat$Generation)[pl.gen]]

  setwd(wd)
  filestats="Proximate_mech_stats.txt"
  fileHeader = c("Statistical_text", "Comparison", "stat_value", "p_value")
  nCol = length(fileHeader)
  write(fileHeader, filestats, ncolumns=nCol)
  t.test.write(end.nz.on, end_pert.nz.on, "Endpoint populations of autoregulation", filestats, nCol)
  t.test.write(end.off,end_pert.off, "Endpoint populations of off-diagonal", filestats, nCol)
  t.test.write(PL.cont,PLA.pert, "Path length at largest difference", filestats, nCol)
}

pl.difference <- function(treatment_dir){
  setwd(treatment_dir[1]) 
  mydat = read.table("data_summary.txt", header=T)
  gens <- unique(mydat$Generation)[order(unique(mydat$Generation))]
  pl <- aggregate(.~Generation, data=mydat, FUN=mean)$path_length
  setwd(treatment_dir[2]) 
  mydat2 = read.table("data_summary.txt", header=T)
  pl2 <- aggregate(.~Generation, data=mydat, FUN=mean)$path_length_ancestral
  pl.diff <- pl - pl2
  return(which.max(abs(pl.diff)))
}

t.test.write <- function(pop1, pop2, comparison, filename, ncol){
  wt <- t.test(pop1,pop2, paired=T)
  write(c("paired t-test", comparison, wt$statistic, signif(wt$p.value,3)), filename, ncol=ncol, append=T)
}

genetic.variation <- function(treatment_dir, wd=getwd()){
  setwd(wd)
  pdf("Genetic_Variation_Boxplots.pdf", width=5, height=7)
  par(mfrow=c(1,1))
  par(mar=c(4,4,2.4,2.4), mgp=c(2.4, 0.8, 0), cex.lab=1.27, cex.axis=1.12)

  setwd(treatment_dir[1])
  mydat <- read.table("data_summary.txt", header=T)
  control_end <- subset(mydat$genetic_variation_euclidian, mydat$Generation == max(mydat$Generation))
  setwd(treatment_dir[2])
  mydat2 <- read.table("data_summary.txt", header=T)
  perturb_end <- subset(mydat2$genetic_variation_euclidian, mydat2$Generation == max(mydat2$Generation))
  setwd(main_dir)

  boxplot(control_end, perturb_end, names=c("control", "experimental"), ylim=c(0,.07), col="grey70", ylab="Pairwise genetic distance", xlab="Environments")
  dev.off()

  setwd(wd)
  filestats="Genetic_variation_stats.txt"
  fileHeader = c("Statistical_text", "Comparison", "stat_value", "p_value")
  nCol = length(fileHeader)
  write(fileHeader, filestats, ncolumns=nCol)
  wilcox.test.write(control_end, perturb_end, "Endpoing genetic variation (control vs. pert)", filestats, nCol)
}


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


##PLOT FIGURES

#Contour plot
setwd("~/GitHub/gene_network/data/Manuscript data/Figure 1")
contour.plot("10000_random_individuals.txt")

#Robustness boxplots
setwd("~/GitHub/gene_network/data/Manuscript data/Figure 2 new/default6")
main_dir = getwd()
treatment_dir <- NULL; treatment_dir[1] <- paste(main_dir,"/control_pop",sep=""); treatment_dir[2] <- paste(main_dir,"/perturb_pop",sep="")

data.summary(treatment_dir)
robustness.boxplots(treament_dir, main_dir, stats=TRUE)
robustness.scatterplot(treatment_dir, main_dir)

proximate.mechanisms(treatment_dir, main_dir)
proximate.mechanisms.stats(treatment_dir, main_dir)

