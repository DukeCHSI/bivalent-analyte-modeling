library(readr)
library(tidyverse)
library(readxl)
library(minpack.lm)
library(zoo)
library(gridExtra)
library(grid)
library(kableExtra)
library(deSolve)
library(reshape)

setwd("~/Summer2022/project/hiv-summer-2022/")
load(file='bivalent_oldResults_t0_short2_nonRegen.Rdata')
nwells <- length(fits_list)

bivalentShort_pars <- NULL
well_idc <- 1:14 #list(5,6,7,10,11,12,13)
for (well_idx in well_idc){
  # if (!(well_idx == 9) && !(well_idx == 10)){
    bivalentShort_pars <- rbind(bivalentShort_pars, fits_list[[well_idx]][["R0"]][["par"]][c(1,2,3,4)])
  # }
}

bivalentShort_pars <- data.frame(bivalentShort_pars)

colnames(bivalentShort_pars) <- c("ka1","ka2","kd1","kd2")

meanShort_ka1 <- mean(bivalentShort_pars$ka1)
meanShort_ka2 <- mean(bivalentShort_pars$ka2)
meanShort_kd1 <- mean(bivalentShort_pars$kd1)
meanShort_kd2 <- mean(bivalentShort_pars$kd2)

sdShort_ka1 <- sd(bivalentShort_pars$ka1)
sdShort_ka2 <- sd(bivalentShort_pars$ka2)
sdShort_kd1 <- sd(bivalentShort_pars$kd1)
sdShort_kd2 <- sd(bivalentShort_pars$kd2)

bivalentShort_pars <- melt(bivalentShort_pars)
colnames(bivalentShort_pars) <- c("ParamName","Value")
bivalentShort_pars$ParamName <- as.factor(bivalentShort_pars$ParamName)

title <- paste("Estimated parameters for standard length of dissociation")

plot <- ggplot(bivalentShort_pars, aes(x=ParamName, y=Value, color=ParamName)) + 
  # ggtitle(title) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
  geom_boxplot(width=0.1, color='black') +
  scale_y_continuous(trans = "log10", breaks = c(1e-7, 1e-4, 1e-1, 1e2, 1e5), limits = c(1e-7, 1e5)) +
  # theme(legend.position="none", text = element_text(size=20),
  #       axis.title.x=element_blank())
plot

cvShort_ka1 <- sdShort_ka1/meanShort_ka1
cvShort_ka2 <- sdShort_ka2/meanShort_ka2
cvShort_kd1 <- sdShort_kd1/meanShort_kd1
cvShort_kd2 <- sdShort_kd2/meanShort_kd2
#####

load(file='bivalent_oldResults_t0_long.Rdata')
nwells <- length(fits_list)

bivalentLong_pars <- NULL

for (well_idx in 1:nwells){
  # if (!(well_idx == 9) && !(well_idx == 10)){
    bivalentLong_pars <- rbind(bivalentLong_pars, fits_list[[well_idx]][["R0"]][["par"]][c(1,2,3,4)])
  # }
}

bivalentLong_pars <- data.frame(bivalentLong_pars)

colnames(bivalentLong_pars) <- c("ka1","ka2","kd1","kd2")

meanLong_ka1 <- mean(bivalentLong_pars$ka1)
meanLong_ka2 <- mean(bivalentLong_pars$ka2)
meanLong_kd1 <- mean(bivalentLong_pars$kd1)
meanLong_kd2 <- mean(bivalentLong_pars$kd2)

sdLong_ka1 <- sd(bivalentLong_pars$ka1)
sdLong_ka2 <- sd(bivalentLong_pars$ka2)
sdLong_kd1 <- sd(bivalentLong_pars$kd1)
sdLong_kd2 <- sd(bivalentLong_pars$kd2)


bivalentLong_pars <- melt(bivalentLong_pars)
colnames(bivalentLong_pars) <- c("ParamName","Value")
bivalentLong_pars$ParamName <- as.factor(bivalentLong_pars$ParamName)

title <- paste("Estimated parameters for extended length of dissociation")
plot <- ggplot(bivalentLong_pars, aes(x=ParamName, y=Value, color=ParamName)) + 
  # ggtitle(title) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
  geom_boxplot(width=0.1, color='black') +
  scale_y_continuous(trans = "log10", breaks = c(1e-7, 1e-4, 1e-1, 1e2, 1e5), limits = c(1e-7, 1e5)) +
  theme(legend.position="none", text = element_text(size=20),
        axis.title.x=element_blank())
plot

cvLong_ka1 <- sdLong_ka1/meanLong_ka1
cvLong_ka2 <- sdLong_ka2/meanLong_ka2
cvLong_kd1 <- sdLong_kd1/meanLong_kd1
cvLong_kd2 <- sdLong_kd2/meanLong_kd2
