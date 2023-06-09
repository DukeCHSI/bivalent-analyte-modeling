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

for (well_idx in 1:nwells){
  bivalentShort_pars <- rbind(bivalentShort_pars, fits_list[[well_idx]][["R0"]][["par"]][c(4)])
}

bivalentShort_pars <- data.frame(bivalentShort_pars)

colnames(bivalentShort_pars) <- c("kd2")
bivalentShort_pars <- melt(bivalentShort_pars)
colnames(bivalentShort_pars) <- c("ParamName","Value")
bivalentShort_pars$ParamName <- as.factor(bivalentShort_pars$ParamName)

plot <- ggplot(bivalentShort_pars, aes(x=ParamName, y=Value, color=ParamName)) + 
  # ggtitle(title) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
  geom_boxplot(width=0.1, color='black') +
  scale_y_continuous(limits = c(NA, 1.5e-4)) +
  theme(legend.position="none", text = element_text(size=20),
        axis.title.x=element_blank())
plot

#####

load(file='bivalent_oldResults_t0_long.Rdata')
nwells <- length(fits_list)

bivalentShort_pars <- NULL

for (well_idx in 1:nwells){
  bivalentShort_pars <- rbind(bivalentShort_pars, fits_list[[well_idx]][["R0"]][["par"]][c(4)])
}

bivalentShort_pars <- data.frame(bivalentShort_pars)

colnames(bivalentShort_pars) <- c("kd2")
bivalentShort_pars <- melt(bivalentShort_pars)
colnames(bivalentShort_pars) <- c("ParamName","Value")
bivalentShort_pars$ParamName <- as.factor(bivalentShort_pars$ParamName)

plot <- ggplot(bivalentShort_pars, aes(x=ParamName, y=Value, color=ParamName)) + 
  # ggtitle(title) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="black") +
  geom_boxplot(width=0.1, color='black') +
  scale_y_continuous(limits = c(NA, 1.5e-4)) +
  theme(legend.position="none", text = element_text(size=20),
        axis.title.x=element_blank())
plot
