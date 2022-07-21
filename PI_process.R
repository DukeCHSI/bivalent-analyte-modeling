library(tidyverse)
library(readxl)
library(minpack.lm)
library(zoo)
library(gridExtra)
library(grid)
library(kableExtra)
library(parallel)
library(svDialogs)
library(deSolve)
setwd("~/Summer2022/project")

ka1 <- 1.4e4
ka2 <- 4.6e3
ka3 <- ka1
ka4 <- ka2
kd1 <- 4.9e-4
kd2 <- 3.4e-5
kd3 <- kd1
kd4 <- kd2

mid_conc1 <- kd1/ka1
mid_conc2 <- kd2/ka2

concs1 <- c(6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07, 1e-6)#c(mid_conc1/4, mid_conc1/2, mid_conc1, 2*mid_conc1, 4*mid_conc1)
concs2 <- c(6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07, 1e-6)#c(mid_conc2/4, mid_conc2/2, mid_conc2, 2*mid_conc2, 4*mid_conc2)
num_conc <- 5
Rmaxs <- c(120, 120, 120, 120, 120)
R0s <- c(0, 0, 0, 0, 0)

pars <- c(ka1, ka2, kd1, kd2, ka3, ka4, kd3, kd4, Rmaxs, R0s)
num_Rmax <- 5
num_R0 <- 5
param_names <- c("ka1", "ka2", "kd1", "kd2", "ka3", "ka4", "kd3", "kd4", rep("Rmax", num_Rmax), rep("R0", num_R0))

time_asc <- seq(from = 0, to = 240, by = 2)
time_dis <- seq(from = 242, to = 1000, by = 2)

time <- c(time_asc, time_dis)

use_regeneration <- 'N'
use_globalRmax <- 'N'
use_secondRebind <- 'N'








lb <- c(1e2, 1e2, 1e-7, 1e-7)
lb2 <- lb
ub <- c(1e7, 1e7, 1e-1, 1e-1)
ub2 <- ub

nonPI_pars <- c(ka1, ka2, kd1, kd2, Rmaxs, R0s)
nonPI_pars_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep("R0", num_R0))

varied_params <- c(ka3, ka4, kd3, kd4)
varied_params2 <- varied_params
num_sample <- 20



load("~/Summer2022/project/SSE_2Pairs.Rdata")

sampled_space <- 10^seq(from = log10(lb[1]), to = log10(ub[1]), length.out=num_sample)
sampled_space <- c(sampled_space, varied_params[1])
plot(log10(sampled_space), log10(SSE_noise[1,]), main="ka3",
     xlab="ka3", ylab="SSE", pch=19)

sampled_space <- 10^seq(from = log10(lb[2]), to = log10(ub[2]), length.out=num_sample)
sampled_space <- c(sampled_space, varied_params[2])
plot(log10(sampled_space), log10(SSE_noise[2,]), main="ka4",
     xlab="ka4", ylab="SSE", pch=19)

sampled_space <- 10^seq(from = log10(lb[3]), to = log10(ub[3]), length.out=num_sample)
sampled_space <- c(sampled_space, varied_params[3])
plot(log10(sampled_space), log10(SSE_noise[3,]), main="kd3",
     xlab="kd3", ylab="SSE", pch=19)

sampled_space <- 10^seq(from = log10(lb[4]), to = log10(ub[4]), length.out=num_sample)
sampled_space <- c(sampled_space, varied_params[4])
plot(log10(sampled_space), log10(SSE_noise[4,]), main="kd4",
     xlab="kd4", ylab="SSE", pch=19)




model <- 'monovalent'

output1 <- run_model(pars, time_asc, time_dis, num_conc, concs1, concs2, param_names, 
                     use_regeneration, 
                     model, 
                     use_globalRmax,
                     use_secondRebind)

model <- 'monovalent'

output2 <- run_model(pars, time_asc, time_dis, num_conc, concs2, concs1, param_names, 
                     use_regeneration, 
                     model, 
                     use_globalRmax,
                     use_secondRebind)


model <- 'biEpitopicLigandHeterogenousAnalyte'


output12 <- run_model(pars, time_asc, time_dis, num_conc, concs1, concs2, param_names, 
                      use_regeneration, 
                      model, 
                      use_globalRmax,
                      use_secondRebind)

df <- output12
# df %>% rename(RU12 = RU) -> df
# df %>% 
#   add_column(RU1 = output1$RU, .before = "Concentration1") %>%
#   add_column(RU2 = output2$RU, .before = "Concentration1") %>%
#   add_column(RU1P2 = output1$RU+output2$RU, .before = "Concentration1") -> df

plot <- ggplot() + 
  geom_line(data =  df, aes(x = Time, y = RU1, group = Concentration1), color = "black", size=1) +
  geom_line(data =  df, aes(x = Time, y = RU2, group = Concentration1), color = "grey", size=1) +
  geom_line(data =  df, aes(x = Time, y = RU12, group = Concentration2), color = "red", linetype = "dotdash", size=1)
  # scale_colour_manual("", breaks = c("New Model", "Sum of 1:1 Outputs"),
  #                     values = c("black", "red"))

plot
