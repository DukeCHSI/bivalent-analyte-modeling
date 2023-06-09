library(readr)
library(tidyverse)
library(readxl)
library(minpack.lm)
library(zoo)
library(gridExtra)
library(grid)
library(kableExtra)
library(deSolve)
library(parallel)
library(svDialogs)

num_cores <- detectCores()/2
options("scipen"=2, "digits"=3)

setwd("~/Summer2022/project/hiv-summer-2022/")

data_file_path <- "paper_data/Downselect 14 sensograms Carterra fits exported data.xlsx"
titration_data <- read_excel(data_file_path, col_names = FALSE, skip = 1, n_max = 1000)

concs <- rbind(1e-6, 5e-7, 2.5e-7, 1.25e-7, 6.25e-8)

nunCol_per_sensorgram <- 20
n_row <- dim(titration_data)[1]
n_col <- dim(titration_data)[2]
num_sensorgrams <- n_col/(nunCol_per_sensorgram)
dissociation_start <- 420
for (i in 1:num_sensorgrams){
  # sensorgram_data <- titration_data[,((i-1)*nunCol_per_sensorgram+1):(i*nunCol_per_sensorgram)]
  # s_n_cols <- dim(sensorgram_data)[2]
  RU_idc <- seq(((i-1)*nunCol_per_sensorgram+1),(i*nunCol_per_sensorgram), 2)
  t_idc <- seq(((i-1)*nunCol_per_sensorgram+2),(i*nunCol_per_sensorgram), 2)
  data_RU_idc <- RU_idc[1:5]
  data_t_idc <- t_idc[1:5]
  fit_RU_idc <- RU_idc[6:10]
  fit_t_idc <- t_idc[6:10]
  
  data_RU <- titration_data[,data_RU_idc]
  data_t <- titration_data[,data_t_idc]
  fit_RU <- titration_data[,fit_RU_idc]
  fit_t <- titration_data[,fit_t_idc]
  
  data_Time_df <- NULL
  data_RU_df <- NULL
  data_concentrations <- NULL
  
  fit_Time_df <- NULL
  fit_RU_df <- NULL
  fit_concentrations <- NULL
  
  data_df <- NULL
  fit_df <- NULL
  
  for (ic in 1:length(concs)){
    data_Time_vals <- na.omit(data.frame(data_t[,ic]))
    data_RU_vals <- na.omit(data.frame(data_RU[,ic]))
    data_Time_df <- bind_rows(data_Time_df, data_Time_vals)
    data_RU_df <- bind_rows(data_RU_df, data_RU_vals)
    
    incl_concentrations <- concs[ic]
    n_time_points <- dim(data_Time_vals)[1]
    data_concentrations <- bind_rows(data_concentrations, data.frame(rep(incl_concentrations, n_time_points)))
    
    fit_Time_vals <- na.omit(data.frame(fit_t[,ic]))
    fit_RU_vals <- na.omit(data.frame(fit_RU[,ic]))
    fit_Time_df <- bind_rows(fit_Time_df, fit_Time_vals)
    fit_RU_df <- bind_rows(fit_RU_df, fit_RU_vals)
    
    incl_concentrations <- concs[ic]
    n_time_points <- dim(fit_Time_vals)[1]
    fit_concentrations <- bind_rows(fit_concentrations, data.frame(rep(incl_concentrations, n_time_points)))
    
    # colnames(data_Time_vals) <- "Time"
    # colnames(fit_Time_vals) <- "Time"
    # intersect(data_Time_vals, fit_Time_vals)
    # common_time_vals <- intersect(data_Time_vals, fit_Time_vals)
  }
  data_df <- cbind(data_Time_df, data_RU_df, data_concentrations)
  colnames(data_df) <- c("Time", "RU", "Concentration")
  data_df$Concentration <- as_factor(data_df$Concentration)

  fit_df <- cbind(fit_Time_df, fit_RU_df, fit_concentrations)
  colnames(fit_df) <- c("Time", "RU", "Concentration")

  plot <- ggplot(data_df, aes(x = Time, y = RU)) + geom_point(size = 1, aes(color = Concentration)) +
    # ggtitle(ligand_desc, subtitle = sub_title) +
    guides(color = guide_legend(reverse=TRUE, override.aes = list(size = 5))) +
    geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1.5) +
    geom_vline(xintercept = dissociation_start, linetype="dashed", color = "black", size=1) +
    theme(text = element_text(size=15)) +
      labs(color = "Concentration (M)") +
      ylab("Response Unit (RU)") + xlab("Time (s)")

  data_file_name <- paste("paper_codes/figures/Carterra/Carterrafit_long",i,".png",sep="")
  ggsave(filename = data_file_name, plot=plot, width = 12, height = 7.31)
}

