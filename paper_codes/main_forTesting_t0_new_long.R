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

# select working directory
# files_directory  <-  rstudioapi::selectDirectory(
#   caption <- "Select Directory",
#   label <- "Select"
# )

setwd("~/Summer2022/project/hiv-summer-2022/")


# source("paper_codes/functions_forTesting_t0.R")
#source("~/Google Drive/My Drive/R/HIV_new_data/code/3odes.R")
#source("~/Google Drive/My Drive/R/HIV_new_data/code/RTFunctions_kyle3.R")

####### Input Files #############################################
#sample_sheet_path <- rstudioapi::selectFile(caption = "Select the sample sheet file",
#                               filter =
#                               existing = TRUE)

sample_sheet_path <- "paper_data/Bivalent Binding 2nd - CH505 CH31 - info sheet.csv"

sample_sheet <- read_csv(sample_sheet_path)

rows_dict <- list(A = 1, E = 2,
                  B = 3, F = 4,
                  C = 5, G = 6,
                  D = 7, H = 8)

sample_info <- process_sample_sheet(sample_sheet_path, rows_dict, sort_order = "pairrow")
#temp <- sample_info$Baseline
#sample_info$Baseline <- sample_info$`Bsl Start`
#sample_info$`Bsl Start` <- temp
#sample_info <- process_sample_sheet(sample_sheet_path)


#data_file_path <- rstudioapi::selectFile(caption = "Select the data file",
#                               filter = "All Files (*)",
#                               existing = TRUE)

data_file_path <- "paper_data/Bivalent Binding 2nd - CH505 CH31 - new output format.xlsx"

titration_data <- read_excel(data_file_path, col_names = TRUE, skip = 2, n_max = 1000)

####### User preferences #############################################

min_allowed_kd <- 10^(-5)
max_iterations <- 100

#min_RU <- 0

min_RU_tol <- 20
max_RU_tol <- 300

####### Process sample sheet #############################################

# There are different numbers of concentrations exported for each well. 
# Here, we delete the observations from the time series and from the ligand_conc data frame 
# that are chosen by the user as not to be included. 
# For some analyses, we will restrict to the chosen concentrations.

#remove wells from ligand each time series.

selected_samples <- select_samples(sample_info, titration_data)


# all time and ru vals for all selected wells
#x_vals <- x_vals[, keep_concentrations_all] 
#y_vals <- y_vals[, keep_concentrations_all]
# Now we only have the ligands chosen in the sample sheet, and we have recorded which concentrations
# to include in fitting

sample_info <- selected_samples$sample_info
keep_concentrations <- selected_samples$keep_concentrations
Time <- selected_samples$Time
RU <- selected_samples$RU
all_concentrations_values <- selected_samples$all_concentrations_values
n_time_points <- selected_samples$n_time_points

rm(selected_samples)

selected_concentrations <- select_concentrations(sample_info, Time, RU)

keep_concentrations <- selected_concentrations$keep_concentrations
sample_info <- selected_concentrations$sample_info
incl_concentrations_values <- selected_concentrations$incl_concentrations_values


####### Baseline correction #############################################

nwells <- dim(sample_info)[1]

# get index of first concentration (selected for analysis) for each well

first_incl_conc_idx_list <- map(.x = 1:nwells, .f = first_conc_indices,
                                sample_info$NumInclConc)
first_conc_idx_list <- map(.x = 1:nwells, .f = first_conc_indices,
                           sample_info$NumConc)

sample_info$FirstInclConcIdx <- as_vector(flatten(first_incl_conc_idx_list))
sample_info$FirstConcIdx <- as_vector(flatten(first_conc_idx_list))

# get baseline averages for each well (first time point)
sample_info$BaselineAverage <- rep(0, nwells)

baseline_info_list <- map_dfr(.x = 1:nwells, .f = get_baseline_indices2,
                              sample_info,
                              Time[, keep_concentrations], RU[, keep_concentrations])
########################################################### Why did we choose BaselineAverage = MinBaseline
sample_info$BaselineAverage <- baseline_info_list$MinBaseline
sample_info$BaselineIdx <- baseline_info_list$BaselineIdx
sample_info$BaselineNegative <- baseline_info_list$BaselineNegative
rm(baseline_info_list)
########################################################### 
# Kyle commented
#sample_info$WellIdx <- 1:nwells
#sample_info %>% arrange(Ligand) %>% select(WellIdx) -> sort_idx

#sort_idx <- as_vector(sort_idx)
########################################################### 
#want to sort output by ligand name. Need to print this list in order of sort_idx
#print these for all selected wells

plot_list_before_baseline <- map(.x = 1:nwells, .f = plot_sensorgrams, sample_info,
                                 Time, RU, 
                                 incl_concentrations_values,
                                 all_concentrations_values,
                                 n_time_points, all = TRUE)

# for (well_idx in 1:nwells){
#   print(well_idx)
#   if (!is.null(plot_list_before_baseline[[well_idx]])){
#     grid.arrange(plot_list_before_baseline[[well_idx]])
#     data_file_name <- paste("figures/data/data",well_idx,".png",sep="")
#     ggsave(filename = data_file_name, plot=plot_list_before_baseline[[well_idx]])
#   }
# }
# dev.off()

#plot_list_before_baseline_sorted <- plot_list_before_baseline[sort_idx]


# baseline correction - only selected concentrations

RU[, keep_concentrations] <- map_dfc(.x = 1:nwells, .f = baseline_correction,
                                     Time[, keep_concentrations],
                                     RU[, keep_concentrations],
                                     sample_info)

plot_list <- map(.x = 1:nwells, .f = plot_sensorgrams_with_nobaseline, sample_info, 
                 Time[, keep_concentrations],
                 RU[, keep_concentrations],
                 incl_concentrations_values,
                 all_concentrations_values,
                 n_time_points, all = FALSE)

# Haven't corrected for duplicated values in incl_conc (only the first to be used)

####### Fits #############################################

fixed_sheet <- NULL

## Fit using non-parallel computing
bivalent_fits_list <- map(.x = 1:nwells, .f = fit_association_dissociation,
                          sample_info,
                          Time[, keep_concentrations],
                          RU[, keep_concentrations],
                          incl_concentrations_values,
                          fixed_sheet,
                          min_allowed_kd,
                          max_iterations,
                          ptol,
                          ftol)

## Fit using parallel computing
fits_list <- mclapply(X = 1:nwells, FUN = fit_association_dissociation, mc.cores = num_cores, sample_info,
                      Time[, keep_concentrations],
                      RU[, keep_concentrations],
                      incl_concentrations_values,
                      fixed_sheet,
                      min_allowed_kd,
                      max_iterations,
                      ptol,
                      ftol)


# save(monovalent_fits_list, file="monovalent.Rdata")
save(fits_list, file="bivalent_newResults_t0_long.Rdata")
# load(file='NoBulkshiftShort.Rdata')
load(file='bivalent_newResults_t0_long.Rdata')
# load(file='bulkshift_long.Rdata')


####### Response curve #############################################
rc_list <- map(.x = 1:nwells, .f = get_response_curve, sample_info, 
               Time, RU,
               all_concentrations_values,
               incl_concentrations_values, n_time_points)
# # ## Bivalent
# plot_bivalent_fit <- map(.x = 1:nwells, .f = plot_sensorgrams_with_monovalent_fits,
#                          sample_info,
#                          fits_list,
#                          Time[, keep_concentrations], RU[, keep_concentrations],
#                          incl_concentrations_values,
#                          n_time_points)
# 
# for (well_idx in 1:nwells){
#   print(well_idx)
#   if (!is.null(plot_bivalent_fit[[well_idx]])){
#     grid.arrange(plot_bivalent_fit[[well_idx]])
#     data_file_name <- paste("figures/bulkshift_long/fitbulkshift_long",well_idx,".png",sep="")
#     ggsave(filename = data_file_name, plot=plot_bivalent_fit[[well_idx]])
#   }
# }
# dev.off()

# Bivalent
plot_bivalent_fit <- map(.x = 1:nwells, .f = plot_sensorgrams_with_fits,
                         sample_info,
                         fits_list,
                         Time[, keep_concentrations], RU[, keep_concentrations],
                         incl_concentrations_values,
                         n_time_points)
for (well_idx in 1:nwells){
  print(well_idx)
  if (!is.null(plot_bivalent_fit[[well_idx]])){
    grid.arrange(plot_bivalent_fit[[well_idx]])
    data_file_name <- paste("paper_codes/figures/bivalent_new_long/fitbivalent_new_long",well_idx,".png",sep="")
    ggsave(filename = data_file_name, plot=plot_bivalent_fit[[well_idx]])
  }
}
dev.off()

# want to order the output alphabetically by ligand
# sort_idx does that
sample_info$PlotSort <- translate_rows_cols_for_sort(sample_info$Row, 
                                                     sample_info$Column,
                                                     rows_dict,
                                                     sort_order="pairblock")


sample_info$WellIdx <- 1:nwells
sample_info %>% arrange(PlotSort, Column, Block) %>% select(WellIdx) -> sort_idx

sort_idx <- as_vector(sort_idx)

# bivalent
pages_list <- map(sort_idx, plot_bivalent_fitting,
                  fits_list, plot_bivalent_fit, sample_info$NumInclConc, rc_list, plot_list_before_baseline)

pdf(file = "paper_codes/figures/bivalent_new_long/Antigen2_out_bivalent_new.pdf")
for (well_idx in 1:nwells){
  print(well_idx)
  if (!is.null(pages_list[[well_idx]]$result)){
    grid.arrange(pages_list[[well_idx]]$result)
  }
}

for (well_idx in 1:nwells){
  plot(print_output(well_idx, pages_list, plot_list_out, sample_info))
}
dev.off()



