library(readr)
library(tidyverse)
library(readxl)
library(minpack.lm)
library(zoo)
library(gridExtra)
library(grid)
library(kableExtra)
library(deSolve)

setwd("~/Google Drive/My Drive/R/HIV_new_data/code")

load(file='monovalent_short.Rdata')
nwells <- length(fits_list)
monoShort_MSE <- NULL

for (well_idx in 1:nwells){
  num_data <- length(fits_list[[well_idx]]$result$R0$fvec)
  # print(num_data)
  monoShort_MSE <- rbind(monoShort_MSE, tail(fits_list[[well_idx]]$result$R0$rsstrace, 1)/num_data)
  # print(monoShort_MSE)
}

load(file='monovalent_long.Rdata')
nwells <- length(monovalent_fits_list)
monoLong_MSE <- NULL
for (well_idx in 1:nwells){
  num_data <- length(monovalent_fits_list[[well_idx]]$R0$fvec)
  monoLong_MSE <- rbind(monoLong_MSE, tail(monovalent_fits_list[[well_idx]]$R0$rsstrace, 1)/num_data)
}

load(file='bulkshift_short.Rdata')
nwells <- length(fits_list)
bulkshiftShort_MSE <- NULL
for (well_idx in 1:nwells){
  num_data <- length(fits_list[[well_idx]]$result$R0$fvec)
  bulkshiftShort_MSE <- rbind(bulkshiftShort_MSE, tail(fits_list[[well_idx]]$result$R0$rsstrace, 1)/num_data)
}

load(file='bulkshift_long.Rdata')
nwells <- length(fits_list)
bulkshiftLong_MSE <- NULL
for (well_idx in 1:nwells){
  num_data <- length(fits_list[[well_idx]]$result$R0$fvec)
  bulkshiftLong_MSE <- rbind(bulkshiftLong_MSE, tail(fits_list[[well_idx]]$result$R0$rsstrace, 1)/num_data)
}

load(file='bivalent_oldResults_short.Rdata')
nwells <- length(bivalent_fits_list)
bivalentShort_MSE <- NULL
for (well_idx in 1:nwells){
  num_data <- length(bivalent_fits_list[[well_idx]]$R0$fvec)
  bivalentShort_MSE <- rbind(bivalentShort_MSE, tail(bivalent_fits_list[[well_idx]]$R0$rsstrace, 1)/num_data)
}

load(file='bivalent_oldResults_long.Rdata')
nwells <- length(bivalent_fits_list)
bivalentLong_MSE <- NULL
for (well_idx in 1:nwells){
  num_data <- length(bivalent_fits_list[[well_idx]]$R0$fvec)
  bivalentLong_MSE <- rbind(bivalentLong_MSE, tail(bivalent_fits_list[[well_idx]]$R0$rsstrace, 1)/num_data)
}

load(file='bivalent_newResults_short.Rdata')
nwells <- length(bivalent_fits_list)
bivalentShort_MSE_new <- NULL
for (well_idx in 1:nwells){
  num_data <- length(bivalent_fits_list[[well_idx]]$R0$fvec)
  bivalentShort_MSE_new <- rbind(bivalentShort_MSE_new, tail(bivalent_fits_list[[well_idx]]$R0$rsstrace, 1)/num_data)
}

load(file='bivalent_newResults_long.Rdata')
nwells <- length(bivalent_fits_list)
bivalentLong_MSE_new <- NULL
for (well_idx in 1:nwells){
  num_data <- length(bivalent_fits_list[[well_idx]]$R0$fvec)
  bivalentLong_MSE_new <- rbind(bivalentLong_MSE_new, tail(bivalent_fits_list[[well_idx]]$R0$rsstrace, 1)/num_data)
}


#### Short
short_comparison <- NULL

MSE_short <- rbind(monoShort_MSE, bulkshiftShort_MSE, bivalentShort_MSE)
mono_name <- rep("1:1 Model", length(monoShort_MSE))
bulk_name <- rep("1:1 Model with Bulk Shift", length(bulkshiftShort_MSE))
biva_name <- rep("BA Model-1", length(bivalentShort_MSE))
biva_name_new <- rep("BA Model-2", length(bivalentShort_MSE))
monoShort_df <- bind_cols('MSE' = monoShort_MSE, 'Model' = mono_name)
bulkShort_df <- bind_cols('MSE' = bulkshiftShort_MSE, 'Model' = bulk_name)
bivaShort_df <- bind_cols('MSE' = bivalentShort_MSE, 'Model' = biva_name)
bivaShort_new_df <- bind_cols('MSE' = bivalentShort_MSE_new, 'Model' = biva_name_new)

short_df <- rbind(monoShort_df, bulkShort_df, bivaShort_df, bivaShort_new_df)

short_df$Model <- as.factor(short_df$Model)

title <- paste("Mean Square Error Comparison for Nominal Length of Dissociation")

p <- ggplot(short_df, aes(x=Model, y=MSE, color=Model)) +
  ggtitle(title) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #scale_fill_brewer(palette="Blues") + theme_classic()
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  theme(legend.position="none", text = element_text(size=20),
        axis.title.x=element_blank()) +
  ylim(2, 16)

p
#### Long
long_comparison <- NULL

MSE_long <- rbind(monoLong_MSE, bulkshiftLong_MSE, bivalentLong_MSE)
mono_name <- rep("1:1 Model", length(monoLong_MSE))
bulk_name <- rep("1:1 Model with Bulk Shift", length(bulkshiftLong_MSE))
biva_name <- rep("BA Model-1", length(bivalentLong_MSE))
biva_name_new <- rep("BA Model-2", length(bivalentLong_MSE))
monoLong_df <- bind_cols('MSE' = monoLong_MSE, 'Model' = mono_name)
bulkLong_df <- bind_cols('MSE' = bulkshiftLong_MSE, 'Model' = bulk_name)
bivaLong_df <- bind_cols('MSE' = bivalentLong_MSE, 'Model' = biva_name)
bivaLong_df_new <- bind_cols('MSE' = bivalentLong_MSE_new, 'Model' = biva_name_new)

long_df <- rbind(monoLong_df, bulkLong_df, bivaLong_df, bivaLong_df_new)

long_df$Model <- as.factor(long_df$Model)

title <- paste("Mean Square Error Comparison for Extended Length of Dissociation")

p <- ggplot(long_df, aes(x=Model, y=MSE, color=Model)) + 
  ggtitle(title) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  # scale_fill_brewer(palette="Blues") + theme_classic()
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  theme(legend.position="none", text = element_text(size=20),
        axis.title.x=element_blank()) +
  ylim(2, 16)

p
# 
# #####
# biva_short_name <- rep("Bivalent (nominal)", length(bivalentShort_MSE))
# biva_short_df <- bind_cols('MSE' = bivalentShort_MSE, 'Model' = biva_short_name)
# biva_long_name <- rep("Bivalent (extended)", length(bivalentLong_MSE))
# biva_long_df <- bind_cols('MSE' = bivalentLong_MSE, 'Model' = biva_long_name)
# 
# long_df <- rbind(mono_df, biva_short_df, biva_long_df)
# 
# p <- ggplot(long_df, aes(x=Model, y=MSE, color=Model)) + 
#   geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
#   #scale_fill_brewer(palette="Blues") + theme_classic()
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   ylim(2, 16)
# p
# 
# #geom_boxplot()
# 
# # load(file='bivalent.Rdata')
# # load(file='bulkshift.Rdata')
