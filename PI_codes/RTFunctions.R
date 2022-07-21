#' Automatically determine dissociation window
#'
#' @param well_idx The corresponding well in the extended sample sheet
#' @param sample_info A tibble. The expanded sample sheet
#' @param x_vals A tibble. Time columns from the carterra output for concentrations chosen for fittting
#' @param y_vals A tibble RU values from the carterra output for concentrations chosen for fittting
#' @param incl_concentrations_ligand A vector. The concentrations chosen for inclusion in fitting.
#' @param max_RU_tol A number. The maximum RU value to be included for determining the dissocation window
#' @param min_RU_tol A number. The minimum RU value to be included for determining the dissocation window
#' @return A number. The end dissociation time required to capture the decay in all of the selected concentrations for the given well
#' @examples
#' find_dissociation_window(well_idx, sample_info, x_vals, y_vals, incl_concentrations_ligand, max_RU_tol, min_RU_tol)

find_dissociation_window <- function(well_idx, sample_info, x_vals, y_vals, 
                                     incl_concentrations_ligand, max_RU_tol, min_RU_tol){
  
  if (sample_info[well_idx,]$`Automate Dissoc. Window` != "Y")
    return(NULL)
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  
  start_time <- baseline + baseline_start + association
  dissociation <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1
  
  end_dissoc_list <- NULL
  df_all <- NULL
  max_idx <- 0
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    RU <- y_vals[, i]
    df <- suppressMessages(bind_cols(Time, RU))
    names(df) <- c("Time", "RU")
    df %>% filter(Time > (start_time - 50)) -> df
    
    # Do not base on low information concentrations
    if (mean(df$RU, na.rm = TRUE) < min_RU_tol | mean(df$RU, na.rm = TRUE) > max_RU_tol)
      next
    df_zoo <- as.zoo(df)
    #  max_idx <- max_idx + 1
    as_tibble(
      rollapply(
        data = df_zoo, FUN =  function(x){
          x_df <- as_tibble(x)
          if (all(is.na(x_df$RU)| is.nan(x_df$RU)))
            return(c(NA,NA))
          else 
            return(coef(lm(RU ~ Time, singular.ok = TRUE,
                           data = x_df)))}, by.column = FALSE,
        
        width = 100)) -> df_out
    names(df_out) <- c("Intercept", "Slope")
    n_vals <- dim(df_out)[1]
    
    #  df_out %>% mutate(x = rep(max_idx, n_vals)) -> df_out
    df_out %>% mutate(RollIndex = 1:n_vals) -> df_out
    max_slope <- max(abs(df_out$Slope))
    target_slope <- .01*max_slope
    window_idx <- which(abs(df_out$Slope) < target_slope)[1]
    
    if (is.na(window_idx)){
      # Once this happens, we are using the entire time series for all concentrations
      end_dissoc_list <- c((df$Time)[length(df$Time)], end_dissoc_list)
      next
    } else
      end_dissoc <- df$Time[window_idx]
    
    end_dissoc_list <- c(end_dissoc, end_dissoc_list)
    
  }
  # Now we have a candidate end of dissoc for each concentration. Thia should be done for selected concentrations
  # Overall window is smallest window that accomodates all the concentrations
  return(max(as.numeric(end_dissoc_list), na.rm = TRUE))  
}

translate_rows_for_sort <- function(x){
  ifelse (x == "A", 1, 
          ifelse(x == "E", 2, 
                 ifelse(x == "B", 3,
                        ifelse(x == "F", 4, 
                               ifelse(x == "C", 5, 
                                      ifelse(x == "G", 6, 
                                             ifelse(x == "D", 7, 
                                                    ifelse(x == "H", 8, 0))))))))
  
}

first_conc_indices <- function(well_idx, num_conc_ligand){
  
  #Compute the correction to baseline. Usually for regenerative case.
  
  # find displacement from last ligand
  if (well_idx == 1){
    first_conc_idx <- 1
  } else 
    first_conc_idx <- sum(num_conc_ligand[1:(well_idx - 1)]) + 1 #number of time series up to this well
  first_conc_idx
}

#' @return A list containing:\tabular{ll}{
#'    \code{pars} \tab A numeric vector of parameter estimates \cr
#'    \tab \cr
#'    \code{std.errs} \tab A numeric vector of standard errors on parameters \cr
#'    \tab \cr
#'    \code{cov.mat} \tab Parameter covariance matrix (excluding mean) \cr
#' }
get_baseline_indices <- function(well_idx, sample_info, x_vals, y_vals){
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  baseline_avg_list <- NULL
  
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    RU <- y_vals[, i]
    df <- suppressMessages(bind_cols("Time" = Time, "RU" = RU))
    colnames(df) <- c("Time", "RU")
    df %>% filter(Time > baseline_start & Time < baseline+ baseline_start)  %>% .$RU -> base_meas
    baseline_avg <- mean(base_meas, na.rm = TRUE) 
    baseline_avg_list <- c(baseline_avg_list, baseline_avg)
  }
  min_baseline <- min(baseline_avg_list)
  baseline_idx <- start_idx + which(baseline_avg_list == min_baseline) - 1
  
  #is highest baseline negative?
  if (baseline_avg < 0)
    baseline_neg <- TRUE else
      baseline_neg <- FALSE
  
  list(BaselineIdx = baseline_idx, MinBaseline = min_baseline, BaselineNegative = baseline_neg)
}

create_dataframe_with_conc <- function(begin_conc_idx, end_conc_idx, x_vals, y_vals,
                                       numerical_concentrations,
                                       n_time_points){
  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)
  
  Time <- x_vals[, begin_conc_idx:end_conc_idx] %>%
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select(value)
  RU <- y_vals[, begin_conc_idx:end_conc_idx] %>%
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select("value")
  
  map_dfr(.x = tibble(numerical_concentrations), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(numerical_concentrations) -> numerical_concentrations
  Concentrations <- numerical_concentrations
  df <- suppressMessages(bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))
  colnames(df) <- c("Time", "RU", "Concentration")
  df
}

#This function is under development
bulkshift_correction <- function(well_idx, x_vals, y_vals, sample_info,
                                 all_concentrations_ligand){
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  association <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- all_concentrations_ligand[well_idx]
  
}

baseline_correction <- function(well_idx, x_vals, y_vals, sample_info){
  
  negative_baseline <- sample_info[well_idx, ]$BaselineNegative
  baseline_average <- sample_info[well_idx, ]$BaselineAverage
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  
  start_idx <- sample_info[well_idx, ]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc  #number of concentrations for this well
  end_idx <- start_idx + num_conc - 1
  
  
  #Need to test if baseline of highest conc is < 0
  
  if (sample_info[well_idx,]$Regen. == "N" & !negative_baseline){
    # add baseline average to all timepoints
    
    y_vals[, start_idx:end_idx] <-
      y_vals[, start_idx:end_idx] - baseline_average
  } else 
  {
    # correct all to mean of zero
    for (i in 1:num_conc){
      # compute baseline average for each concentration
      # subtract from baseline average from all RU vals for that concentration
      
      Time <- x_vals[, start_idx + (i-1)]
      RU <- y_vals[, start_idx + (i-1)]
      
      df <- suppressMessages(bind_cols("Time" = Time, "RU" = RU))
      colnames(df) <- c("Time", "RU")
      
      df %>% filter(Time > baseline_start & Time < baseline + baseline_start) %>% .$RU -> base_meas # select for defined baseline time period
      # This command is split up because mean(.$RU) would not parse properly
      base_corr <- mean(base_meas, na.rm = TRUE)
      y_vals[, start_idx + (i-1)] <- y_vals[, start_idx + (i-1)] - base_corr
    }
  }
  y_vals[, start_idx:end_idx]
  
}
get_response_curve <- function(well_idx, sample_info, x_vals, y_vals,
                               all_concentrations_values,
                               incl_concentrations_values, n_time_points){
  
  start_incl_idx <- sample_info[well_idx,]$FirstInclConcIdx
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1
  
  num_incl_conc <- sample_info[well_idx,]$NumInclConc
  end_incl_idx <- start_incl_idx + num_incl_conc - 1
  incl_conc_values <- incl_concentrations_values[start_incl_idx:end_incl_idx]
  
  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  
  df <- create_dataframe_with_conc(start_idx, end_idx, x_vals, y_vals,
                                   all_concentrations_values[start_idx:end_idx],
                                   n_time_points)
  
  df %>% filter((Time >= baseline + baseline_start + association - 10)
                & (Time <= baseline + baseline_start + association - 5)) %>%
    group_by(Concentration) %>%
    summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC
  df_RC %>% mutate(Included =
                     as_factor(ifelse(Concentration %in% incl_conc_values,
                                      "Yes", "No"))) -> df_RC
  ggplot(df_RC, aes(x = Concentration,
                    y = AverageRU)) +
    geom_point(aes(color = Included)) +
    geom_line() +
    scale_x_log10() +
    ggtitle(ligand_desc)
}



get_best_window <- function(well_idx, sample_info, x_vals, y_vals, 
                            num_conc, concentrations, displacement_per_ligand){
  
  association <- sample_info[well_idx,]$Association
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  end_assoc_time <- association + baseline + baseline_start
  start_idx <- displacement_per_ligand + 1
  end_idx <- start_idx + num_conc - 1
  end_assoc_resp <- NULL
  
  association_end <- baseline + baseline_start + association
  
  df <- create_dataframe_with_conc(start_idx, end_idx, x_vals, y_vals,
                                   concentrations,
                                   n_time_points)
  
  #this fails for well 61 because there are no times that meet the criteria
  df %>% filter((Time >= baseline + baseline_start + association - 10) 
                & (Time <= baseline + baseline_start + association - 5)) %>% 
    group_by(Concentration) %>% 
    summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC
  
  # check to see if any results for df_RC
  
  if (dim(df_RC)[1] == 0)
    return(NULL)
  
  # record the differences between consecutive responses
  sum_diff <- NULL
  for (i in 1:(num_conc - 1)){
    sum_diff <- suppressMessages(bind_cols(sum_diff, df_RC[i+1,]$AverageRU - df_RC[i,]$AverageRU))
  }
  # add sums for each 5 cycle window.
  cum_sum <- rollapply(as_vector(flatten(sum_diff)), 4, FUN = sum)
  start_conc_idx <- which(cum_sum == max(cum_sum))
  
  # check start and end slopes, may be better to fit 4 instead of five
  
  slopes <- as_vector(sum_diff[start_conc_idx:(start_conc_idx+3)])
  
  remove_concentration <- ifelse(slopes < 0.2*mean(slopes), 1, 0)
  
  # return 5 best consecutive concentrations 
  if(sum(remove_concentration) == 0)
    return(concentrations[start_conc_idx:(start_conc_idx + 4)])
  
  # if some slopes are less than 20% of the mean, remove one concentration 
  # low end or high end, depending on which slope is smaller
  
  # remove_concentration has at least one '1' value
  # if both first and last are tagged, we remove the smallest
  # if only one is tagged, it will still be the smallest
  if (slopes[1] < slopes[4])
    return(concentrations[(start_conc_idx+1):(start_conc_idx + 4)])
  else
    return(concentrations[(start_conc_idx):(start_conc_idx + 3)])
  
}

plot_sensorgrams <- function(well_idx, 
                             sample_info,
                             x_vals,
                             y_vals, 
                             incl_conc_values, 
                             all_concentrations_values,
                             n_time_points,
                             all_concentrations = FALSE){
  if (!all_concentrations){
    start_idx <- sample_info[well_idx,]$FirstInclConcIdx
    num_conc <- sample_info[well_idx,]$NumInclConc
  } else
  {
    start_idx <- sample_info[well_idx,]$FirstConcIdx
    num_conc <- sample_info[well_idx,]$NumConc
    incl_conc_values <- all_concentrations_values
  }
  end_idx <- start_idx + num_conc - 1
  
  incl_conc_values <- incl_conc_values[start_idx:end_idx]
  
  ligand_desc <- sample_info[well_idx,]$Ligand
  
  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)
  
  
  Time <- x_vals[, start_idx:end_idx] %>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select(value)
  RU <- y_vals[, start_idx:end_idx]%>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select("value")
  
  
  map_dfr(.x = tibble(incl_conc_values), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(incl_conc_values) -> incl_conc_values
  
  Concentrations <- incl_conc_values 
  
  df <- suppressMessages(bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  df$Concentration <- as_factor(formatC(df$Concentration, format = "e",digits = 2))
  
  
  sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row", sample_info[well_idx,]$Row,
                     "Column", sample_info[well_idx,]$Column)
  
  ggplot(df, aes(x = Time, y = RU, color = Concentration)) + geom_point(size = 0.5) +
    ggtitle(ligand_desc, subtitle = sub_title) 
}

summary_fit_with_constraints <- function(fit_object){
  hessian <- fit_object$hessian
  pars <- fit_object$par
  info <- fit_object$info
  
  n <- nrow(hessian)
  test_zeroes <- apply(hessian, 1 , function(x) sum(x==0))
  nonsingular_rows <- which(test_zeroes != n)
  if (length(nonsingular_rows) == n & info != 5)
    return(summary(fit_object)$coefficients)
  hessian <- hessian[nonsingular_rows, nonsingular_rows]
  
  std_err_full <- rep(NA, n)
  
  # get table directly. Code is pulled from summary.nls.lm
  
  if (info != 5) {            # when info is 5, that means the iterations maxed out and fit is not valid
    ibb <- chol(hessian)
    ih <- chol2inv(ibb)
    p <- length(pars)
    rdf <- length(fit_object$fvec) - p
    resvar <- deviance(fit_object)/rdf
    se <- sqrt(diag(ih) * resvar)
  }
  
  std_err_full[nonsingular_rows] <- se
  
  df <- data.frame(Estimate = pars, "Std. Error" = std_err_full)
  colnames(df) <- c("Estimate", "Std. Error")
  df
}


fit_kd <- function(pars, df, incl_concentrations, num_conc, kd, t0 = t0){
  #pars ("Rmax" one for each concentration,"ka", "tstart" one for each concentration)
  
  R0 <- pars[1:num_conc]
  kd <- pars[(num_conc+1)]
  
  err_assoc <- NULL
  err_dissoc <- NULL
  
  for (i in 1:num_conc){
    
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    
    dissoc_formula <- R0[i]*exp(-kd*(df_dissoc$Time - t0))
    
    err_dissoc <-   c(err_dissoc, df_dissoc$RU - dissoc_formula)
    
  }
  err_dissoc
  
}

########################################################################################
# ODE Models
########################################################################################
heterogenousLigand_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dB1  <- -ka1*Am*B1 + kd1*AB1
    dB2  <- -ka2*Am*B2 + kd2*AB2
    dAB1 <-  ka1*Am*B1 - kd1*AB1
    dAB2 <-  ka2*Am*B2 - kd2*AB2
    
    #return the rate of change
    list(c(dB1, dB2, dAB1, dAB2))
  })
}

bivalentAnalyte_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dL  <- -(2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2)
    dX1 <-  (2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2)
    dX2 <-  ka2*X1*L - 2*kd2*X2
    
    #return the rate of change
    list(c(dL,dX1,dX2))
  })
}

monovalent_model <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dL <- -ka1*Am*L + kd1*X1
    dX1 <- ka1*Am*L - kd1*X1
    
    # return the rate of change
    list(c(dL, dX1))
  })
}

heterogenousAnalyte_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dB   <- -(ka1*Am1*mw1*B - kd1*A1B)/mw1*n1 - (ka2*Am2*mw2*B - kd2*A2B)/mw2*n2
    dA1B <- ka1*Am1*mw1*B - kd1*A1B
    dA2B <- ka2*Am2*mw2*B - kd2*A2B
    
    #return the rate of change
    list(c(dB, dA1B, dA2B))
  })
}

twoState_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dB   <- -(ka1*Am*B - kd1*AB)
    dAB  <- ka1*Am*B - kd1*AB - (ka2*AB - kd2*ABx)
    dABx <- ka2*AB - kd2*ABx
    
    #return the rate of change
    list(c(dB,dAB,dABx))
  })
}

biEpitopicLigandHeterogenousAnalyte_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dL <- - ka1*A1m*L + kd1*A1L - ka2*A2m*L + kd2*A2L
    dA1L <- ka1*A1m*L - kd1*A1L - ka2*A1L*A2m + kd2*A1LA2
    dA2L <- ka2*A2m*L - kd2*A2L - ka1*A1m*A2L + kd1*A1LA2
    dA1LA2 <- ka1*A1m*A2L - kd1*A1LA2 + ka2*A1L*A2m - kd2*A1LA2
    
    #return the rate of change
    # print(list(c(t, L, A1L, A2L, A1LA2)))
    # print(dL + dA1L + dA2L + dA1LA2)
    list(c(dL, dA1L, dA2L, dA1LA2))
  })
}


objective_function <- function(pars, df, incl_concentrations, num_conc, param_names, 
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind){
  
  #pars ("Rmax" one for each concentration,"ka", "tstart" one for each concentration)
  
  err <- NULL
  RU_data <- NULL
  
  for (i in 1:num_conc){
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(AssocIndicator == 1) -> df_assoc
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    df_dissoc <- df_dissoc[seq(1, dim(df_dissoc)[1], 1),]
    
    RU_temp <- c(df_assoc$RU, df_dissoc$RU)
    RU_data <- c(RU_data, RU_temp)
  }
  tryCatch(
    expr = {
      # Your code...
      # goes here...
      # ...
      sim_results <- run_model(pars, df, num_conc, incl_concentrations, param_names, 
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind)
      solved_RU <- sim_results$RU
      err <- na.omit(solved_RU - RU_data)
    },
    error = function(e){
      # (Optional)
      # Do this if an error is caught...
      err <- na.omit(rep(-1e3, length(RU_data)) - RU_data)
    }
  )
  err
}

fit_association_dissociation <- function(well_idx, sample_info, x_vals, y_vals,
                                         incl_concentrations_values,
                                         min_allowed_kd = 10^(-5),
                                         max_iterations = 500,
                                         ptol = 10^(-10),
                                         ftol = 10^(-10)){
  
  # this function will fit all concentrations for one well
  use_bulkShift <- sample_info[well_idx,]$Bulkshift
  use_regeneration <- sample_info[well_idx,]$Regen.
  use_globalRmax <- sample_info[well_idx,]$`Global Rmax`
  use_secondRebind <- sample_info[well_idx,]$`Sec. Rebinding`
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  
  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1
  
  association <- sample_info[well_idx,]$Association
  dissociation <- sample_info[well_idx,]$Dissociation
  
  assoc_start <- baseline + baseline_start
  assoc_end <- assoc_start + association
  
  dissoc_start <- assoc_end
  dissoc_end <- assoc_end + dissociation
  
  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)
  
  Time <- x_vals[, start_idx:end_idx] %>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select(value)
  RU <- y_vals[, start_idx:end_idx]%>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select("value")
  
  incl_concentrations <-
    incl_concentrations_values[start_idx:end_idx]
  
  map_dfr(.x = tibble(incl_concentrations), 
          .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(incl_concentrations) -> incl_concentrations_rep
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = incl_concentrations_rep)
  colnames(df) <- c("Time", "RU", "Concentration") #force correct names - dplyr is doing weird things
  df <- na.omit(df)
  #df$Time <- df$Time - baseline - baseline_start
  
  #do both dissociation and association
  df %>% mutate(AssocIndicator = 
                  ifelse((Time >= assoc_start & Time < assoc_end), 1, 0), 
                DissocIndicator = ifelse(Time > dissoc_start & Time < dissoc_end, 1, 0)) -> df
  #df %>% filter(DissocIndicator == 1) -> df
  
  df %>% group_by(Concentration) %>% summarise(max = max(RU, na.rm = TRUE)) -> Rmax_start
  df %>% filter(AssocIndicator == 1) -> df_assoc
  df_assoc %>% group_by(Concentration) %>% summarise(min = min(RU, na.rm = TRUE)) -> R0_start
  
  ###########################  ###########################  ###########################
  ###########################  ###########################  ###########################
  
  model <- sample_info[well_idx,]$Model
  
  if (use_globalRmax == "N"){
    num_Rmax <- num_conc
  }
  else{
    num_Rmax <- 1
  }
  
  if (use_regeneration == "N"){
    num_R0 <- num_conc
    init_R0s <- R0_start$min
  }
  else{
    num_R0 <- 0
    init_R0s <- NULL
  }
  
  if (model == 'heterogenousLigand'){
    lower <- c(0, 0, 0, 0, rep(0, 2*num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", 2*num_Rmax), rep("R0", num_R0))
  }
  else if (model == 'bivalentAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep("R0", num_R0))
  }
  else if (model == 'monovalent'){
    lower <- c(0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "kd1", rep("Rmax", num_Rmax), rep("R0", num_R0))
  }
  else if (model == 'heterogenousAnalyte'){
    lower <- c(0, 0, 0, 0, 0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "ka2", "kd1", "kd2", "mw1", "mw2", "n1", "n2", rep("Rmax", num_Rmax), rep("R0", num_R0))
  }
  else if (model == 'twoState'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep("R0", num_R0))
  }
  else if (model == 'biEpitopicLigandHeterogenousAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep("R0", num_R0))
  }
  
  init_params_list <- NULL
  
  best_error <- 1e64
  t0_start <- dissoc_start
  estimated_params_list <- NULL
  estimated_fval <- NULL
  
  time_start <- Sys.time()
  
  num_run <- sample_info[well_idx,]$`Num Run`
  
  for (i in 1:2){
    
    if (model == 'heterogenousLigand'){
      lower_powers <- c(1, -7, -7, -7)
      upper_powers <- c(7, -1, -1, -1)
      init_ks <- runif(4, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, 2*num_Rmax)
      
      # ka1, ka2, kd1, kd2, Rmax_1-10, RU0_1-5
      init_params <- c(init_ks, init_Rmaxs, init_R0s)
    }
    else if (model == 'bivalentAnalyte'){
      lower_powers <- c(1, -7, -7, -7)
      upper_powers <- c(7, -1, -1, -1)
      init_ks <- runif(4, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, num_Rmax)
      
      # ka1, ka2, kd1, kd2, Rmax_1-5, RU0_1-5
      init_params <- c(init_ks, init_Rmaxs, init_R0s)
    }
    else if (model == 'monovalent'){
      lower_powers <- c(1, -7)
      upper_powers <- c(7, -1)
      init_ks <- runif(2, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, num_Rmax)
      
      # ka1, kd1, Rmax_1-5, RU0_1-5
      init_params <- c(init_ks, init_Rmaxs, init_R0s)
    }
    else if (model == 'heterogenousAnalyte'){
      lower_powers <- c(1, -7, -7, -7)
      upper_powers <- c(7, -1, -1, -1)
      init_ks <- runif(4, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, num_Rmax)
      
      # ka1, ka2, kd1, kd2, mw1, mw2, n1, n2, Rmax_1-5, RU0_1-5
      init_params <- c(init_ks, 0.5, 0.5, 1, 1, init_Rmaxs, init_R0s)
    }
    else if (model == 'twoState'){
      lower_powers <- c(1, -7, -7, -7)
      upper_powers <- c(7, -1, -1, -1)
      init_ks <- runif(4, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, num_Rmax)
      
      # ka1, ka2, kd1, kd2, Rmax_1-5, RU0_1-5
      init_params <- c(init_ks, init_Rmaxs, init_R0s)
    }
    else if (model == 'biEpitopicLigandHeterogenousAnalyte'){
      lower_powers <- c(1, -7, -7, -7)
      upper_powers <- c(7, -1, -1, -1)
      init_ks <- runif(4, min = 0, max = 1)*(upper_powers - lower_powers) + lower_powers
      init_ks <- 10^init_ks
      
      # 
      init_Rmaxs <- max(Rmax_start$max)*rep(1, num_Rmax)
      
      # ka1, ka2, kd1, kd2, Rmax_1-5, RU0_1-5
      init_params <- c(init_ks, init_Rmaxs, init_R0s)
    }
    
    
    res_R0 <- NULL
    error_val <- 2e64
    
    print("Initial Parameters:")
    print(init_params)
    
    tryCatch(
      expr = {
        # Your code...
        # goes here...
        # ...
        res_R0 <- nls.lm(init_params, fn = objective_function, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc, 
                         param_names = param_names,
                         use_regeneration = use_regeneration, model = model, use_globalRmax = use_globalRmax,
                         use_secondRebind = use_secondRebind,
                         control = nls.lm.control(maxiter = 1000),
                         lower = lower)
        
        estimated_fval <- rbind(estimated_fval, res_R0$deviance)
        estimated_params_list <- rbind(estimated_params_list, coefficients(res_R0))
        error_val <- res_R0$deviance
        
        # print("Estimated Parameters:")
        # print(coefficients(res_R0))
      },
      error = function(e){
        # (Optional)
        # Do this if an error is caught...
        res_R0$deviance <- 2e64
        estimated_fval <- rbind(estimated_fval, res_R0$deviance)
        estimated_params_list <- rbind(estimated_params_list, init_params)
        error_val <- res_R0$deviance
      }
    )
    
    
    if (error_val < best_error){
      # best_res_R0 <- NULL
      # best_error <- 0
      best_res_R0 <- res_R0
      best_error <- error_val
    }
  }
  
  # full_param_table <- cbind(init_params_list, estimated_params_list, estimated_fval)
  # full_param_table <- full_param_table %>% as_tibble() #&>& setNames()
  # write_csv(full_param_table, paste('~/Google Drive/My Drive/R/HIV_new_data/code/param_table_biexop/full_param',as.character(well_idx),'.csv',sep=''))
  # 
  # kd_param_table <- as.matrix(full_param_table[,c(1:4,15:18,29)])
  # kd_param_table <- kd_param_table %>% as_tibble() %>% setNames(c("ka1_0", "ka2_0", "kd1_0", "kd2_0", "ka1", "ka2", "kd1", "kd2", "SSE"))
  # write_csv(kd_param_table, paste('~/Google Drive/My Drive/R/HIV_new_data/code/param_table_biexop/kd_param',
  #                                 as.character(well_idx),'.csv',
  #                                 sep=''))
  
  time_end <- Sys.time()
  elapsed_time <- time_end - time_start
  print("Time:")
  print(elapsed_time)
  
  pars <- coefficients(best_res_R0)
  print("Best Params:")
  print(pars)
  
  fit_outcomes <- run_model(pars, df, num_conc, incl_concentrations, param_names, 
                            use_regeneration, 
                            model, 
                            use_globalRmax,
                            use_secondRebind)
  
  list("R0" = best_res_R0, "FitOutcomes" = fit_outcomes)
}

run_model <- function(pars, df, num_conc, incl_concentrations, param_names, 
                      use_regeneration, 
                      model, 
                      use_globalRmax,
                      use_secondRebind){
  
  full_output_RU <- NULL
  
  for (i in 1:num_conc){
    
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(AssocIndicator == 1) -> df_assoc
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    df_dissoc <- df_dissoc[seq(1, dim(df_dissoc)[1], 1),]
    
    
    if (model == 'heterogenousLigand'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 4
      
      if (use_globalRmax == "N"){
        state <- c(B1  = Rmaxs[i],
                   B2  = Rmaxs[num_conc+i],
                   AB1 = 0,
                   AB2 = 0)
      }
      else{
        state <- c(B1  = Rmaxs[1],
                   B2  = Rmaxs[2],
                   AB1 = 0,
                   AB2 = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = Concentration[i])
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = heterogenousLigand_model, parms = parameters)
      df_assoc$RU <- out[,4] + out[,5] + RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc), df_dissoc$Time])
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = 0)
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4],
                 out[length(t_asc),5])
      
      out <- ode(y = state, times = t_dis, func = heterogenousLigand_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),4] + out[2:length(t_dis),5] + RI
      
    }
    else if (model == 'bivalentAnalyte'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 4
      
      if (use_globalRmax == "N"){
        state <- c(L  = Rmaxs[i],
                   X1 = 0,
                   X2 = 0)
      }
      else{
        state <- c(L  = Rmaxs[1],
                   X1 = 0,
                   X2 = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = Concentration[i])
      
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], df_dissoc$Time)
      
      if (use_secondRebind == "Y"){
        parameters <- c(ka1 = ka1,
                        kd1 = kd1,
                        ka2 = ka2,
                        kd2 = kd2,
                        Am = 0)
      }
      else{
        parameters <- c(ka1 = 0,
                        kd1 = kd1,
                        ka2 = 0,
                        kd2 = kd2,
                        Am = 0)
      }
      
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = bivalentAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),3] + out[2:length(t_dis),4] + RI
      
    }
    else if (model == 'monovalent'){
      
      ka1 <- pars[param_names == "ka1"]
      kd1 <- pars[param_names == "kd1"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 2
      
      if (use_globalRmax == "N"){
        state <- c(L  = Rmaxs[i],
                   X1 = 0)
      }
      else{
        state <- c(L  = Rmaxs[1],
                   X1 = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      Am = Concentration[i])
      
      
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = monovalent_model, parms = parameters)
      df_assoc$RU <- out[,3] + RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)],df_dissoc$Time)
      
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      Am = 0)
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3])
      
      out <- ode(y = state, times = t_dis, func = monovalent_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),3] + RI
      
    }
    else if (model == 'heterogenousAnalyte'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      mw1 <- pars[param_names == "mw1"]
      mw2 <- pars[param_names == "mw2"]
      n1  <- pars[param_names == "n1"]
      n2  <- pars[param_names == "n2"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 8
      
      if (use_globalRmax == "N"){
        state <- c(B   = Rmaxs[i],
                   A1B = 0,
                   A2B = 0)
      }
      else{
        state <- c(B   = Rmaxs[1],
                   A1B = 0,
                   A2B = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      mw1 = mw1,
                      mw2 = mw2,
                      n1  = n1,
                      n2  = n2,
                      Am1 = Concentration[i],
                      Am2 = Concentration[i]) ############## Need to change this to reflect the info file
      
      
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = heterogenousAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] +  RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], df_dissoc$Time)
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      mw1 = mw1,
                      mw2 = mw2,
                      n1  = n1,
                      n2  = n2,
                      Am1 = 0, ############## Need to change this to reflect the info file
                      Am2 = 0) ############## Need to change this to reflect the info file
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = heterogenousAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),3] + out[2:length(t_dis),4] +  RI
      
    }
    else if (model == 'twoState'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 4
      
      if (use_globalRmax == "N"){
        state <- c(B   = Rmaxs[i],
                   AB  = 0,
                   ABx = 0)
      }
      else{
        state <- c(B   = Rmaxs[1],
                   AB  = 0,
                   ABx = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = Concentration[i])
      
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = twoState_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], df_dissoc$Time)
      
      if (use_secondRebind == "Y"){
        parameters <- c(ka1 = ka1,
                        kd1 = kd1,
                        ka2 = ka2,
                        kd2 = kd2,
                        Am = 0)
      }
      else{
        parameters <- c(ka1 = 0,
                        kd1 = kd1,
                        ka2 = 0,
                        kd2 = kd2,
                        Am = 0)
      }

      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = twoState_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),3] + out[2:length(t_dis),4] + RI
      
    }
    else if (model == 'biEpitopicLigandHeterogenousAnalyte'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      Rmaxs <- pars[param_names == "Rmax"]
      
      num_pars <- 4
      
      if (use_globalRmax == "N"){
        state <- c(L   = Rmaxs[i],
                   A1L = 0,
                   A2L = 0,
                   A1LA2 = 0)
      }
      else{
        state <- c(L   = Rmaxs[1],
                   A1L = 0,
                   A2L = 0,
                   A1LA2 = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      A1m = Concentration[i],
                      A2m = Concentration[i]) ############## Need to change this to reflect the info file
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = biEpitopicLigandHeterogenousAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + 2*out[,5] +  RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], df_dissoc$Time)
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      ka2 = 0,
                      kd2 = kd2,
                      A1m = Concentration[i],
                      A2m = Concentration[i])
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4],
                 out[length(t_asc),5])
      
      out <- ode(y = state, times = t_dis, func = biEpitopicLigandHeterogenousAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[2:length(t_dis),3] + out[2:length(t_dis),4] + 2*out[2:length(t_dis),5] +  RI
    }
    
    full_output_RU <- bind_rows(full_output_RU, df_assoc, df_dissoc)
  }
  full_output_RU %>% select(Time, RU, Concentration)
}

plot_sensorgrams_with_fits <- function(well_idx, sample_info, fits, x_vals, y_vals, 
                                       incl_conc_values, n_time_points){
  
  if (!is.null(fits[[well_idx]]$error))
    return(NULL)
  
  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1
  
  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  association <- sample_info[well_idx,]$Association
  dissociation <- sample_info[well_idx,]$Dissociation
  
  fit_df <- fits[[well_idx]]$FitOutcomes
  
  n_vals <- dim(x_vals)[2]
  names(x_vals) <- 1:n_vals
  names(y_vals) <- 1:n_vals
  
  
  Time <- x_vals[, start_idx:end_idx] %>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select(value)
  RU <- y_vals[, start_idx:end_idx]%>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select("value")
  
  numerical_concentration <- incl_conc_values[start_idx:end_idx]
  
  map_dfr(.x = tibble(numerical_concentration), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(numerical_concentration) -> numerical_concentration
  
  Concentrations <- numerical_concentration 
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  df %>% filter(Time > baseline & Time < baseline + association + dissociation) -> df
  #df %>% filter(AssocIndicator == 1 | DissocIndicator == 1) -> df
  df <- na.omit(df)
  #bind_cols(df,FittedRU = fit_RU) -> df
  
  colnames(df) <- c("Time", "RU", "Concentration", "FittedRU")
  
  df$Concentration <- as_factor(df$Concentration)
  
  sub_title <- paste("1:1 Langmuir Model Fitting with Nominal Length of Dissociation")
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
  #                    sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  dissociation_start <- sample_info$Baseline[[well_idx]] + sample_info$`Bsl Start`[[well_idx]] + sample_info$Association[[well_idx]]
  
  # plot <- ggplot(df, aes(x = Time, y = RU)) + geom_point(size = 0.09, aes(color = Concentration)) #+
  #   #geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1) +
  #   geom_vline(xintercept = dissociation_start, linetype="dashed",
  #              color = "black", size=1)
  
  ggplot(df, aes(x = Time, y = RU)) + geom_point(size = 0.09, aes(color = Concentration)) +
    ggtitle(ligand_desc, subtitle = sub_title) +
    geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1) +
    geom_vline(xintercept = dissociation_start, linetype="dashed", color = "black", size=0.5)
}



plot_bivalent_fitting <- function(well_idx, fits_list, plot_list_out, 
                                  num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  # Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x)) # For same Rmax
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  "Rmax", R0_label))) #R0_label "kd2"
  # par_names <- as_vector(flatten(c("ka1", "kd1", Rmax_label, R0_label))) #R0_label "kd2"
  
  pars <- coefficients(fits_list[[well_idx]]$R0)
  
  result_summary <- summary(fits_list[[well_idx]]$R0)
  summary_names <- colnames(result_summary$coefficients)
  result_summary$coefficients %>% as_tibble -> par_err_table
  
  colnames(par_err_table) <- summary_names
  par_err_table <- bind_cols(Names = par_names, par_err_table)
  
  par_err_table %>% filter(!str_detect(Names,"R_0")) -> par_err_table
  par_names <- par_err_table$Names
  par_err_table %>% select(Estimate, `Std. Error`)  %>%
    mutate(Estimate = format(Estimate,big.mark=",",decimal.mark=".")) %>%
    tableGrob(rows = par_names, theme = ttheme_minimal()) -> tb1
  
  residuals(fits_list[[well_idx]]$R0) -> RU_resid
  fits_list[[well_idx]]$FitOutcomes$Time -> Time_resid
  fits_list[[well_idx]]$FitOutcomes$RU -> RU
  fits_list[[well_idx]]$FitOutcomes$Concentration -> Concentration
  resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = RU, Concentration = as.factor(Concentration)), 
                       aes(x = RU, y = Residuals)) + geom_point(size = 0.01, aes(color = Concentration)) +
    ggtitle(label = "Residuals")
  
  plot_before_baseline <- plot_list_before_baseline[[well_idx]] + ggtitle("Raw Sensorgram")
  # grid.arrange(plot_list_out[[well_idx]])
  grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
               rc_list[[well_idx]], ncol=2)
}









combine_output <- function(well_idx, fits_list, plot_list_out, rc_list, sample_info){
  
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  
  if (sample_info[well_idx,]$`Global Rmax` == "Y")
    global_rmax <- TRUE else
      global_rmax <- FALSE
    
    if (sample_info[well_idx,]$`Bulkshift` == "Y")
      bulkshift <- TRUE else
        bulkshift <- FALSE
      
      num_conc <- sample_info[well_idx,]$NumInclConc
      
      if (global_rmax)
        Rmax_label <- "Rmax"
      else
        Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
      
      
      R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
      bulkshift_label <- map_dfr(tibble(1:num_conc), function(x) paste("Bulkshift", x))
      
      
      pars <- coefficients(fits_list[[well_idx]]$result$R0)
      
      if (bulkshift)
        par_names <- as_vector(flatten(c(Rmax_label, "ka", R0_label, "kd", bulkshift_label)))
      else
        par_names <- as_vector(flatten(c(Rmax_label, "ka", R0_label, "kd")))
      
      # result_summary <- summary(fits_list[[well_idx]]$result$R0)
      #R's built-in summary method doesn't play nicely when the some of the parameters hit their limiting values (the hessian is singular)
      # I've adapted the function to return NA's for std error when the limits are reached.
      
      result_summary <- summary_fit_with_constraints(fits_list[[well_idx]]$result$R0)
      #summary_fit_with_constraints returns the coefficients table from summary.nls.lm
      
      summary_names <- colnames(result_summary)
      result_summary %>% as_tibble -> par_err_table
      
      colnames(par_err_table) <- summary_names
      par_err_table <- suppressMessages(bind_cols(Names = par_names, par_err_table))
      
      par_err_table %>% filter(!str_detect(Names,"R_0")) -> par_err_table
      par_err_table %>% filter(!str_detect(Names,"Bulkshift")) -> par_err_table
      
      
      par_names <- par_err_table$Names
      
      par_err_table %>% 
        filter(Names == "ka" | Names == "kd") %>%
        select(Estimate, `Std. Error`)  %>%
        mutate(Estimate = format(signif(Estimate, 3),big.mark=",",decimal.mark=".", scientific = TRUE)) %>%
        mutate(`Std. Error` = format(signif(`Std. Error`, 3),big.mark=",",decimal.mark=".", scientific = TRUE)) -> kakd_out
      
      par_err_table %>% 
        filter(!(Names == "ka" | Names == "kd")) %>% 
        select(Estimate, `Std. Error`)  %>%
        mutate(Estimate = format(round(Estimate,2),big.mark=",",decimal.mark=".", scientific = FALSE))  %>%
        mutate(`Std. Error` = format(round(`Std. Error`, 2),big.mark=",",decimal.mark=".", scientific = FALSE)) -> rest_out
      
      bind_rows(rest_out, kakd_out) %>%
        tableGrob(rows = par_names, theme = ttheme_minimal()) -> tb1
      
      residuals(fits_list[[well_idx]]$result$R0) -> RU_resid
      fits_list[[well_idx]]$result$FitOutcomes$Time -> Time_resid
      fits_list[[well_idx]]$result$FitOutcomes$Concentration -> Concentration_resid
      
      
      resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = Time_resid, Concentration = as_factor(Concentration_resid)), 
                           aes(x = Time, y = Residuals, color = Concentration)) + geom_point(size = 0.01) +
        ggtitle(label = "Residuals")
      
      grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
                   rc_list[[well_idx]], ncol=2)
}

print_output <- function(well_idx, pages_list, plot_list_out, sample_info){
  
  if (is.null(pages_list[[well_idx]]$error) & !is.null(pages_list[[well_idx]]$result))
    return(arrangeGrob(pages_list[[well_idx]]$result))
  
  err_msg <- paste("The following well has an unrecoverable error:",
                   well_idx, "Block ", sample_info$Block, "Row", sample_info$Row)
  err_msg <- paste0(err_msg, sample_info$Column)
  
  if (is.null(plot_list_out[[well_idx]]))
    return(p1 = arrangeGrob(textGrob(err_msg), plot_list[[well_idx]]))
  else
    return(p1 = arrangeGrob(textGrob(err_msg), plot_list_out[[well_idx]]))
  
}     
get_response_curve <- function(well_idx, sample_info, x_vals, y_vals, 
                               all_concentrations_values,
                               incl_concentrations_values, 
                               n_time_points){
  
  start_incl_idx <- sample_info[well_idx,]$FirstInclConcIdx
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  num_incl_conc <- sample_info[well_idx,]$NumInclConc
  end_incl_idx <- start_incl_idx + num_incl_conc - 1
  end_idx <- start_idx + num_conc - 1
  
  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  
  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)
  
  
  Time <- x_vals[, start_idx:end_idx] %>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select(value)
  RU <- y_vals[, start_idx:end_idx]%>% 
    pivot_longer(cols = everything()) %>% arrange(as.numeric(name)) %>%
    select("value")
  
  
  numerical_concentration <- all_concentrations_values[start_idx:end_idx]
  numerical_concentration_incl <-
    incl_concentrations_values[start_incl_idx:end_incl_idx]
  
  map_dfr(.x = tibble(numerical_concentration), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(numerical_concentration) -> numerical_concentration
  
  Concentrations <- numerical_concentration 
  
  df <- suppressMessages(bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  
  df %>% filter((Time >= baseline + baseline_start + association - 10) 
                & (Time <= baseline + baseline_start + association - 5)) %>% 
    group_by(Concentration) %>% 
    summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC
  
  df_RC %>% mutate(Included = 
                     as_factor(ifelse(Concentration %in% numerical_concentration_incl,
                                      "Yes", "No"))) -> df_RC
  
  
  ggplot(df_RC, aes(x = Concentration, 
                    y = AverageRU)) + 
    geom_point(aes(color = Included)) +
    geom_line() +
    scale_x_log10() +
    ggtitle(ligand_desc)
}

get_csv <- function(well_idx, fits_list, sample_info){
  
  num_conc <- sample_info_fits[well_idx,]$NumInclConc
  
  if (sample_info[well_idx,]$`Global Rmax` == "Y")
    global_rmax <- TRUE else
      global_rmax <- FALSE
    
    if (!global_rmax){
      Rmax <- rep(NA, 5)
      Rmax_se <- rep(NA,5)
    }
    
    Bulkshift <- rep(NA, 5)
    Bulkshift_se <- rep(NA,5)
    R0 <- rep(NA,5)
    R0_se <- rep(NA,5)
    
    pars <- coefficients(fits_list[[well_idx]]$result$R0)
    
    if (sample_info[well_idx, ]$Bulkshift == "Y")
      bulkshift <- TRUE else
        bulkshift <- FALSE
    
    
    
    result_summary <- summary_fit_with_constraints(fits_list[[well_idx]]$result$R0)
    
    # first column of summary is the estimate. Second is standard error
    
    summary_names <- colnames(result_summary)
    result_summary %>% as_tibble %>% select(Estimate, `Std. Error`) -> par_err_table
    
    colnames(par_err_table) <- summary_names[1:2]
    
    if (global_rmax){
      Rmax <- par_err_table[1,]$Estimate
      Rmax_se <- par_err_table[1,]$`Std. Error`
      ka <- par_err_table[2,]$Estimate
      ka_se <- par_err_table[2,]$`Std. Error`
      R0[1:num_conc] <- par_err_table[3:(num_conc + 2), ]$Estimate
      R0_se[1:num_conc] <- par_err_table[3:(num_conc + 2), ]$`Std. Error`
      kd <- par_err_table[2*num_conc + 3,]$Estimate
      kd_se <- par_err_table[2*num_conc + 3,]$`Std. Error`
      if (bulkshift){
        Bulkshift[1:num_conc] <- par_err_table[(2*num_conc + 4):(3*num_conc+3), ]$Estimate
        Bulkshift_se[1:num_conc] <- par_err_table[(2*num_conc + 4):(3*num_conc+3), ]$`Std. Error`
      }
      
      
    } else {
      Rmax[1:num_conc] <- par_err_table[1:num_conc,]$Estimate
      Rmax_se[1:num_conc] <- par_err_table[1:num_conc,]$`Std. Error`
      ka <- par_err_table[num_conc + 1,]$Estimate
      ka_se <- par_err_table[num_conc + 1,]$`Std. Error`
      R0[1:num_conc] <- par_err_table[(num_conc+2):(2*num_conc + 1), ]$Estimate
      R0_se[1:num_conc] <- par_err_table[(num_conc+2):(2*num_conc + 1), ]$`Std. Error`
      kd <- par_err_table[2*num_conc + 2,]$Estimate
      kd_se <- par_err_table[2*num_conc + 2,]$`Std. Error`
      if (bulkshift){
        Bulkshift[1:num_conc] <- par_err_table[(2*num_conc + 3):(3*num_conc+2), ]$Estimate
        Bulkshift_se[1:num_conc] <- par_err_table[(2*num_conc + 3):(3*num_conc+2), ]$`Std. Error`
      }
      
    }
    
    c(Rmax = Rmax, Rmax_se = Rmax_se, ka = ka, ka_se = ka_se, kd = kd, kd_se = kd_se, Bulkshift = Bulkshift, Bulkshift_se = Bulkshift_se, R0 = R0, R0_se = R0_se)
}