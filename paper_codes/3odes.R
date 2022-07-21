bivalent_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dL  <- -(2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2)
    dX1 <-  (2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2)
    dX2 <-  ka2*X1*L - 2*kd2*X2;
    
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

bivalent_full_ka0_objective <- function(pars, df, incl_concentrations, num_conc){
  
  Rmax <- pars[1:num_conc]
  
  #pars ("Rmax" one for each concentration,"ka", "tstart" one for each concentration)
  
  ka1 <- pars[1]
  ka2 <- pars[2]
  kd1 <- pars[3]
  kd2 <- pars[4]
  num_pars <- 4
  
  err <- NULL
  
  for (i in 1:num_conc){
    
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(AssocIndicator == 1) -> df_assoc
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    df_dissoc <- df_dissoc[seq(1, dim(df_dissoc)[1], 1),]
    # Association
    parameters <- c(ka1 = ka1,
                    kd1 = kd1,
                    ka2 = ka2,
                    kd2 = kd2,
                    Am = Concentration[i])
    
    state <- c(L  = pars[num_pars+i],
               X1 = 0,
               X2 = 0)
    RI <- pars[num_pars+num_conc+i]
    # state <- c(L  = pars[num_pars+1],
    #            X1 = 0,
    #            X2 = 0)
    # RI <- pars[num_pars+1+i]
    
    
    t_asc <- df_assoc$Time
    out <- ode(y = state, times = t_asc, func = bivalent_model, parms = parameters)
    y_asc <- out[,3] + out[,4] + RI

    # Dissociation
    t_dis <- df_dissoc$Time
    ##############################################################################
    parameters <- c(ka1 = ka1,
                    kd1 = kd1,
                    ka2 = ka2,
                    kd2 = kd2,
                    Am = 0)

    state <- c(out[length(t_asc),2],
               out[length(t_asc),3],
               out[length(t_asc),4])


    out <- ode(y = state, times = t_dis, func = bivalent_model, parms = parameters)
    y_dis <- out[,3] + out[,4] + RI
    ##############################################################################
    ##############################################################################
    solved_RU = c(y_asc, y_dis)
    RU <- c(df_assoc$RU, df_dissoc$RU)
    err <- c(err, RU - solved_RU)
    
  }
  err <- na.omit(err)
  err
}

fit_bivalent_full <- function(well_idx, sample_info, x_vals, y_vals,
                                                incl_concentrations_values){
  
  # this function will fit all concentrations for one well
  
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
  
  lower <- c(0, 0, 0, 0, rep(0, num_conc), rep(-Inf, num_conc))
  # lower <- c(0, 0, 0, 0, 0, rep(-Inf, num_conc))
  #upper <- c(1e6, 1e-1, 1e-2, 1e-2, Rmax_start$max)
  # ka1, ka2, kd1, kd2, Rmax1-5, R01-5
  
  ka1_IGs <- c(1e2, 1e3, 1e4)
  ka2_IGs <- c(1e-6, 1e-5, 1e-4)
  kd1_IGs <- c(1e-4, 1e-3, 1e-2)
  kd2_IGs <- c(1e-6, 1e-5, 1e-4)
  
  
  init_params_list <- NULL
  
  for (ka1_IG in ka1_IGs){
    for (ka2_IG in ka2_IGs){
      for (kd1_IG in kd1_IGs){
        for (kd2_IG in kd2_IGs){
          temp_IGs <- c(ka1_IG, ka2_IG, kd1_IG, kd2_IG, max(Rmax_start$max)*rep(1, num_conc), R0_start$min)
          # temp_IGs <- c(ka1_IG, ka2_IG, kd1_IG, kd2_IG, max(Rmax_start$max), R0_start$min)
          init_params_list <- rbind(init_params_list, temp_IGs)
        }
      }
    }
  }
  
  # fval <- 1e64
  # t0_start <- dissoc_start
  # estimated_params_list <- NULL
  # estimated_fval <- NULL
  
  best_error <- 1e64
  t0_start <- dissoc_start
  estimated_params_list <- NULL
  estimated_fval <- NULL
  
  res_R0 <- NULL
  error_val <- 2e64  
  
  # for (i in 1:dim(init_params_list)[1]){
  for (i in 1:5){
    # print(i)
    init_params <- init_params_list[i,]
    
    tryCatch(
      expr = {
        # Your code...
        # goes here...
        # ...
        
        # init_params <- c(1e3, 1e-5, 1e-3, 1e-5, max(Rmax_start$max)*rep(1, num_conc), R0_start$min)
        res_R0 <- nls.lm(init_params,fn = bivalent_full_ka0_objective, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc,
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
      best_res_R0 <- res_R0
      best_error <- error_val
    }
  }
  
  full_param_table <- estimated_fval
  full_param_table <- full_param_table %>% as_tibble() #&>& setNames()
  write_csv(full_param_table, paste('~/Summer2022/project/paper_codes/param_table_new/full_param',as.character(well_idx),'.csv',sep=''))
  
  # kd_param_table <- as.matrix(full_param_table[,c(1:4,15:18,29)])
  # kd_param_table <- kd_param_table %>% as_tibble() %>% setNames(c("ka1_0", "ka2_0", "kd1_0", "kd2_0", "ka1", "ka2", "kd1", "kd2", "SSE"))
  # write_csv(kd_param_table, paste('~/Summer2022/project/paper_codes/param_table_new/kd_param',
  #                                 as.character(well_idx),'.csv',
  #                                 sep=''))
  
  pars <- coefficients(best_res_R0)
  temp_pars <- pars
  pars[5:9] <- rep(temp_pars[5], num_conc)
  pars[10:14] <- temp_pars[6:length(temp_pars)]
  
  print(pars)
  ka <- pars[1:2]
  kd <- pars[3:4]
  Rmax <- pars[5:length(pars)]
  fit_outcomes <- get_bivalent_fit_outcomes(Rmax, ka, kd, df, num_conc,
                                   incl_concentrations,
                                   t0 = t0_start)
  
  # list("R0" = res_R0, "tstart" = res_tstart, "FitOutcomes" = fit_outcomes)
  list("R0" = best_res_R0, "FitOutcomes" = fit_outcomes)
}

get_bivalent_fit_outcomes <- function(Rmax, ka, kd, df, num_conc,
                             incl_concentrations,
                             t0){
  ka1 <- ka[1]
  ka2 <- ka[2]
  kd1 <- kd[1]
  kd2 <- kd[2]
  num_pars <- 4
  
  full_output_RU <- NULL
  
  for (i in 1:num_conc){
    
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(AssocIndicator == 1) -> df_assoc
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    df_dissoc <- df_dissoc[seq(1, dim(df_dissoc)[1], 1),]
    # Association
    parameters <- c(ka1 = ka1,
                    kd1 = kd1,
                    ka2 = ka2,
                    kd2 = kd2,
                    Am = Concentration[i])
    
    state <- c(L  = Rmax[i],
               X1 = 0,
               X2 = 0)
    RI <- Rmax[num_conc+i]
    t_asc <- df_assoc$Time
    out <- ode(y = state, times = t_asc, func = bivalent_model, parms = parameters)
    df_assoc$RU <- out[,3] + out[,4] + RI

    # Dissociation
    t_dis <- df_dissoc$Time
    ##############################################################################
    parameters <- c(ka1 = ka1,
                    kd1 = kd1,
                    ka2 = ka2,
                    kd2 = kd2,
                    Am = 0)

    state <- c(out[length(t_asc),2],
               out[length(t_asc),3],
               out[length(t_asc),4])
    out <- ode(y = state, times = t_dis, func = bivalent_model, parms = parameters)
    df_dissoc$RU <- out[,3] + out[,4] + RI
    ##############################################################################
    ##############################################################################
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
  #fit_df <- fits[[well_idx]]$FitOutcomes
  
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
  df <- na.omit(df)

  colnames(df) <- c("Time", "RU", "Concentration", "FittedRU")
  
  df$Concentration <- as_factor(df$Concentration)
  
  sub_title <- paste("Model Fitting")
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
  #                    sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  dissociation_start <- sample_info$Baseline[[well_idx]] + sample_info$`Bsl Start`[[well_idx]] + sample_info$Association[[well_idx]]
  
  ggplot(df, aes(x = Time, y = RU)) + geom_point(size = 0.09, aes(color = Concentration)) +
    ggtitle(ligand_desc, subtitle = sub_title) +
    geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1) +
    geom_vline(xintercept = dissociation_start, linetype="dashed", color = "black", size=0.5)
}