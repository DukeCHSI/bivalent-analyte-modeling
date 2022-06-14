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

bispecificLigandMixture_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dL1   <- - ka1*A1m*L1 + kd1*A1L1
    dL2   <- - ka2*A1m*L2 + kd2*A2L2
    dA1L1 <-   ka1*A2m*L1 - kd1*A1L1 + k12*L1*AL2
    dA2L2 <-   ka2*A2m*L2 - kd2*A2L2 + k21*L2*AL1
    
    #return the rate of change
    list(c(dL1, dL2, dA1L1, dA2L2))
  })
}

bispecificLigandHeterogenousAnalyte_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    # dL1   <- - ka1*A1m*L1 + kd1*A1L1
    # dL2   <- - ka2*A1m*L2 + kd2*A2L2
    # dA1L1 <-   ka1*A2m*L1 - kd1*A1L1
    # dA2L2 <-   ka2*A2m*L2 - kd2*A2L2
    #return the rate of change
    # list(c(dL1, dL2, dA1L1, dA2L2))
    
    dL <- - ka1*A1m*L + kd1*A1L - ka2*A2m*L + kd2*A2L
    dA1L <- ka1*A1m*L - kd1*A1L - ka2*A1L*A2m + kd2*A1LA2
    dA2L <- ka2*A2m*L - kd2*A2L - ka1*A2L*A1m + ka1*A1LA2
    dA1LA2 <- ka1*A1m*A2L - kd1*A1LA2 + ka2*A1L*A2m - kd2*A1LA2
    
    #return the rate of change
    list(c(dL, dA1L, dA2L, dA1LA2))
  })
}


objective_function <- function(pars, df, incl_concentrations, num_conc, use_regeneration, model){
  
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
      sim_results <- run_model(pars, df, num_conc, incl_concentrations, use_regeneration, model)
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

fit_heterogenousLigand_full <- function(well_idx, sample_info, x_vals, y_vals,
                                        incl_concentrations_values){
  
  # this function will fit all concentrations for one well
  use_bulkShift <- sample_info[well_idx,]$Bulkshift
  use_regeneration <- sample_info[well_idx,]$Regen.
  use_globalRmax <- sample_info[well_idx,]$`Global Rmax`
  
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
  }
  else if (model == 'bivalentAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'monovalent'){
    lower <- c(0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'heterogenousAnalyte'){
    lower <- c(0, 0, 0, 0, 0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'twoState'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'bivalentAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'bispecificLigandMixture'){
    lower <- c(0, 0, 0, 0, -Inf, -Inf, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  else if (model == 'bispecificLigandHeterogenousAnalyte_model'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(-Inf, num_R0))
  }
  
  
  
  
  init_params_list <- NULL
  
  best_error <- 1e64
  t0_start <- dissoc_start
  estimated_params_list <- NULL
  estimated_fval <- NULL
  
  time_start <- Sys.time()
  
  num_run <- sample_info[well_idx,]$`Num Run`
  
  for (i in 1:num_run){
    
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
      
      # ka1, kd2, kd1, kd2, mw1, mw2, n1, n2, Rmax_1-5, RU0_1-5
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
    else if (model == 'bispecificLigandHeterogenousAnalyte'){
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
    
    # print("Initial Parameters:")
    # print(init_params)
    
    tryCatch(
      expr = {
        # Your code...
        # goes here...
        # ...
        res_R0 <- nls.lm(init_params, fn = objective_function, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc, 
                         use_regeneration = use_regeneration, model = model,
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
  
  fit_outcomes <- run_model(pars, df, num_conc,
                            incl_concentrations,
                            use_regeneration,
                            model)
  
  list("R0" = best_res_R0, "FitOutcomes" = fit_outcomes)
}

run_model <- function(pars, df, num_conc, incl_concentrations, use_regeneration, model){
  
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
      
      ka1 <- pars[1]
      ka2 <- pars[2]
      kd1 <- pars[3]
      kd2 <- pars[4]
      num_pars <- 4
      if (use_regeneration == "N"){
        RI <- pars[num_pars+2*num_conc+i]
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
      
      state <- c(B1  = pars[num_pars+i],
                 B2  = pars[num_pars+num_conc+i],
                 AB1 = 0,
                 AB2 = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = heterogenousLigand_model, parms = parameters)
      df_assoc$RU <- out[,4] + out[,5] + RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      ka2 = 0,
                      kd2 = kd2,
                      Am = Concentration[i])
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4],
                 out[length(t_asc),5])
      
      out <- ode(y = state, times = t_dis, func = heterogenousLigand_model, parms = parameters)
      df_dissoc$RU <- out[,4] + out[,5] + RI
      
    }
    else if (model == 'bivalentAnalyte'){
      
      ka1 <- pars[1]
      ka2 <- pars[2]
      kd1 <- pars[3]
      kd2 <- pars[4]
      num_pars <- 4
      if (use_regeneration == "N"){
        RI <- pars[num_pars+num_conc+i]
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
      
      state <- c(L  = pars[num_pars+i],
                 X1 = 0,
                 X2 = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = 0)
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = bivalentAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[,3] + out[,4] + RI
      
    }
    else if (model == 'monovalent'){
      
      ka1 <- pars[1]
      kd1 <- pars[2]
      num_pars <- 2
      if (use_regeneration == "N"){
        RI <- pars[num_pars+num_conc+i]
      }
      else{
        RI <- 0
      }
      
      # Association
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      Am = Concentration[i])
      
      state <- c(L  = pars[num_pars+i],
                 X1 = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = monovalent_model, parms = parameters)
      df_assoc$RU <- out[,3] + RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      Am = Concentration[i])
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3])
      
      out <- ode(y = state, times = t_dis, func = monovalent_model, parms = parameters)
      df_dissoc$RU <- out[,3] + RI
      
    }
    else if (model == 'heterogenousAnalyte'){
      
      ka1 <- pars[1]
      ka2 <- pars[2]
      kd1 <- pars[3]
      kd2 <- pars[4]
      mw1 <- pars[5]
      mw2 <- pars[6]
      n1  <- pars[7]
      n2  <- pars[8]
      num_pars <- 8
      
      if (use_regeneration == "N"){
        RI <- pars[num_pars+num_conc+i]
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
      
      state <- c(B   = pars[num_pars+i],
                 A1B = 0,
                 A2B = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = heterogenousAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] +  RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      ka2 = 0,
                      kd2 = kd2,
                      mw1 = mw1,
                      mw2 = mw2,
                      n1  = n1,
                      n2  = n2,
                      Am1 = Concentration[i], ############## Need to change this to reflect the info file
                      Am2 = Concentration[i]) ############## Need to change this to reflect the info file
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = heterogenousAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[,3] + out[,4] +  RI
      
    }
    else if (model == 'twoState'){
      
      ka1 <- pars[1]
      ka2 <- pars[2]
      kd1 <- pars[3]
      kd2 <- pars[4]
      num_pars <- 4
      if (use_regeneration == "N"){
        RI <- pars[num_pars+num_conc+i]
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
      
      state <- c(B   = pars[num_pars+i],
                 AB  = 0,
                 ABx = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = twoState_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      ka2 = 0,
                      kd2 = kd2,
                      Am = Concentration[i])
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4])
      
      out <- ode(y = state, times = t_dis, func = twoState_model, parms = parameters)
      df_dissoc$RU <- out[,3] + out[,4] + RI
      
    }
    else if (model == 'bispecificLigandHeterogenousAnalyte'){
      
      ka1 <- pars[1]
      ka2 <- pars[2]
      kd1 <- pars[3]
      kd2 <- pars[4]
      
      num_pars <- 4
      
      if (use_regeneration == "N"){
        RI <- pars[num_pars+num_conc+i]
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
      
      state <- c(L   = pars[num_pars+i],
                 A1L = 0,
                 A2L = 0,
                 A1LA2 = 0)
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = bispecificLigandHeterogenousAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + out[,5] +  RI
      
      # Dissociation
      t_dis <- df_dissoc$Time
      parameters <- c(ka1 = 0,
                      kd1 = kd1,
                      ka2 = 0,
                      kd2 = kd2,
                      Am1 = Concentration[i],
                      Am2 = Concentration[i])
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4],
                 out[length(t_asc),5])
      
      out <- ode(y = state, times = t_dis, func = bispecificLigandHeterogenousAnalyte_model, parms = parameters)
      df_dissoc$RU <- out[,3] + out[,4] + out[,5] +  RI
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