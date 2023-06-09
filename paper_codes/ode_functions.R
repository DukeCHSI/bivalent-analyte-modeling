########################################################################################
# ODE Models
########################################################################################
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

objective_function <- function(pars, df, incl_concentrations, num_conc, param_names, 
                               use_RI,
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind,
                               fixed_info=NULL){
  
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
      sim_results <- run_model(pars, df, num_conc, incl_concentrations, param_names,
                               use_RI,
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind,
                               fixed_info)
      
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
                                         fixed_sheet=NULL,
                                         min_allowed_kd = 10^(-5),
                                         max_iterations = 500,
                                         ptol = 10^(-10),
                                         ftol = 10^(-10)){
  
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
  df %>% group_by(Concentration) %>% summarise(min = min(RU, na.rm = TRUE)) -> R0_start
  
  ###########################  ###########################  ###########################
  # Manual modification for options
  ###########################  ###########################  ###########################
  use_globalRmax <- sample_info[well_idx,]$`Global Rmax`
  use_bulkShift <- sample_info[well_idx,]$Bulkshift
  use_regeneration <- sample_info[well_idx,]$Regen.
  use_secondRebind <- sample_info[well_idx,]$`Sec. Rebinding`
  model <- sample_info[well_idx,]$Model
  
  #
  use_RI <- "N"
  fixed_info <- NULL
  
  ###########################  ###########################  ###########################
  ###########################  ###########################  ###########################
  
  if (use_globalRmax == "N"){
    num_Rmax <- num_conc
  }
  else{
    num_Rmax <- 1
  }
  
  if (use_RI == "Y"){
    num_R0 <- num_conc
    init_R0s <- R0_start$min
  }
  else{
    num_R0 <- 0
    init_R0s <- NULL
  }
  
  if (use_regeneration == "N"){
    num_tstar <- num_conc
  }
  else{
    num_tstar <- 0
  }
  
  if (model == 'bivalentAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(0, num_R0), rep(0, num_tstar))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep('R0', num_R0), rep('t_star', num_tstar))
  }
  
  ka1_IGs <- c(1e2, 1e3, 1e4)
  ka2_IGs <- c(1e-6, 1e-5, 1e-4)
  kd1_IGs <- c(1e-4, 1e-3, 1e-2)
  kd2_IGs <- c(1e-6, 1e-5, 1e-4)
  
  
  init_params_list <- NULL
  t0_count <- 1
  for (ka1_IG in ka1_IGs){
    for (ka2_IG in ka2_IGs){
      for (kd1_IG in kd1_IGs){
        for (kd2_IG in kd2_IGs){
          # Kinetics parameters, Global/local Rmaxs, Refractive index, Non-regenerative
          # (ka1, ka2, kd1, kd2), (Rmax's), (R0s), (tstar)
          temp_IGs <- c(ka1_IG, ka2_IG, kd1_IG, kd2_IG, max(Rmax_start$max)*rep(1, num_Rmax), init_R0s, runif(num_tstar, min=0, max=120))
          init_params_list <- rbind(init_params_list, temp_IGs)
        }
      }
    }
  }
  
  best_error <- 1e64
  t0_start <- dissoc_start
  estimated_params_list <- NULL
  estimated_fval <- NULL
  
  res_R0 <- NULL
  error_val <- 2e64 
  solved_success <- FALSE
  
  time_start <- Sys.time()
  
  num_run <- dim(init_params_list)[1]
  
  for (i in 1:num_run){
    
    init_params <- init_params_list[i,]
    
    print("Initial Params:")
    print(init_params)
    
    tryCatch(
      expr = {
        res_R0 <- nls.lm(init_params, fn = objective_function, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc,
                         param_names = param_names,
                         use_RI = use_RI,
                         use_regeneration = use_regeneration, 
                         model = model, 
                         use_globalRmax = use_globalRmax,
                         use_secondRebind = use_secondRebind,
                         fixed_info = fixed_info,
                         control = nls.lm.control(maxiter = 1000),
                         lower = lower)
        
        print("Estimated Params:")
        print(coefficients(res_R0))
        
        solved_success <- TRUE
      },
      error = function(cond){
        solved_success <- FALSE
      }
    )
    
    
    if (solved_success == TRUE){
      estimated_fval <- rbind(estimated_fval, res_R0$deviance)
      estimated_params_list <- rbind(estimated_params_list, coefficients(res_R0))
      error_val <- res_R0$deviance
    }
    else{
      res_R0$deviance <- 2e64
      estimated_fval <- rbind(estimated_fval, res_R0$deviance)
      estimated_params_list <- rbind(estimated_params_list, init_params)
      error_val <- res_R0$deviance
    }
    
    if (error_val < best_error){
      best_res_R0 <- res_R0
      best_error <- error_val
    }
    solved_success <- FALSE
  }
  
  # full_param_table <- cbind(init_params_list, estimated_params_list, estimated_fval)
  # full_param_table <- full_param_table %>% as_tibble() #&>& setNames()
  # write_csv(full_param_table, paste('~/Summer2022/project/hiv-summer-2022/paper_codes/param_table_old_short2/full_param',as.character(well_idx),'.csv',sep=''))
  
  # kd_param_table <- as.matrix(full_param_table[,c(1:4,7:12,23)])
  # kd_param_table <- kd_param_table %>% as_tibble() %>% setNames(c("ka1_0", "ka2_0", "kd1_0", "kd2_0", "ka1", "ka2", "kd1", "kd2", "SSE"))
  # write_csv(kd_param_table, paste('~/Summer2022/project/paper_codes/param_table_old/k_param',
  #                                 as.character(well_idx),'.csv',
  #                                 sep=''))
  
  time_end <- Sys.time()
  elapsed_time <- time_end - time_start
  print("Time:")
  print(elapsed_time)
  
  pars <- coefficients(best_res_R0)
  print("Best Params:")
  print(pars)
  
  
  # load(file='bivalent_oldResults_t0_short2_nonRegen.Rdata')
  # pars <- fits_list[[well_idx]]$R0$par
  # best_res_R0 <- fits_list[[well_idx]]$R0
  
  fit_outcomes <- run_model(pars, df, num_conc,
                            incl_concentrations,
                            param_names,
                            use_RI,
                            use_regeneration,
                            model,
                            use_globalRmax,
                            use_secondRebind,
                            fixed_info)
  
  list("R0" = best_res_R0, "FitOutcomes" = fit_outcomes)
}



run_model <- function(pars, df, num_conc, incl_concentrations, param_names, 
                      use_RI,
                      use_regeneration, 
                      model, 
                      use_globalRmax,
                      use_secondRebind,
                      fixed_info=NULL){
  
  full_output_RU <- NULL
  # print("Print t_star")
  # print(pars)
  
  for (i in 1:num_conc){
    
    df_i <- df %>% filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration
    
    df_i %>% filter(AssocIndicator == 1) -> df_assoc
    df_i %>% filter(DissocIndicator == 1) -> df_dissoc
    df_dissoc <- df_dissoc[seq(1, dim(df_dissoc)[1], 1),]
    
    if (model == 'bivalentAnalyte'){
      
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
      
      if (use_RI == "Y"){
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
                      Am = incl_concentrations[i])
      asc_indicator <- df_assoc$AssocIndicator
      
      if (use_regeneration == "N"){
        t_pred_asc <- NULL
        t_stars <- pars[param_names == "t_star"]
        t_star <- t_stars[i]
        if (t_star >= 1){
          t0 <- df_assoc$Time[1] - t_star
          t_pred_asc <- seq(from = t0, df_assoc$Time[1]-1, by = 1)
          
          t_asc <- c(t_pred_asc, df_assoc$Time)
          t_asc <- t_asc - t0
          
          asc_indicator <- c(rep(0, length(t_pred_asc)), asc_indicator)
          
          out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
          df_assoc$RU <- out[asc_indicator == 1, 3] + out[asc_indicator == 1, 4] + RI
        }
        else{
          t0 <- df_assoc$Time[1]
          t_asc <- c(t_pred_asc, df_assoc$Time)
          t_asc <- t_asc - t0
          
          asc_indicator <- c(rep(0, length(t_pred_asc)), asc_indicator)
          
          out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
          df_assoc$RU <- out[asc_indicator == 1, 3] + out[asc_indicator == 1, 4] + RI
        }
      }
      else{
        t0 <- df_assoc$Time[1]
        t_asc <- df_assoc$Time
        t_asc <- t_asc - t0
        
        out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
        df_assoc$RU <- out[, 3] + out[, 4] + RI
      }
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], df_dissoc$Time - t0)
      
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
  
  
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
  #                    sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  dissociation_start <- sample_info$Baseline[[well_idx]] + sample_info$`Bsl Start`[[well_idx]] + sample_info$Association[[well_idx]]
  
  ### For sensorgrams
  sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
                     sample_info[well_idx,]$Row,
                     "Column", sample_info[well_idx,]$Column)
  ggplot(df, aes(x = Time, y = RU)) + geom_point(size = 1, aes(color = Concentration)) +
    # ggtitle(ligand_desc, subtitle = sub_title) +
    guides(color = guide_legend(reverse=TRUE, override.aes = list(size = 5))) +
    geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1.5) +
    geom_vline(xintercept = dissociation_start, linetype="dashed", color = "black", size=1) +
    theme(text = element_text(size=10)) +
    labs(color = "Concentration (M)") +
    ylab("Response Unit (RU)") + xlab("Time (s)")
  
  ### For manuscript plots
  # sub_title <- paste("Bivalent Analyte Model with Standard Length of Dissociation")
  # ggplot(df, aes(x = Time, y = RU)) + geom_point(size = 1, aes(color = Concentration)) +
  #   guides(color = guide_legend(reverse=TRUE, override.aes = list(size = 5))) +
  #   # ggtitle("CH31 mAb binding to CH505 T/F gp120", subtitle = sub_title) +
  #   geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1.5) +
  #   geom_vline(xintercept = dissociation_start, linetype="dashed",
  #              color = "black", size=1) +
  #   theme(text = element_text(size=15)) +
  #   labs(color = "Concentration (M)") +
  #   ylab("Response Unit (RU)") + xlab("Time (s)")
}