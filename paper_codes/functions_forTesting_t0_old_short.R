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
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind,
                               fixed_info=NULL){
  
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
  
  # this function will fit all concentrations for one well
  use_bulkShift <- sample_info[well_idx,]$Bulkshift
  use_regeneration <- sample_info[well_idx,]$Regen.
  use_globalRmax <- sample_info[well_idx,]$`Global Rmax`
  # use_secondRebind <- sample_info[well_idx,]$`Sec. Rebinding`
  
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
  use_globalRmax <- "N"
  use_bulkShift <- "N"
  use_regeneration <- "Y" #sample_info[well_idx,]$Regen.
  
  # model <- sample_info[well_idx,]$Model
  model <- 'bivalentAnalyte_tstar'
  
  fixed_info <- NULL
  # fixed_info$CH31 <- fixed_sheet$CH31[well_idx,]
  # fixed_info$PGT121 <- fixed_sheet$PGT121[well_idx,]
  
  use_secondRebind <- "N" #sample_info[well_idx,]$`Sec. Rebinding`
  ###########################  ###########################  ###########################
  ###########################  ###########################  ###########################
  
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
  
  if (model == 'bivalentAnalyte'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(0, num_R0), rep(0, num_conc))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep('t_star', num_conc))
  }
  else if (model == 'bivalentAnalyte_tstar'){
    lower <- c(0, 0, 0, 0, rep(0, num_Rmax), rep(0, num_R0), rep(0, num_conc))
    param_names <- c("ka1", "ka2", "kd1", "kd2", rep("Rmax", num_Rmax), rep('t_star', num_conc))
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
          temp_IGs <- c(ka1_IG, ka2_IG, kd1_IG, kd2_IG, max(Rmax_start$max)*rep(1, num_Rmax), runif(num_conc, min=0, max=120))
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
  
  num_run <- dim(init_params_list)[1] #Esample_info[well_idx,]$`Num Run`
  
  # num_run <- 2
  # init_params_list <- init_params_list[1:num_run,]
  
  for (i in 1:num_run){
    
    init_params <- init_params_list[i,]
    
    print("Initial Params:")
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
  
  full_param_table <- cbind(init_params_list, estimated_params_list, estimated_fval)
  full_param_table <- full_param_table %>% as_tibble() #&>& setNames()
  write_csv(full_param_table, paste('~/Summer2022/project/hiv-summer-2022/paper_codes/param_table_old_short/full_param',as.character(well_idx),'.csv',sep=''))
  
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
  # load(file='Results_2Fixed_FixedRmax.Rdata')
  
  
  # pars <- fits_list[[well_idx]]$R0
  
  fit_outcomes <- run_model(pars, df, num_conc,
                            incl_concentrations,
                            param_names,
                            use_regeneration,
                            model,
                            use_globalRmax,
                            use_secondRebind,
                            fixed_info)
  
  list("R0" = best_res_R0, "FitOutcomes" = fit_outcomes)
}



run_model <- function(pars, df, num_conc, incl_concentrations, param_names, 
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
      }else{
        state <- c(L  = Rmaxs[1],
                   X1 = 0,
                   X2 = 0)
      }
      
      if (use_regeneration == "N"){
        RIs <- pars[param_names == "R0"]
        RI <- RIs[i]
      }else{
        RI <- 0
      }
      
      
      ########### Pre-Association ########### 
      # For Pre-association part, if t_star/t_adjusted > 1, then solve ODE from 0 to t_star/t_adjusted
      # the obtained output from this pre-Association phase will be used as the initial conditions for Association phase
      # 
      t_adjusted <- pars[param_names == "t_star"]
      if (t_adjusted[i] > 1){
        
        parameters <- c(ka1 = ka1,
                        kd1 = kd1,
                        ka2 = ka2,
                        kd2 = kd2,
                        Am = incl_concentrations[i])
        
        # t0 <- df_assoc$Time[1] - t_adjusted
        t_pred_assoc <- seq(from = 0, t_adjusted[i], by = 1)
        # t_pred_assoc <- t_pred_assoc - t0
        
        out <- ode(y = state, times = t_pred_assoc, func = bivalentAnalyte_model, parms = parameters)
        
        state <- c(out[length(t_pred_assoc),2],
                   out[length(t_pred_assoc),3],
                   out[length(t_pred_assoc),4])
      }
      
      #
      
      ########### Association ########### 
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      Am = incl_concentrations[i])
      
      t_asc <- df_assoc$Time
      out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
      df_assoc$RU <- out[,3] + out[,4] + RI
      
      ########### Dissociation ########### 
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
    else if (model == 'bivalentAnalyte_tstar'){
      
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
                      Am = incl_concentrations[i])
      
      
      t_pred_asc <- NULL
      asc_indicator <- df_assoc$AssocIndicator
      t_stars <- pars[param_names == "t_star"]
      t_star <- t_stars[i]
      if (t_star != 0){
        t0 <- df_assoc$Time[1] - t_star
        # t0 <- floor(t0)
        t_pred_asc <- seq(from = t0, df_assoc$Time[1]-1, by = 1)
        # t_pred_asc <- t_pred_asc - t0
        
        t_asc <- c(t_pred_asc, df_assoc$Time)
        t_asc <- t_asc - t0
        
        asc_indicator <- c(rep(0, length(t_pred_asc)), asc_indicator)
        
        out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
        df_assoc$RU <- out[asc_indicator == 1, 3] + out[asc_indicator == 1, 4] + RI
      }
      else{
        t_asc <- c(t_pred_asc, df_assoc$Time)
        t_asc <- t_asc - t0
        
        asc_indicator <- c(rep(0, length(t_pred_asc)), asc_indicator)
        
        out <- ode(y = state, times = t_asc, func = bivalentAnalyte_model, parms = parameters)
        df_assoc$RU <- out[asc_indicator == 1, 3] + out[asc_indicator == 1, 4] + RI
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
  
  sub_title <- paste("Bivalent Analyte Model-1 with Extended Length of Dissociation")
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













########################################################################################
# Helper function
########################################################################################

process_sample_sheet <- function(sample_sheet_path, rows_dict, sort_order="PairRow"){
  # Process Sample Sheet
  
  # Position/Channel/Sensor has multiple wells. Need to expand to one row per well.
  sample_sheet <- read_csv(sample_sheet_path)
  sample_sheet %>% separate(`Position/Channel/Sensor`, into = c("BeginPosition", "EndPosition", sep = "-")) %>%
    #  filter(!is.na(EndPosition)) %>%
    separate(BeginPosition, into = c("Row", "Column"), sep = 1) %>%
    separate(EndPosition, into = c("EndRow", "EndColumn"), sep = 1) %>%
    mutate(Column = as.integer(Column)) %>%
    mutate(EndColumn = as.integer(EndColumn)) -> RowExpansionTemplate
  
  length_ss <- dim(RowExpansionTemplate)[1]
  ss_exp <- NULL
  
  for(i in 1:length_ss){
    
    old_row <- RowExpansionTemplate[i,]
    old_row$Row <- str_to_upper(old_row$Row)                                                    
    old_row$EndRow <- str_to_upper(old_row$EndRow)
    
    start_letter <- old_row$Row
    end_letter   <- old_row$EndRow
    
    if (is.na(end_letter))
      end_letter <- start_letter
    
    ####################################################################### For cases, A1B2??? How would it look like?
    if(start_letter != end_letter) {
      start_idx <- which(LETTERS == start_letter) #                         Was LETTERS defined somewhere?                        
      end_idx <- which(LETTERS == end_letter)     #                         start_letter, end_letteer?      
      for (j in start_idx+1:(end_idx-1)){
        new_row <- old_row
        new_row$Row <- LETTERS[j]
        if (j == start_idx+1)
          ss_exp <- rbind(old_row, new_row)
        else
          ss_exp <- rbind(ss_exp, new_row)
        print(ss_exp)
      }
    }
    else
      ss_exp <- rbind(ss_exp, old_row)
  }
  
  ### Now to expand columns if necessary
  RowExpansionTemplate <- ss_exp
  
  length_ss <- dim(RowExpansionTemplate)[1]
  ss_exp <- NULL
  
  for(i in 1:length_ss){
    old_row <- RowExpansionTemplate[i,]
    
    start_col <- old_row$Column
    end_col   <- old_row$EndColumn
    
    if (is.na(end_col))
      end_col <- start_col
    
    if (start_col != end_col){
      for (j in (start_col+1):end_col){
        new_row <- old_row
        new_row$Column++
          if (j == start_idx+1)
            ss_exp <- rbind(old_row, new_row)
        else
          ss_exp <- rbind(ss_exp, new_row)
      }
    }
    else
      ss_exp <- rbind(ss_exp, RowExpansionTemplate[i,])
  }
  
  ss_exp %>% select(-EndRow, -EndColumn, -"-") -> ss_exp
  
  # Now expand blocks
  current_idx <- 0
  new_current_idx <- 0
  nsamples <- dim(ss_exp)[1]
  sample_info_expanded <- NULL
  
  for(i in 1:nsamples){
    
    num_blocks <- str_count(ss_exp$`Block/Chip/Tray`[i], ";") + 1
    
    blocks <- as.numeric(flatten(str_split(ss_exp$`Block/Chip/Tray`[i], ";")))
    
    # create an expanded sample sheet with one row for each block
    
    tmp_sample_info_expanded <- purrr::map_dfr(seq_len(num_blocks), ~ss_exp[i,])
    tmp_sample_info_expanded$Block <- blocks
    #  tmp_sample_info_expanded <- data.frame(tmp_sample_info_expanded)
    sample_info_expanded <- bind_rows(sample_info_expanded, tmp_sample_info_expanded)
    
  }
  # Rows are sorted A,E,B,F,C,G,D,H
  sample_info_expanded$RowSort <- translate_rows_cols_for_sort(sample_info_expanded$Row, 
                                                               sample_info_expanded$Column,
                                                               rows_dict,
                                                               sort_order)
  
  #sample_info_expanded %>% mutate(RowSort = translate_rows_for_sort(Row)) ->
  #  sample_info_expanded
  
  sample_info_expanded %>% arrange(RowSort, Column, Block) -> sample_info_expanded
  
  sample_info_expanded
}

#this is probably to be deleted. We don't actually use the ligand_conc data at all.
#process_carterra_output <- function(data_file_path){

# Each ligand/conc has two columns x and y (time and RU)
# The ligand/conc info is the first header line. The second header is just x,y repeated for 
# each ligand/conc combo.

# ligand_conc <- as.character(read_excel(data_file_path, col_names = FALSE, n_max = 1)[1,])
# ligand_conc <- tibble(Name = ligand_conc[which(ligand_conc != "NA")])

#ligand_conc %>% separate(Name, into =  c("Ligand", "Rest"), sep = "-") %>%
#  separate(Rest, into = c("quote", "Antigen", "Conc_str", "Concentration", "Rest"), sep = " ") %>% 
# select(-"quote", -"Conc_str", -"Rest") %>% 
#  separate(Ligand, into = c("Rest", "Ligand", "Rest2"), sep = " ") %>% 
#  select(-Rest, -Rest2) %>% separate(Concentration, into = c("Concentration", "Rest")) %>% 
#    select(-Rest) -> ligand_conc

#  ligand_conc
#}

select_samples <- function(sample_info, titration_data, remove_dissociation=NULL){
  
  sum_sum <- 0
  
  remove_ligands <- which(sample_info$Incl. == "N")
  keep_ligands <- which(sample_info$Incl. == "Y")
  
  ####################################################################### -starts_with("X"): What does "-" do?
  titration_data %>% select(everything(), -starts_with("Y")) -> x_vals
  titration_data %>% select(everything(), -starts_with("X")) -> y_vals
  
  n_time_points <- dim(x_vals)[1]
  nsamples <- dim(sample_info)[1]
  
  #loop over all samples to be included. Exclude concentrations not selected for analysis.
  
  displacement_per_ligand <- 0
  keep_concentrations_all <- NULL
  all_concentrations_ligand <- NULL
  all_concentrations_values <- NULL
  
  # baseline_dur <- unique(sample_info$Baseline)
  # association_dur <- unique(sample_info$Association)
  # dissociation_dur <- unique(sample_info$Dissociation)
  
  # if 
  
  for(i in 1:nsamples){
    # Count the number of blocks
    sum_sum <- sum_sum + nchar(sample_info$Block[i])
    
    num_conc <- str_count(sample_info$`All Concentrations`[i], ",") + 1
    
    
    if (sample_info$Incl.[i] == "N"){
      displacement_per_ligand <- i*num_conc
      next
    }
    
    concentrations <- 
      as.numeric(flatten(str_split(sample_info$`All Concentrations`[i], ",")))
    all_concentrations_ligand <- c(all_concentrations_ligand, num_conc)
    all_concentrations_values <- c(all_concentrations_values, concentrations)
    
    #Use this to match to ligand_conc and to x_vals, y_vals
    keep_conc_idx <- 1:num_conc + displacement_per_ligand
    
    keep_concentrations_all <- c(keep_concentrations_all, keep_conc_idx)
    
    displacement_per_ligand <- i*num_conc
    
  }
  # 
  # baseline_time <- sample_info$Baseline
  
  # time and ru values for selected wells
  x_vals_select <- x_vals[ , keep_concentrations_all]
  y_vals_select <- y_vals[ , keep_concentrations_all]
  
  # if (is.null(remove_dissociation) == FALSE)
  # {
  #   for(i in 1:nsamples){
  #       # Count the number of blocks
  #       
  #   }
  # }
  
  sample_info <- sample_info[keep_ligands, ]
  sample_info$NumConc <- all_concentrations_ligand  ############################# Why do we include all concentration ligand?
  
  list(Time = x_vals_select, RU = y_vals_select, 
       sample_info = sample_info,
       all_concentrations_values = all_concentrations_values,
       n_time_points = n_time_points)
}

select_concentrations <- function(sample_info, x_vals, y_vals){
  #after we select the wells, we select the concentrations out of the remaining time series
  
  n_time_points <- dim(x_vals)[1]
  nsamples <- dim(sample_info)[1]
  
  #loop over all samples to be included. Exclude concentrations not selected for analysis.
  
  displacement_per_ligand <- 0
  keep_concentrations <- NULL
  incl_concentrations_ligand <- NULL
  incl_concentrations_values <- NULL
  
  for(i in 1:nsamples){
    
    num_conc <- str_count(sample_info$`All Concentrations`[i], ",") + 1
    
    if (sample_info$Incl.[i] == "N"){
      displacement_per_ligand <- i*num_conc
      next
    }
    
    concentrations <- 
      as.numeric(flatten(str_split(sample_info$`All Concentrations`[i], ",")))
    
    if (str_to_upper(sample_info$`Incl. Conc.`[i]) == "ALL"){
      incl_concentrations <- concentrations
    } else {
      incl_concentrations <- 
        as.numeric(flatten(str_split(sample_info$`Incl. Conc.`[i], ",")))
    }
    num_incl <- length(incl_concentrations)
    
    if (num_incl > 5){
      incl_concentrations <- 
        get_best_window(i, sample_info, x_vals, y_vals,
                        num_incl, concentrations,
                        displacement_per_ligand)
      num_incl <- length(incl_concentrations)
    }
    
    incl_concentrations_ligand <- c(incl_concentrations_ligand, num_incl)
    incl_concentrations_values <- c(incl_concentrations_values, incl_concentrations)
    
    #Use this to match to ligand_conc and to x_vals, y_vals
    incl_idx <- which(concentrations %in% incl_concentrations) + displacement_per_ligand
    keep_concentrations <- c(keep_concentrations, incl_idx)
    
    # the following keeps all of the concentrations for each of the wells selected to include, but
    # excludes the wells not selected
    
    displacement_per_ligand <- i*num_conc
    
  }
  
  sample_info$NumInclConc <- incl_concentrations_ligand
  
  list(keep_concentrations = keep_concentrations, 
       sample_info = sample_info, incl_concentrations_values = incl_concentrations_values)
  
}


find_dissociation_window <- function(well_idx, sample_info, x_vals, y_vals, 
                                     all_concentrations_ligand, max_RU_tol, min_RU_tol){
  
  baseline <- sample_info[well_idx,]$Baseline
  association <- sample_info[well_idx,]$Association
  
  start_time <- baseline + association
  dissociation <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- all_concentrations_ligand[well_idx]
  end_idx <- start_idx + num_conc - 1
  
  df_all <- NULL
  max_idx <- 0
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    #    RU <- log(as_vector(y_vals[, i]))
    RU <- y_vals[, i]
    df <- bind_cols(Time, RU)
    names(df) <- c("Time", "RU")
    df %>% filter(Time > (start_time - 100)) -> df
    
    # Do not base on low information concentrations
    if (mean(df$RU, na.rm = TRUE) < min_RU_tol | mean(df$RU, na.rm = TRUE) > max_RU_tol)
      next
    df_zoo <- as.zoo(df)
    max_idx <- max_idx + 1
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
    df_out %>% mutate(x = rep(max_idx, n_vals)) -> df_out
    df_out %>% mutate(RollIndex = 1:n_vals) -> df_out
    df_all <- bind_rows(df_out, df_all)
    
    df %>% ggplot(aes(x = Time, y = RU)) + geom_point() + 
      geom_abline(slope = df_out$Slope[1], intercept = df_out$Intercept[1]) +
      geom_abline(slope = df_out$Slope[10], intercept = df_out$Intercept[10], color = "purple") +
      geom_abline(slope = df_out$Slope[100], intercept = df_out$Intercept[100], color = "yellow") +
      geom_abline(slope = df_out$Slope[200], intercept = df_out$Intercept[200], color = "green") +
      geom_abline(slope = df_out$Slope[300], intercept = df_out$Intercept[300], color = "blue") +
      geom_abline(slope = df_out$Slope[400], intercept = df_out$Intercept[400], color = "red") +
      geom_abline(slope = df_out$Slope[410], intercept = df_out$Intercept[410], color = "orange") 
  }
  
  
  max_list[[max_idx]] %>% as_tibble() %>% mutate(x = 1:list_length) %>%
    ggplot(aes(x = x, y = abs(Time))) + geom_point()
  
  df <- bind_cols("Time" = Time, "RU" = RU)
  colnames(df) <- c("Time", "RU")
  
}

translate_rows_cols_for_sort <- function(rows_list, columns_list, dict, sort_order="PairRow"){
  sort_list <- NULL
  # Get unique column values
  unique_col <- unique(columns_list)
  # Get the largest column value
  num_columns <- max(unique_col)
  
  
  
  dict <- dict[order(unlist(dict))]
  dict_values <- as.numeric(dict)
  dict_names <- names(dict)
  
  # Create a 1-2 dictionary
  odd_even_dict <- dict
  
  for (i in 1:length(dict)){
    if (dict_values[i] %% 2 == 1)
      odd_even_dict[i] = 1
    else
      odd_even_dict[i] = 2
  }
  
  if (tolower(sort_order) == "roi"){
    # Assume wells go in pairs [A1,E1],...,[A4,E4],[B1,F1],..,[B4,F4], etc.
    n_jump <- num_columns*2
    
    # Create sort values for rows based on their orders
    A_list <- seq(odd_even_dict$A+(ceiling(dict$A/2)-1)*n_jump,
                  (ceiling(dict$A/2))*n_jump,
                  2)
    E_list <- seq(odd_even_dict$E+(ceiling(dict$E/2)-1)*n_jump,
                  (ceiling(dict$E/2))*n_jump,
                  2)
    B_list <- seq(odd_even_dict$B+(ceiling(dict$B/2)-1)*n_jump,
                  (ceiling(dict$B/2))*n_jump,
                  2)
    F_list <- seq(odd_even_dict$F+(ceiling(dict$F/2)-1)*n_jump,
                  (ceiling(dict$F/2))*n_jump,
                  2)
    C_list <- seq(odd_even_dict$C+(ceiling(dict$C/2)-1)*n_jump,
                  (ceiling(dict$C/2))*n_jump,
                  2)
    G_list <- seq(odd_even_dict$G+(ceiling(dict$G/2)-1)*n_jump,
                  (ceiling(dict$G/2))*n_jump,
                  2)
    D_list <- seq(odd_even_dict$D+(ceiling(dict$D/2)-1)*n_jump,
                  (ceiling(dict$D/2))*n_jump,
                  2)
    H_list <- seq(odd_even_dict$H+(ceiling(dict$H/2)-1)*n_jump,
                  (ceiling(dict$H/2))*n_jump,
                  2)
    
    for (i in 1:length(rows_list)){
      if (rows_list[i] == "A")
        sort_list[i] <- A_list[columns_list[i]]
      else if (rows_list[i] == "B")
        sort_list[i] <- B_list[columns_list[i]]
      else if (rows_list[i] == "C")
        sort_list[i] <- C_list[columns_list[i]]
      else if (rows_list[i] == "D")
        sort_list[i] <- D_list[columns_list[i]]
      else if (rows_list[i] == "E")
        sort_list[i] <- E_list[columns_list[i]]
      else if (rows_list[i] == "F")
        sort_list[i] <- F_list[columns_list[i]]
      else if (rows_list[i] == "G")
        sort_list[i] <- G_list[columns_list[i]]
      else
        sort_list[i] <- H_list[columns_list[i]]
    }
  }
  else if (tolower(sort_order) == "pairrow"){
    n_jump <- num_columns
    
    # Create sort values for rows based on their orders
    A_list <- seq(1+(dict$A-1)*n_jump,
                  (dict$A)*n_jump,
                  1)
    E_list <- seq(1+(dict$E-1)*n_jump,
                  (dict$E)*n_jump,
                  1)
    B_list <- seq(1+(dict$B-1)*n_jump,
                  (dict$B)*n_jump,
                  1)
    F_list <- seq(1+(dict$F-1)*n_jump,
                  (dict$F)*n_jump,
                  1)
    C_list <- seq(1+(dict$C-1)*n_jump,
                  (dict$C)*n_jump,
                  1)
    G_list <- seq(1+(dict$G-1)*n_jump,
                  (dict$G)*n_jump,
                  1)
    D_list <- seq(1+(dict$D-1)*n_jump,
                  (dict$D)*n_jump,
                  1)
    H_list <- seq(1+(dict$H-1)*n_jump,
                  (dict$H)*n_jump,
                  1)
    
    for (i in 1:length(rows_list)){
      if (rows_list[i] == "A")
        sort_list[i] <- A_list[columns_list[i]]
      else if (rows_list[i] == "B")
        sort_list[i] <- B_list[columns_list[i]]
      else if (rows_list[i] == "C")
        sort_list[i] <- C_list[columns_list[i]]
      else if (rows_list[i] == "D")
        sort_list[i] <- D_list[columns_list[i]]
      else if (rows_list[i] == "E")
        sort_list[i] <- E_list[columns_list[i]]
      else if (rows_list[i] == "F")
        sort_list[i] <- F_list[columns_list[i]]
      else if (rows_list[i] == "G")
        sort_list[i] <- G_list[columns_list[i]]
      else
        sort_list[i] <- H_list[columns_list[i]]
    }
  }
  else if (tolower(sort_order) == "pairblock"){
    n_jump <- num_columns*2
    # Create sort values for rows based on their orders
    A_list <- seq(odd_even_dict$A+(ceiling(dict$A/2)-1)*n_jump,
                  (ceiling(dict$A/2))*n_jump,
                  2)
    E_list <- seq(odd_even_dict$E+(ceiling(dict$E/2)-1)*n_jump,
                  (ceiling(dict$E/2))*n_jump,
                  2)
    B_list <- seq(odd_even_dict$B+(ceiling(dict$B/2)-1)*n_jump,
                  (ceiling(dict$B/2))*n_jump,
                  2)
    F_list <- seq(odd_even_dict$F+(ceiling(dict$F/2)-1)*n_jump,
                  (ceiling(dict$F/2))*n_jump,
                  2)
    C_list <- seq(odd_even_dict$C+(ceiling(dict$C/2)-1)*n_jump,
                  (ceiling(dict$C/2))*n_jump,
                  2)
    G_list <- seq(odd_even_dict$G+(ceiling(dict$G/2)-1)*n_jump,
                  (ceiling(dict$G/2))*n_jump,
                  2)
    D_list <- seq(odd_even_dict$D+(ceiling(dict$D/2)-1)*n_jump,
                  (ceiling(dict$D/2))*n_jump,
                  2)
    H_list <- seq(odd_even_dict$H+(ceiling(dict$H/2)-1)*n_jump,
                  (ceiling(dict$H/2))*n_jump,
                  2)
    
    if (odd_even_dict$A == 2)
      A_list <-  A_list - 1
    if (odd_even_dict$E == 2)
      E_list <-  E_list - 1
    if (odd_even_dict$B == 2)
      B_list <-  B_list - 1
    if (odd_even_dict$F == 2)
      F_list <-  F_list - 1
    if (odd_even_dict$C == 2)
      C_list <-  C_list - 1
    if (odd_even_dict$G == 2)
      G_list <-  G_list - 1
    if (odd_even_dict$D == 2)
      D_list <-  D_list - 1
    if (odd_even_dict$G == 2)
      G_list <-  G_list - 1
    for (i in 1:length(rows_list)){
      if (rows_list[i] == "A")
        sort_list[i] <- A_list[columns_list[i]]
      else if (rows_list[i] == "B")
        sort_list[i] <- B_list[columns_list[i]]
      else if (rows_list[i] == "C")
        sort_list[i] <- C_list[columns_list[i]]
      else if (rows_list[i] == "D")
        sort_list[i] <- D_list[columns_list[i]]
      else if (rows_list[i] == "E")
        sort_list[i] <- E_list[columns_list[i]]
      else if (rows_list[i] == "F")
        sort_list[i] <- F_list[columns_list[i]]
      else if (rows_list[i] == "G")
        sort_list[i] <- G_list[columns_list[i]]
      else
        sort_list[i] <- H_list[columns_list[i]]
    }
  }
  sort_list
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

get_baseline_indices2 <- function(well_idx, sample_info, x_vals, y_vals){
  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_avg_list <- NULL
  
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    RU <- y_vals[, i]
    df <- bind_cols("Time" = Time, "RU" = RU)
    colnames(df) <- c("Time", "RU")
    df %>% filter(Time < baseline)  %>% .$RU -> base_meas
    base_meas <- as.numeric(base_meas)
    baseline_avg <- mean(base_meas, na.rm = TRUE)
    #if (is.nan(baseline_avg))
    #  baseline_avg <- 0
    baseline_avg_list <- c(baseline_avg_list, baseline_avg)
  }
  min_baseline <- min(baseline_avg_list)
  if (is.nan(min_baseline)){
    baseline_idx <- start_idx
    min_baseline <- 0
    baseline_avg <- 0
  }
  else
    baseline_idx <- start_idx + which(baseline_avg_list == min_baseline) - 1
  
  #is highest baseline negative?
  if (min_baseline < 0)
    baseline_neg <- TRUE
  else
    baseline_neg <- FALSE
  
  list(BaselineIdx = baseline_idx, MinBaseline = min_baseline, BaselineNegative = baseline_neg)
}

get_baseline_indices <- function(well_idx, sample_info, x_vals, y_vals){
  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1
  
  baseline <- sample_info[well_idx,]$Baseline
  baseline_avg_list <- NULL
  
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    RU <- y_vals[, i]
    df <- bind_cols("Time" = Time, "RU" = RU)
    colnames(df) <- c("Time", "RU")
    df %>% filter(Time < baseline)  %>% .$RU -> base_meas
    baseline_avg <- mean(base_meas, na.rm = TRUE) 
    baseline_avg_list <- c(baseline_avg_list, baseline_avg)
  }
  min_baseline <- min(baseline_avg_list)
  baseline_idx <- start_idx + which(baseline_avg_list == min_baseline) - 1
  
  if (is.nan(min_baseline)){
    baseline_idx <- start_idx
    min_baseline <- 0
    baseline_avg <- 0
  }
  
  #is highest baseline negative?
  if (baseline_avg < 0)
    baseline_neg <- TRUE 
  else
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
  
  map_dfr(.x = tibble(numerical_concentration), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    arrange(numerical_concentration) -> numerical_concentration
  Concentrations <- numerical_concentration
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  colnames(df) <- c("Time", "RU", "Concentration")
  df
}

bulkshift_correction <- function(well_idx, x_vals, y_vals, sample_info,
                                 all_concentrations_ligand){
  
  baseline <- sample_info[well_idx,]$Baseline
  association <- sample_info[well_idx,]$Association
  association <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- all_concentrations_ligand[well_idx]
}

baseline_correction <- function(well_idx, x_vals, y_vals, sample_info){
  
  negative_baseline <- sample_info[well_idx, ]$BaselineNegative
  
  baseline_average <- sample_info[well_idx, ]$BaselineAverage
  start_idx <- sample_info[well_idx, ]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc  #number of concentrations for this well
  end_idx <- start_idx + num_conc - 1
  
  
  #Need to test if baseline of highest conc is < 0
  
  if (sample_info[well_idx,]$Regen. == "N" & !negative_baseline){
    # add baseline average to all timepoints
    
    y_vals[, start_idx:end_idx] <-
      y_vals[, start_idx:end_idx] #- baseline_average
  } else 
  {
    # correct all to mean of zero
    for (i in 1:num_conc){
      # compute baseline average for each concentration
      # subtract from baseline average from all RU vals for that concentration
      
      baseline <- sample_info[well_idx,]$Baseline
      
      if (baseline != 0)
      {
        Time <- x_vals[, start_idx + (i-1)]
        RU <- y_vals[, start_idx + (i-1)]
        
        df <- bind_cols("Time" = Time, "RU" = RU)
        colnames(df) <- c("Time", "RU")
        
        df %>% filter(Time < baseline) %>% .$RU -> base_meas # select for defined baseline time period
        # This command is split up because mean(.$RU) would not parse properly
        base_meas <- as.numeric(base_meas)
        base_corr <- mean(base_meas, na.rm = TRUE)
        if (is.nan(base_corr))
        {
          base_corr <- 0
        }
        ############################################ Need to fix this
        y_vals[, start_idx + (i-1)] <- y_vals[, start_idx + (i-1)] #- base_corr
      }
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
  association <- sample_info[well_idx,]$Association
  
  df <- create_dataframe_with_conc(start_idx, end_idx, x_vals, y_vals,
                                   all_concentrations_values[start_idx:end_idx],
                                   n_time_points)
  
  df %>% filter((Time >= baseline + association - 10)
                & (Time <= baseline + association - 5)) %>%
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
  for (i in start_idx:end_idx){
    end_assoc_resp <- 
      bind_cols(end_assoc_resp, y_vals[which(x_vals[,i] >= association_end)[1], i])
  }
  # record the differences between consecutive responses
  sum_diff <- NULL
  for (i in 1:(num_conc - 1)){
    sum_diff <- bind_cols(sum_diff, end_assoc_resp[1,i+1] - end_assoc_resp[1,i])
  }
  # add sums for each 5 cycle window.
  cum_sum <- rollapply(as_vector(flatten(sum_diff)), 5, FUN = sum)
  start_conc_idx <- which(cum_sum == max(cum_sum))
  
  # check start and end slopes, may be better to fit 4 instead of five
  
  slopes <- as_vector(sum_diff[start_conc_idx:(start_conc_idx+3)])
  
  remove_concentration <- ifelse(slopes < 0.2*mean(slopes), 1, 0)
  
  if(sum(remove_concentration) == 0)
    # return 5 best consecutive concentrations
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

plot_sensorgrams <- function(well_idx, sample_info, x_vals, y_vals, 
                             incl_conc_values, 
                             all_concentrations_values,
                             n_time_points,
                             all = FALSE){
  if (!all){
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
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  df$Concentration <- as_factor(df$Concentration)
  
  sub_title <- paste("Data")
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row", sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  ggplot(df, aes(x = Time, y = RU, color = Concentration)) + geom_point(size = 0.5) +
    guides(color = guide_legend(reverse=TRUE, override.aes = list(size = 5))) +
    ggtitle(ligand_desc, subtitle = sub_title) + 
    geom_vline(xintercept = 420, linetype="dashed", 
               color = "black", size=1) + 
    geom_vline(xintercept = 120, linetype="dashed", 
               color = "grey50", size=1) +
    theme(text = element_text(size=20))
}

plot_sensorgrams_with_nobaseline <- function(well_idx, sample_info, x_vals, y_vals, 
                                             incl_conc_values, 
                                             all_concentrations_values,
                                             n_time_points,
                                             all = FALSE){
  if (!all){
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
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  df$Concentration <- as_factor(df$Concentration)
  
  sub_title <- paste("Bivalent Analyte Data")
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
  #                    sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  ggplot(df, aes(x = Time, y = RU, color = Concentration)) + geom_point(size = 0.5) +
    ggtitle(ligand_desc, subtitle = sub_title) + 
    geom_vline(xintercept = 420, linetype="dashed", color = "black", size=0.5) + 
    xlim(120,2250)
}

combine_output <- function(well_idx, fits_list, plot_list_out, 
                           num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  #R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label)))
  
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
  resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = RU, Concentration = plot_list_out[[well_idx]]$data$Concentration), 
                       aes(x = RU, y = Residuals)) + geom_point(size = 0.01, aes(color = Concentration)) +
    ggtitle(label = "Residuals")
  
  plot_before_baseline <- plot_list_before_baseline[[well_idx]] + ggtitle("Raw Sensorgram")
  grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
               rc_list[[well_idx]], ncol=2)
  #grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
  #             rc_list[[well_idx]], plot_list_before_baseline[[well_idx]], ncol=3)
}

plot_fitting <- function(well_idx, fits_list, plot_list_out, 
                         num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  # par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label, R0_label))) #R0_label "kd2"
  par_names <- as_vector(flatten(c("ka1", "kd1", Rmax_label, R0_label))) #R0_label "kd2"
  
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
  
  # if (info != 5) {            # when info is 5, that means the iterations maxed out and fit is not valid
  ibb <- chol(hessian)
  ih <- chol2inv(ibb)
  p <- length(pars)
  rdf <- length(fit_object$fvec) - p
  resvar <- deviance(fit_object)/rdf
  se <- sqrt(diag(ih) * resvar)
  # }
  
  std_err_full[nonsingular_rows] <- se
  
  df <- data.frame(Estimate = pars, "Std. Error" = std_err_full)
  colnames(df) <- c("Estimate", "Std. Error")
  df
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

plot_bivalent_fitting <- function(well_idx, fits_list, plot_list_out, 
                                  num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  # Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("t0", x)) # For same Rmax
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  "Rmax", R0_label))) #R0_label "kd2"
  # par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2", Rmax_label, R0_label))) #R0_label "kd2"
  
  pars <- coefficients(fits_list[[well_idx]]$R0)
  
  # result_summary <- summary(fits_list[[well_idx]]$R0)
  
  result_summary <- summary_fit_with_constraints(fits_list[[well_idx]]$R0)
  print(well_idx)
  summary_names <- colnames(result_summary)
  
  # rownames(result_summary) <- par_names
  
  result_summary %>% as_tibble -> par_err_table
  
  colnames(par_err_table) <- summary_names
  par_err_table <- bind_cols(Names = par_names, par_err_table)
  
  # par_err_table %>% filter(!str_detect(Names,"R_0")) -> par_err_table
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

plot_monovalent_fitting <- function(well_idx, fits_list, plot_list_out, 
                                    num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  # par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label, R0_label))) #R0_label "kd2"
  par_names <- as_vector(flatten(c("ka1", "kd1", Rmax_label, R0_label))) #R0_label "kd2"
  
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

get_both_tables <- function(well_idx, monovalent_fits_list, bivalent_fits_list, plot_list_out, 
                            num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(monovalent_fits_list[[well_idx]]$error))
    return(NULL)
  if (!is.null(bivalent_fits_list[[well_idx]]$error))
    return(NULL)
  
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  
  # Monovalent
  par_names <- as_vector(flatten(c("ka1", "kd1", Rmax_label, R0_label)))
  
  pars <- coefficients(monovalent_fits_list[[well_idx]]$R0)
  
  result_summary <- summary(monovalent_fits_list[[well_idx]]$R0)
  summary_names <- colnames(result_summary$coefficients)
  result_summary$coefficients %>% as_tibble -> par_err_table
  
  colnames(par_err_table) <- summary_names
  par_err_table <- bind_cols(Names = par_names, par_err_table)
  
  par_err_table %>% filter(!str_detect(Names,"R_0")) -> par_err_table
  par_names <- par_err_table$Names
  par_err_table %>% select(Estimate, `Std. Error`)  %>%
    mutate(Estimate = format(Estimate,big.mark=",",decimal.mark=".")) %>%
    tableGrob(rows = par_names, theme = ttheme_minimal()) -> tb1
  
  # title <- textGrob("Monovalent", gp = gpar(fontsize = 10))
  # 
  # tb1 <- gtable_add_grob(
  #   tb1, list(title),
  #   t=10, l=10, r =10
  # )
  
  # Bivalent
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label, R0_label))) #R0_label "kd2"
  
  pars <- coefficients(bivalent_fits_list[[well_idx]]$R0)
  
  result_summary <- summary(bivalent_fits_list[[well_idx]]$R0)
  summary_names <- colnames(result_summary$coefficients)
  result_summary$coefficients %>% as_tibble -> par_err_table
  
  colnames(par_err_table) <- summary_names
  par_err_table <- bind_cols(Names = par_names, par_err_table)
  
  par_err_table %>% filter(!str_detect(Names,"R_0")) -> par_err_table
  par_names <- par_err_table$Names
  par_err_table %>% select(Estimate, `Std. Error`)  %>%
    mutate(Estimate = format(Estimate,big.mark=",",decimal.mark=".")) %>%
    tableGrob(rows = par_names, theme = ttheme_minimal()) -> tb2
  
  # title <- textGrob("Bivalent", gp = gpar(fontsize = 10))
  # 
  # tb2 <- gtable_add_grob(
  #   tb2, list(title),
  #   pos = -1
  # )
  
  grid.arrange(tb1, tb2, ncol=2)
}



plot_data <- function(well_idx, fits_list, plot_list_out, 
                      num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  #R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  # par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label)))
  par_names <- as_vector(flatten(c("ka1", "kd1", Rmax_label)))
  
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
  resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = RU, Concentration = plot_list_out[[well_idx]]$data$Concentration), 
                       aes(x = RU, y = Residuals)) + geom_point(size = 0.01, aes(color = Concentration)) +
    ggtitle(label = "Residuals")
  
  plot_before_baseline <- plot_list_before_baseline[[well_idx]] + ggtitle("Raw Sensorgram")
  grid.arrange(plot_list_before_baseline[[well_idx]])
  #grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
  #             rc_list[[well_idx]], plot_list_before_baseline[[well_idx]], ncol=3)
}

plot_data_highest_conc <- function(well_idx, fits_list, plot_list_out, 
                                   num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- num_conc[well_idx]
  Rmax_label <- map_dfr(tibble(1:num_conc), function(x) paste("Rmax", x))
  #R0_label <- map_dfr(tibble(1:num_conc), function(x) paste("R_0", x))
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label)))
  
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
  resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = RU, Concentration = plot_list_out[[well_idx]]$data$Concentration), 
                       aes(x = RU, y = Residuals)) + geom_point(size = 0.01, aes(color = Concentration)) +
    ggtitle(label = "Residuals")
  
  plot_before_baseline <- plot_list_before_baseline[[well_idx]] + ggtitle("Raw Sensorgram")
  grid.arrange(plot_list_before_baseline[[well_idx]])
  #grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
  #             rc_list[[well_idx]], plot_list_before_baseline[[well_idx]], ncol=3)
}

plot_sensorgrams_highest_conc <- function(well_idx, sample_info, x_vals, y_vals, 
                                          incl_conc_values, 
                                          all_concentrations_values,
                                          n_time_points,
                                          all = FALSE){
  if (!all){
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
  baseline <- sample_info[well_idx,]$Baseline
  association <- sample_info[well_idx,]$Association
  dissociation <- sample_info[well_idx,]$Dissociation
  
  
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
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  
  colnames(df) <- c("Time", "RU", "Concentration")
  df %>% filter(Time > baseline & Time < baseline + association + dissociation) -> df
  df <- na.omit(df)
  df$Concentration <- as_factor(df$Concentration)
  
  df %>% filter(Concentration == all_concentrations_values[length(all_concentrations_values)]) -> df_highest
  
  sub_title <- paste("Bivalent Analyte Data")
  # sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row", sample_info[well_idx,]$Row,
  #                    "Column", sample_info[well_idx,]$Column)
  
  ggplot(df_highest, aes(x = Time, y = RU, color = Concentration)) + geom_point(size = 0.5) +
    #geom_line(data =  fit_df, aes(x = Time, y = RU, group = Concentration), color = "black", size=1) +
    ggtitle(ligand_desc, subtitle = sub_title) #+ #geom_vline(xintercept = 420, linetype="dashed", 
  #           color = "black", size=1) + xlim(420,730) + ylim(1, 1.4)
}

plot_fitting_X1X2 <- function(well_idx, fits_list, plot_list_out, 
                              num_conc, rc_list, plot_list_before_baseline){
  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)
  num_conc <- 1:num_conc
  Rmax_label <- map_dfr(tibble(1:1), function(x) paste("Rmax", x))
  R0_label <- map_dfr(tibble(1:1), function(x) paste("R_0", x))
  par_names <- as_vector(flatten(c("ka1", "ka2", "kd1", "kd2",  Rmax_label)))
  
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
  fits_list[[well_idx]]$FitOutcomes$X1 -> X1
  fits_list[[well_idx]]$FitOutcomes$X2 -> X2
  
  min_Time <- min(Time_resid)
  max_Time <- max(Time_resid)
  
  RUdata <- plot_list_before_baseline[[well_idx]]$data$RU
  Timedata <- plot_list_before_baseline[[well_idx]]$data$Time
  resid_plot <- ggplot(data = tibble(Residuals = RU_resid, Time = RU),
                       aes(x = RU, y = Residuals)) + geom_point(size = 0.01, color='black') +
    ggtitle(label = "Residuals")
  
  plot_before_baseline <- plot_list_before_baseline[[well_idx]] + ggtitle("Raw Sensorgram")
  
  df <- tibble(RU = RU, Time = Time_resid, X1 = X1, X2 = X2)
  df_data <- tibble(data = RUdata, Time = Timedata)
  df_data %>% filter(Time > min_Time & Time < max_Time) -> df_data
  highest_conc <- max(as.numeric(as.character((unique(plot_list_out[[well_idx]]$data$Concentration)))), na.rm=TRUE)
  #df %>% filter(Concentration == all_concentrations_values[length(all_concentrations_values)]) -> df_highest
  df$Concentration <- highest_conc
  sim_plot <- ggplot() +
    geom_point(data = df_data, aes(x = Time, y = data, colour = "data"), size=0.5) +
    geom_line(data = df, aes(x = Time, y = RU, colour = "RU"), size=1) + 
    geom_line(data = df, aes(x = Time, y = X1, colour = "AL1"), size=1) + 
    geom_line(data = df, aes(x = Time, y = X2, colour = "AL2"), size=1) + 
    geom_vline(xintercept = 420, linetype="dashed", color = "black", size=1) +
    ggtitle(label = "Simulation for Highest Analyte Concentration")
  
  # grid.arrange(sim_plot)
  grid.arrange(sim_plot, tb1, resid_plot,
               rc_list[[well_idx]], ncol=2)
}


print_output <- function(well_idx, pages_list, plot_list_out, sample_info){
  
  if (is.null(pages_list[[well_idx]]$error) & !is.null(pages_list[[well_idx]]))
    return(arrangeGrob(pages_list[[well_idx]]))
  
  print(well_idx)
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
  
  df <- bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations)
  
  colnames(df) <- c("Time", "RU", "Concentration")
  
  
  df %>% filter((Time >= baseline + association - 10) 
                & (Time <= baseline + association - 5)) %>% 
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