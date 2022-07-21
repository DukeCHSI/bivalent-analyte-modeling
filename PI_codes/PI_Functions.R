biEpitopicLigandHeterogenousAnalyte_model <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dL <- - ka1*A1m*L + kd1*A1L - ka2*A2m*L + kd2*A2L
    dA1L <- ka1*A1m*L - kd1*A1L - ka4*A1L*A2m + kd4*A1LA2
    dA2L <- ka2*A2m*L - kd2*A2L - ka3*A2L*A1m + kd3*A1LA2
    dA1LA2 <- ka3*A2L*A1m - kd3*A1LA2 + ka4*A1L*A2m - kd4*A1LA2
    
    #return the rate of change
    # print(list(c(t, L, A1L, A2L, A1LA2)))
    # print(dL + dA1L + dA2L + dA1LA2)
    list(c(dL, dA1L, dA2L, dA1LA2))
  })
}

monovalent_model <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dL <- -ka1*Am*L + kd1*X1
    dX1 <- ka1*Am*L - kd1*X1
    
    # return the rate of change
    # print(list(c(t, L, X1)))
    list(c(dL, dX1))
  })
}

objective_function <- function(varied_pars, 
                               fixed_idx, fixed_param, 
                               nonPI_pars, 
                               data, t_asc, t_dis, 
                               num_conc, incl_concentrations1, incl_concentrations2, 
                               param_names, 
                               use_regeneration, 
                               model, 
                               use_globalRmax,
                               use_secondRebind){
  temp_pars <- append(varied_pars, fixed_param, after=fixed_idx-1)
  
  pars <- append(nonPI_pars, temp_pars, after=4)
  # print(pars)
  sim_results <- run_model(pars, t_asc, t_dis, num_conc, incl_concentrations1, incl_concentrations2, param_names, 
               use_regeneration, 
               model, 
               use_globalRmax,
               use_secondRebind)
  
  sim_RU <- sim_results$RU
  data_RU <- data$RU
  
  error <- sim_RU - data_RU
  
  error
}
run_model <- function(pars, t_asc, t_dis, num_conc, incl_concentrations1, incl_concentrations2, param_names, 
                      use_regeneration, 
                      model, 
                      use_globalRmax,
                      use_secondRebind){
  
  df_RU <- NULL
  df_RU1 <- NULL
  df_RU2 <- NULL
  df_RU12 <- NULL
  df_Time <- NULL
  df_Concentration1 <- NULL
  df_Concentration2 <- NULL
  
  for (i in 1:num_conc){
    
    Time <- c(t_asc, t_dis)
    Concentration1 <- incl_concentrations1
    Concentration2 <- incl_concentrations2
    
    if (model == 'biEpitopicLigandHeterogenousAnalyte'){
      
      ka1 <- pars[param_names == "ka1"]
      ka2 <- pars[param_names == "ka2"]
      ka3 <- pars[param_names == "ka3"]
      ka4 <- pars[param_names == "ka4"]
      
      kd1 <- pars[param_names == "kd1"]
      kd2 <- pars[param_names == "kd2"]
      kd3 <- pars[param_names == "kd3"]
      kd4 <- pars[param_names == "kd4"]
      
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
                      ka3 = ka3,
                      kd3 = kd3,
                      ka4 = ka4,
                      kd4 = kd4,
                      A1m = Concentration1[i],
                      A2m = Concentration2[i]) ############## Need to change this to reflect the info file
      
      t_asc <- t_asc
      out <- ode(y = state, times = t_asc, func = biEpitopicLigandHeterogenousAnalyte_model, parms = parameters)
      df_assoc_RU <- out[,3] + out[,4] + 2*out[,5] +  RI
      df_assoc_RU1 <- out[,3]
      df_assoc_RU2 <- out[,4]
      df_assoc_RU12 <- 2*out[,5]
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], t_dis)
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      ka2 = ka2,
                      kd2 = kd2,
                      ka3 = ka3,
                      kd3 = kd3,
                      ka4 = ka4,
                      kd4 = kd4,
                      A1m = 0,
                      A2m = 0)
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3],
                 out[length(t_asc),4],
                 out[length(t_asc),5])
      
      out <- ode(y = state, times = t_dis, func = biEpitopicLigandHeterogenousAnalyte_model, parms = parameters)
      df_dissoc_RU <- out[2:length(t_dis),3] + out[2:length(t_dis),4] + 2*out[2:length(t_dis),5] +  RI
      df_dissoc_RU1 <- out[2:length(t_dis),3]
      df_dissoc_RU2 <- out[2:length(t_dis),4]
      df_dissoc_RU12 <- 2*out[2:length(t_dis),5]
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
                      Am = Concentration1[i])
      
      t_asc <- t_asc
      out <- ode(y = state, times = t_asc, func = monovalent_model, parms = parameters)
      df_assoc_RU <- out[,3] + RI
      
      # Dissociation
      t_dis <- c(t_asc[length(t_asc)], t_dis)
      
      parameters <- c(ka1 = ka1,
                      kd1 = kd1,
                      Am = 0)
      
      state <- c(out[length(t_asc),2],
                 out[length(t_asc),3])
      
      out <- ode(y = state, times = t_dis, func = monovalent_model, parms = parameters)
      df_dissoc_RU <- out[2:length(t_dis),3] + RI
    }
    
    
    
    df_RU <- c(df_RU, df_assoc_RU, df_dissoc_RU)
    df_RU1 <- c(df_RU1, df_assoc_RU1, df_dissoc_RU1)
    df_RU2 <- c(df_RU2, df_assoc_RU2, df_dissoc_RU2)
    df_RU12 <- c(df_RU12, df_assoc_RU12, df_dissoc_RU12)
    df_Time <- c(df_Time, Time)
    df_Concentration1 <- c(df_Concentration1, rep(Concentration1[i], length(Time)))
    df_Concentration2 <- c(df_Concentration2, rep(Concentration2[i], length(Time)))
    
    full_output_RU <- t(rbind("Time"=df_Time, "RU"= df_RU, "RU1"= df_RU1, 
                              "RU2" = df_RU2, "RU12" = df_RU12, "Concentration1"=df_Concentration1, "Concentration2"=df_Concentration2))
  }
  as_tibble(full_output_RU)
}