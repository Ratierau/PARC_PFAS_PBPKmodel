library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_3")

#===============
# Generate predictions
#===============

#  absolute average fold error
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}


#load("model.RData")
# Body weight 
BW_F <- 0.2 #kg
BW_M <- 0.3 #kg
admin.time <- 0


sex <- "F"
#Female, IV, 1mg/kg dose
admin.type <-"iv"
admin.dose <- 1 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1), seq(1,72,1))
solution_F_iv_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

#Female, oral 1mg/kg dose
admin.type <-"oral"
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,100,1))
solution_F_oral_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

#Female, oral 5mg/kg dose
admin.dose <- 5 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,72,1))
solution_F_oral_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Female, oral 25mg/kg dose
admin.dose <- 25 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,96,1))
solution_F_oral_25 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))


sex <- "M"
#Male, IV, 1mg/kg dose
admin.type <-"iv"
admin.time <- 0
admin.dose <- 1 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,528,1))
solution_M_iv_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                           y = inits, parms = parameters, events = events,
                                           method="lsodes",rtol = 1e-7, atol = 1e-7))

#Male, oral 1mg/kg dose
admin.type <-"oral"
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,529,1))
solution_M_oral_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Male, oral 5mg/kg dose
admin.dose <- 5 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,527,1))
solution_M_oral_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Male, oral 25mg/kg dose
admin.dose <- 25 * BW*1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <-  c( seq(0,0.9,0.1),seq(1,528,1))
solution_M_oral_25 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))


##############################
### Load experimental data ###
##############################

# Custom rounding function
custom_round <- function(time) {
  if (time < 1) {
    return(round(time, 1))
  } else {
    return(round(time, 0))
  }
}

score <- rep(NA,8)
df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_F_iv_1 <- df_serum[df_serum$Dose_mg_per_kg == 1 & df_serum$Type == "iv" & df_serum$Sex == "F",]
rounded_time <- sapply(obs_serum_F_iv_1$Time_h, custom_round)
rounded_soltime <- sapply(solution_F_iv_1$time, custom_round)
preds_serum_F_iv_1 <- solution_F_iv_1[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[1] <- AAFE(preds_serum_F_iv_1, obs_serum_F_iv_1$`Concentration_ug/g`)
results_serum_F_iv_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_F_iv_1$Dose_mg_per_kg,
                          "Tissue" = obs_serum_F_iv_1$Tissue ,
                          "Type" = "iv", sex = "F",
                          "Observed" = obs_serum_F_iv_1$`Concentration_ug/g`,
                          "Predicted" = preds_serum_F_iv_1, "Time" = obs_serum_F_iv_1$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_F_oral_1 <- df_serum[df_serum$Dose_mg_per_kg == 1 & df_serum$Type == "oral" & df_serum$Sex == "F",]
rounded_time <- sapply(obs_serum_F_oral_1$Time_h, custom_round)
rounded_soltime <- sapply(solution_F_oral_1$time, custom_round)
preds_serum_F_oral_1 <- solution_F_oral_1[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[2] <- AAFE(preds_serum_F_oral_1, obs_serum_F_oral_1$`Concentration_ug/g`)
results_serum_F_oral_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_F_oral_1$Dose_mg_per_kg,
                                    "Tissue" = obs_serum_F_oral_1$Tissue ,
                                    "Type" = "oral", sex = "F",
                                    "Observed" = obs_serum_F_oral_1$`Concentration_ug/g`,
                                    "Predicted" = preds_serum_F_oral_1, "Time" = obs_serum_F_oral_1$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_F_oral_5 <- df_serum[df_serum$Dose_mg_per_kg == 5 & df_serum$Type == "oral" & df_serum$Sex == "F",]
rounded_time <- sapply(obs_serum_F_oral_5$Time_h, custom_round)
rounded_soltime <- sapply(solution_F_oral_5$time, custom_round)
preds_serum_F_oral_5 <- solution_F_oral_5[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[3] <- AAFE(preds_serum_F_oral_5, obs_serum_F_oral_5$`Concentration_ug/g`)
results_serum_F_oral_5<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_F_oral_5$Dose_mg_per_kg,
                                    "Tissue" = obs_serum_F_oral_5$Tissue ,
                                    "Type" = "oral", sex = "F",
                                    "Observed" = obs_serum_F_oral_5$`Concentration_ug/g`,
                                    "Predicted" = preds_serum_F_oral_5, "Time" = obs_serum_F_oral_5$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_F_oral_25 <- df_serum[df_serum$Dose_mg_per_kg == 25 & df_serum$Type == "oral" & df_serum$Sex == "F",]
rounded_time <- sapply(obs_serum_F_oral_25$Time_h, custom_round)
rounded_soltime <- sapply(solution_F_oral_25$time, custom_round)
preds_serum_F_oral_25 <- solution_F_oral_25[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[4] <- AAFE(preds_serum_F_oral_25, obs_serum_F_oral_25$`Concentration_ug/g`)
results_serum_F_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_F_oral_25$Dose_mg_per_kg,
                                     "Tissue" = obs_serum_F_oral_25$Tissue ,
                                     "Type" = "oral", sex = "F",
                                     "Observed" = obs_serum_F_oral_25$`Concentration_ug/g`,
                                     "Predicted" = preds_serum_F_oral_25, "Time" = obs_serum_F_oral_25$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_M_iv_1 <- df_serum[df_serum$Dose_mg_per_kg == 1 & df_serum$Type == "iv" & df_serum$Sex == "M",]
rounded_time <- sapply(obs_serum_M_iv_1$Time_h, custom_round)
rounded_soltime <- sapply(solution_M_iv_1$time, custom_round)
preds_serum_M_iv_1 <- solution_M_iv_1[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[5] <- AAFE(preds_serum_M_iv_1, obs_serum_M_iv_1$`Concentration_ug/g`)
results_serum_M_iv_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_M_iv_1$Dose_mg_per_kg,
                                  "Tissue" = obs_serum_M_iv_1$Tissue ,
                                  "Type" = "iv", sex = "M",
                                  "Observed" = obs_serum_M_iv_1$`Concentration_ug/g`,
                                  "Predicted" = preds_serum_M_iv_1, "Time" = obs_serum_M_iv_1$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_M_oral_1 <- df_serum[df_serum$Dose_mg_per_kg == 1 & df_serum$Type == "oral" & df_serum$Sex == "M",]
rounded_time <- sapply(obs_serum_M_oral_1$Time_h, custom_round)
rounded_soltime <- sapply(solution_M_oral_1$time, custom_round)
preds_serum_M_oral_1 <- solution_M_oral_1[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[6] <- AAFE(preds_serum_M_oral_1, obs_serum_M_oral_1$`Concentration_ug/g`)
results_serum_M_oral_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_M_oral_1$Dose_mg_per_kg,
                                    "Tissue" = obs_serum_M_oral_1$Tissue ,
                                    "Type" = "oral", sex = "M",
                                    "Observed" = obs_serum_M_oral_1$`Concentration_ug/g`,
                                    "Predicted" = preds_serum_M_oral_1, "Time" = obs_serum_M_oral_1$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_M_oral_5 <- df_serum[df_serum$Dose_mg_per_kg == 5 & df_serum$Type == "oral" & df_serum$Sex == "M",]
rounded_time <- sapply(obs_serum_M_oral_5$Time_h, custom_round)
rounded_soltime <- sapply(solution_M_oral_5$time, custom_round)
preds_serum_M_oral_5 <- solution_M_oral_5[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[7] <- AAFE(preds_serum_M_oral_5, obs_serum_M_oral_5$`Concentration_ug/g`)
results_serum_M_oral_5<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_M_oral_5$Dose_mg_per_kg,
                                    "Tissue" = obs_serum_M_oral_5$Tissue ,
                                    "Type" = "oral", sex = "M",
                                    "Observed" = obs_serum_M_oral_5$`Concentration_ug/g`,
                                    "Predicted" = preds_serum_M_oral_5, "Time" = obs_serum_M_oral_5$Time_h)


df_serum <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Serum.xlsx")
obs_serum_M_oral_25 <- df_serum[df_serum$Dose_mg_per_kg == 25 & df_serum$Type == "oral" & df_serum$Sex == "M",]
rounded_time <- sapply(obs_serum_M_oral_25$Time_h, custom_round)
rounded_soltime <- sapply(solution_M_oral_25$time, custom_round)
preds_serum_M_oral_25 <- solution_M_oral_25[rounded_soltime %in% rounded_time, "Cplasma"]/1000
score[8] <- AAFE(preds_serum_M_oral_25, obs_serum_M_oral_25$`Concentration_ug/g`)
results_serum_M_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_serum_M_oral_25$Dose_mg_per_kg,
                                     "Tissue" = obs_serum_M_oral_25$Tissue ,
                                     "Type" = "oral", sex = "M",
                                     "Observed" = obs_serum_M_oral_25$`Concentration_ug/g`,
                                     "Predicted" = preds_serum_M_oral_25, "Time" = obs_serum_M_oral_25$Time_h)

results_df <- rbind(results_serum_F_iv_1, results_serum_F_oral_1, results_serum_F_oral_5, results_serum_F_oral_25,
                    results_serum_M_iv_1, results_serum_M_oral_1, results_serum_M_oral_5, results_serum_M_oral_25)
AAFE_Kemper_serum <- mean(score)
print(paste0("The AAFE on the serum data of Kemper et al. (2003) was ", AAFE_Kemper_serum))
write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Kemper_2003_serum_results.csv",
          "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_3/Kemper_2003_serum_results.csv",
          row.names =F)
