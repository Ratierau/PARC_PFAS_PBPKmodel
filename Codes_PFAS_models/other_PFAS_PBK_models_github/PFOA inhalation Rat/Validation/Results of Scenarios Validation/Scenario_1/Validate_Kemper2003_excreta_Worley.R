library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_1")

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

#===============
# Generate predictions
#===============


#load("model.RData")
# Body weight 
BW_F <- 0.2 #kg
BW_M <- 0.3 #kg
sample_time <- seq(0,673,1)
admin.type <-"oral"
admin.time <- 0
sex <- "F"
#Female, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_F_oral_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

#Female, oral 5mg/kg dose
admin.dose <- 5 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_F_oral_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Female, oral 25mg/kg dose
admin.dose <- 25 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_F,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_F_oral_25 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))


sex <- "M"
#Male, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_M_oral_1 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Male, oral 5mg/kg dose
admin.dose <- 5 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_M_oral_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))

#Male, oral 25mg/kg dose
admin.dose <- 25 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW_M,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_M_oral_25 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                             y = inits, parms = parameters, events = events,
                                             method="lsodes",rtol = 1e-7, atol = 1e-7))


##############################
### Load experimental data ###
##############################
score <- rep(NA, 6)
df_excreta <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Excreta_Worley.xlsx")

admin.dose <- 1 * BW * 1000 #ug
obs_F_1_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 1 & 
                              df_excreta$Tissue == "Urine" & df_excreta$Sex == "F",]
rounded_time <- round(obs_F_1_urine$Time_h)
rounded_soltime <- round(solution_F_oral_1$time)
preds_urine_F_oral_1 <- solution_F_oral_1[rounded_soltime %in% rounded_time, "Murine"]/1000
score[1] <- AAFE(preds_urine_F_oral_1, (obs_F_1_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_F_oral_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_F_1_urine$Dose_mg_per_kg,
                                    "Tissue" = obs_F_1_urine$Tissue ,
                                    "Type" = obs_F_1_urine$Type, sex = obs_F_1_urine$Sex,
                                    "Observed" = (obs_F_1_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                    "Predicted" = preds_urine_F_oral_1,  "Time" = obs_F_1_urine$Time_h)

admin.dose <- 5 * BW * 1000 #ug
obs_F_5_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 5 & 
                              df_excreta$Tissue == "Urine" & df_excreta$Sex == "F",]
rounded_time <- round(obs_F_5_urine$Time_h)
rounded_soltime <- round(solution_F_oral_5$time)
preds_urine_F_oral_5 <- solution_F_oral_5[rounded_soltime %in% rounded_time, "Murine"]/1000
score[2] <- AAFE(preds_urine_F_oral_5, (obs_F_5_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_F_oral_5<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_F_5_urine$Dose_mg_per_kg,
                                    "Tissue" = obs_F_5_urine$Tissue ,
                                    "Type" = obs_F_5_urine$Type, sex = obs_F_5_urine$Sex,
                                    "Observed" = (obs_F_5_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                    "Predicted" = preds_urine_F_oral_5,  "Time" = obs_F_5_urine$Time_h)


admin.dose <- 25 * BW * 1000 #ug
obs_F_25_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                                     df_excreta$Tissue == "Urine" & df_excreta$Sex == "F",]
rounded_time <- round(obs_F_25_urine$Time_h)
rounded_soltime <- round(solution_F_oral_25$time)
preds_urine_F_oral_25 <- solution_F_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[3] <- AAFE(preds_urine_F_oral_25, (obs_F_25_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_F_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_F_25_urine$Dose_mg_per_kg,
                          "Tissue" = obs_F_25_urine$Tissue ,
                          "Type" = obs_F_25_urine$Type, sex =obs_F_25_urine$Sex,
                          "Observed" = (obs_F_25_urine$`Cum_dose_%`/100)*admin.dose/1000,
                          "Predicted" = preds_urine_F_oral_25,  "Time" = obs_F_25_urine$Time_h)


admin.dose <- 1 * BW * 1000 #ug
obs_M_1_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 1 & 
                              df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_1_urine$Time_h)
rounded_soltime <- round(solution_M_oral_1$time)
preds_urine_M_oral_1 <- solution_M_oral_1[rounded_soltime %in% rounded_time, "Murine"]/1000
score[4] <- AAFE(preds_urine_M_oral_1, (obs_M_1_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_M_oral_1<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_M_1_urine$Dose_mg_per_kg,
                                    "Tissue" = obs_M_1_urine$Tissue ,
                                    "Type" = obs_M_1_urine$Type, sex = obs_M_1_urine$Sex,
                                    "Observed" = (obs_M_1_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                    "Predicted" = preds_urine_M_oral_1,  "Time" = obs_M_1_urine$Time_h)

admin.dose <- 5 * BW * 1000 #ug
obs_M_5_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 5 & 
                              df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_5_urine$Time_h)
rounded_soltime <- round(solution_M_oral_5$time)
preds_urine_M_oral_5 <- solution_M_oral_5[rounded_soltime %in% rounded_time, "Murine"]/1000
score[5] <- AAFE(preds_urine_M_oral_5, (obs_M_5_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_M_oral_5<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_M_5_urine$Dose_mg_per_kg,
                                    "Tissue" = obs_M_5_urine$Tissue ,
                                    "Type" = obs_M_5_urine$Type, sex = obs_M_5_urine$Sex,
                                    "Observed" = (obs_M_5_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                    "Predicted" = preds_urine_M_oral_5,  "Time" = obs_M_5_urine$Time_h)

admin.dose <- 25 * BW * 1000 #ug
obs_M_25_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                               df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_25_urine$Time_h)
rounded_soltime <- round(solution_M_oral_25$time)
preds_urine_M_oral_25 <- solution_M_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[6] <- AAFE(preds_urine_M_oral_25, (obs_M_25_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_M_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_M_25_urine$Dose_mg_per_kg,
                                     "Tissue" = obs_M_25_urine$Tissue ,
                                     "Type" = obs_M_25_urine$Type, sex = obs_M_25_urine$Sex,
                                     "Observed" = (obs_M_25_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                     "Predicted" = preds_urine_M_oral_25,  "Time" = obs_M_25_urine$Time_h)



results_df <- rbind(results_urine_F_oral_1, results_urine_F_oral_5, 
                    results_urine_F_oral_25, results_urine_M_oral_1,
                    results_urine_M_oral_5,  results_urine_M_oral_25)


AAFE_Worley <- mean(score)
print(paste0("The AAFE on the excreta data of Kemper et al. (2003) from Worley is ", AAFE_Worley))



write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Kemper_2003_excreta_Worley_results.csv",
          "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_1/Kemper_2003_excreta_Worley_results.csv",
          row.names =F)
