library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_3")


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
sample_time <- seq(0,672,1)
admin.type <-"oral"
admin.dose <- 25 * BW*1000 #ug
admin.time <- 0


sex <- "F"
#Female, oral 25mg/kg dose
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

score <- rep(NA, 4)
df_excreta <- openxlsx::read.xlsx("Raw_Data/Kemper_2003/Kemper_Excreta_Loccisano.xlsx")
obs_F_25_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                                     df_excreta$Tissue == "Urine" & df_excreta$Sex == "F",]
rounded_time <- round(obs_F_25_urine$Time_h)
rounded_soltime <- round(solution_F_oral_25$time)
preds_urine_F_oral_25 <- solution_F_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[1] <- AAFE(preds_urine_F_oral_25, (obs_F_25_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_F_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_F_25_urine$Dose_mg_per_kg,
                          "Tissue" = obs_F_25_urine$Tissue ,
                          "Type" = obs_F_25_urine$Type, sex = obs_F_25_urine$Sex,
                          "Observed" = (obs_F_25_urine$`Cum_dose_%`/100)*admin.dose/1000,
                          "Predicted" = preds_urine_F_oral_25,  "Time" = obs_F_25_urine$Time_h)


obs_F_25_feces <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                               df_excreta$Tissue == "Feces" & df_excreta$Sex == "F",]
rounded_time <- round(obs_F_25_feces$Time_h)
rounded_soltime <- round(solution_F_oral_25$time)
preds_feces_F_oral_25 <- solution_F_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[2] <- AAFE(preds_feces_F_oral_25, (obs_F_25_feces$`Cum_dose_%`/100)*admin.dose/1000)
results_feces_F_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_F_25_feces$Dose_mg_per_kg,
                                     "Tissue" = obs_F_25_feces$Tissue ,
                                     "Type" = obs_F_25_feces$Type, sex = obs_F_25_feces$Sex,
                                     "Observed" = (obs_F_25_feces$`Cum_dose_%`/100)*admin.dose/1000,
                                     "Predicted" = preds_feces_F_oral_25, "Time" = obs_F_25_feces$Time_h)



obs_M_25_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                               df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_25_urine$Time_h)
rounded_soltime <- round(solution_M_oral_25$time)
preds_urine_M_oral_25 <- solution_M_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[3] <- AAFE(preds_urine_M_oral_25,(obs_M_25_urine$`Cum_dose_%`/100)*admin.dose/1000)
results_urine_M_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_M_25_urine$Dose_mg_per_kg,
                                     "Tissue" = obs_M_25_urine$Tissue ,
                                     "Type" = obs_M_25_urine$Type, sex = obs_M_25_urine$Sex,
                                     "Observed" = (obs_M_25_urine$`Cum_dose_%`/100)*admin.dose/1000,
                                     "Predicted" = preds_urine_M_oral_25, "Time" = obs_M_25_urine$Time_h)

obs_M_25_feces <- df_excreta[df_excreta$Dose_mg_per_kg == 25 & 
                               df_excreta$Tissue == "Feces" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_25_feces$Time_h)
rounded_soltime <- round(solution_M_oral_25$time)
preds_feces_M_oral_25 <- solution_M_oral_25[rounded_soltime %in% rounded_time, "Murine"]/1000
score[4] <- AAFE(preds_feces_M_oral_25,(obs_M_25_feces$`Cum_dose_%`/100)*admin.dose/1000)
results_feces_M_oral_25<- data.frame("Study" = "Kemper_2003", "Dose" =  obs_M_25_feces$Dose_mg_per_kg,
                                     "Tissue" = obs_M_25_feces$Tissue ,
                                     "Type" = obs_M_25_feces$Type, sex = obs_M_25_feces$Sex,
                                     "Observed" = (obs_M_25_feces$`Cum_dose_%`/100)*admin.dose/1000,
                                     "Predicted" = preds_feces_M_oral_25, "Time" = obs_M_25_feces$Time_h)

results_df <- rbind(results_urine_F_oral_25, results_feces_F_oral_25, results_urine_M_oral_25, results_feces_M_oral_25)


AAFE_Loccisano <- mean(score)
print(paste0("The AAFE on the excreta data of Kemper et al. (2003) from Loccisano is ", AAFE_Loccisano))


write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Kemper_2003_excreta_Loccisano_results.csv",
          "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_3/Kemper_2003_excreta_Loccisano_results.csv",
          row.names =F)
