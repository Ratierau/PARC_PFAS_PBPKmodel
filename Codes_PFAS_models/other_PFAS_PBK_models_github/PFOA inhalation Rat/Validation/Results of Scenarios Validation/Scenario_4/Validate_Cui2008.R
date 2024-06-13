library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_4")

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
BW <- 0.2 #kg
sex <- "M"
admin.type <-"oral"
#5mg/kg dose
admin.time <- seq(0.01,27*24.01,24)
dose <- rep(5,28) #mg/kg
admin.dose <- dose * BW*1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <- seq(0,28*24,1)
solution_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

# Carefully select the order so that it matches the experimental data presented in the xlsx file
pred_comps <- c( "Cart", "Cliver", "Ckidney","Clungs","Cheart",
                  "Cspleen", "Cgonads", "Cbrain" )
solution_5 <- solution_5[solution_5$time == 28*24,pred_comps]/1000 #[ug/L]/1000-->[ug/g]


#20mg/kg dose
dose <- rep(20,28) #mg/kg
admin.dose <- dose * BW * 1000 #ug
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_20 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                       y = inits, parms = parameters, events = events,
                                       method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_20 <- solution_20[solution_20$time == 28*24,pred_comps]/1000 #[ug/L]/1000-->[ug/g]

##############################
### Load experimental data ###
##############################

# Load the mean values 
df <- openxlsx::read.xlsx("Raw_Data/Cui_2008/Cui_tissues.xlsx")
tissues_5 <- df[df$Dose_mg_per_kg == 5,]
tissues_20 <- df[df$Dose_mg_per_kg == 20,]

score <- rep(NA,2)
# Gather the required data for the x-y plot
score[1] <- AAFE(solution_5,tissues_5$`Concentration_ug/g`)
results_df_5<- data.frame("Study" = "Cui_2010", "Dose" =  tissues_5$Dose_mg_per_kg,
                          "Tissue" = tissues_5$Tissue ,
                          "Type" = "oral",
                          "Observed" =tissues_5$`Concentration_ug/g`,
                          "Predicted" = unname(t(solution_5)), "Time" = tissues_5$Time_h )


results_df_20<- data.frame("Study" = "Cui_2010", "Dose" =  tissues_20$Dose_mg_per_kg,
                           "Tissue" = tissues_20$Tissue,
                           "Type" = "oral",
                           "Observed" = tissues_20$`Concentration_ug/g`,
                           "Predicted" = unname(t(solution_20)), "Time" = tissues_20$Time_h)

score[2] <- AAFE(solution_20,tissues_20$`Concentration_ug/g`)

results_df <- rbind(results_df_5, results_df_20)

AAFE_Cui <- mean(score)
print(paste0("The AAFE on the tissue data of Cui et al. (2010) was ", AAFE_Cui))

write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Cui_2008_results.csv",
          "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_4/Cui_2008_results.csv",
          row.names =F)
