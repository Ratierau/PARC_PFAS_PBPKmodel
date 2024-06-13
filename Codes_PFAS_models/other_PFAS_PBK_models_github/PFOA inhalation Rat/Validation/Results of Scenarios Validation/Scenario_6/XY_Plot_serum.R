library(ggplot2)

#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_6")

# Load Results 

# 
Kemper_serum_results <- read.csv("Kemper_2003_serum_results.csv", header = T)
Kemper_iv_F_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "iv" &
                               Kemper_serum_results$sex == "F", ]
Kemper_iv_F_1$Experiment <- "Kemper | 1mg/kg | iv | F"

Kemper_oral_F_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "F", ]
Kemper_oral_F_1$Experiment <- "Kemper | 1mg/kg | oral | F"

Kemper_oral_F_5 <- Kemper_serum_results[Kemper_serum_results$Dose == 5 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "F", ]
Kemper_oral_F_5$Experiment <- "Kemper | 5mg/kg | oral | F"
Kemper_oral_F_25 <- Kemper_serum_results[Kemper_serum_results$Dose == 25 & Kemper_serum_results$Type == "oral" &
                                  Kemper_serum_results$sex == "F", ]
Kemper_oral_F_25$Experiment <- "Kemper | 25mg/kg | oral | F"

Kemper_iv_M_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "iv" &
                               Kemper_serum_results$sex == "M", ]
Kemper_iv_M_1$Experiment <- "Kemper | 1mg/kg | iv | M"

Kemper_oral_M_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "M", ]
Kemper_oral_M_1$Experiment <- "Kemper | 1mg/kg | oral | M"

Kemper_oral_M_5 <- Kemper_serum_results[Kemper_serum_results$Dose == 5 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "M", ]
Kemper_oral_M_5$Experiment <- "Kemper | 5mg/kg | oral | M"
Kemper_oral_M_25 <- Kemper_serum_results[Kemper_serum_results$Dose == 25 & Kemper_serum_results$Type == "oral" &
                                  Kemper_serum_results$sex == "M", ]
Kemper_oral_M_25$Experiment <- "Kemper | 25mg/kg | oral | M"
# 
# Kudo_high <- Kudo_results[Kudo_results$Dose == 16.560, ] 
# Kudo_high$Experiment <-  "Kudo | 16mg/kg | iv"
# 
# Kim_results <- read.csv("Kim_results.csv", header = T)
# Kim_iv <- Kim_results[Kim_results$Type == "iv",]
# Kim_iv$Experiment <-  "Kim | 1mg/kg | iv"
# 
# Kim_oral <- Kim_results[Kim_results$Type == "oral",]
# Kim_oral$Experiment <-  "Kim | 1mg/kg | oral"


Experiment <- scales::hue_pal()(8)

names(Experiment) <- unique(c( Kemper_iv_F_1$Experiment, Kemper_iv_M_1$Experiment,
                              Kemper_oral_F_1$Experiment,Kemper_oral_M_1$Experiment, Kemper_oral_F_5$Experiment, 
                               
                               Kemper_oral_M_5$Experiment,  Kemper_oral_F_25$Experiment, Kemper_oral_M_25$Experiment))



scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5, alpha = 0.7) +  # Identity line in log10 scale
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # +3-fold error line
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # -3-fold error line
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # +10-fold error line
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # -10-fold error line
  geom_point(data = Kemper_iv_F_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_iv_M_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_manual(values = Experiment)+
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))+
  theme_light()+
  labs(y = "Predicted PFOA (mg/L tissue)",
       x = "Observed PFOA ( mg/L tissue)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.75, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=12,  hjust = 0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        )

#print(scatter_plot)

ggsave("validation_plot_PFOA_serum.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


