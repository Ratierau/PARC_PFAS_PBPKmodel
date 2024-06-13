library(ggplot2)

#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Results of Scenarios Validation/Scenario_3")

# Load Results 
Cui_results <- read.csv("Cui_2008_results.csv", header = T)
Cui_low <- Cui_results[Cui_results$Dose == 5, ]
Cui_low$Experiment <- "Cui | 5mg/kg | oral | M"
Cui_high <- Cui_results[Cui_results$Dose == 20, ]
Cui_high$Experiment <- "Cui | 20mg/kg | oral | M"
# 
Tissue_markers <-  c(0:10,14)
names(Tissue_markers) <- c( "Lung",   "Spleen",  "Liver",   "Kidney", "Plasma", "Heart",
                            "Brain","Testis", "Stomach", "Intestines", "Carcass")

Experiment <- scales::hue_pal()(2)

names(Experiment) <- unique(c(Cui_low$Experiment, Cui_high$Experiment))



scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5, alpha = 0.7) +  # Identity line in log10 scale
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # +3-fold error line
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # -3-fold error line
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # +10-fold error line
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # -10-fold error line
  geom_point(data = Cui_low, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Cui_high, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_manual(values = Experiment)+
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))+
  scale_shape_manual(values = Tissue_markers)+
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

ggsave("validation_plot_PFOA_tissues.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


