# Comparison of the free/unbound fractions and partition coefficients
# Original free fractions were from Trine; as previously calculated in Loccisano
# New fraction unbound is taken from the latest Fischer et al. (fu plasma https://doi.org/10.1021/acs.est.3c07415 ) and Rye et al. (fu tissues https://doi.org/10.1021/acs.est.4c04050) papers. 


rm(list=ls()) # to clear out the global environment


# set work directory
HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Script"
setwd(HOME)

# creating a new folder to store the results 
WhatAmITesting <- "COMPARE_Free.fu_PCs"
Saveoutput <- file.path('C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Results', Sys.Date(), WhatAmITesting)
dir.create(Saveoutput, WhatAmITesting, recursive = TRUE)

# bring my plotting function
source("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Script/Theme_CP_function.R")

# load packages
library(lubridate)
library(ggplot2)
library(writexl)
library(ggpubr)
library(tidyverse)
library(broom)
library(patchwork)

# OLD ####

Free = 0.02 # Free fraction of PFOA in plasma

PL_old =2.2  # Plasma/liver partition coefficient
PF_old = 0.04  # Plasma/fat partition coefficient
PK_old = 1.05  # Plasma/kidney partition coefficient
PSk_old = 0.1  # Plasma/skin partition coefficient
PR_old = 0.12  # Plasma/rest of the body partition coefficient
PG_old = 0.05  # Plasma/gut partition coefficient

## Free fraction of chemical in tissues ##

FreeL = Free/PL_old  # liver
# FreeF = Free/PF  # fat
FreeK = Free/PK_old  # kidney
# FreeSk = Free/PSk  # skin
# FreeR = Free/PR  # rest of the body
# FreeG = Free/PG  # gut


# NEW ####

fu = 0.061 # Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415

fuL = 0.0360 # Rye et al. 2024 https://doi.org/10.1021/acs.est.4c04050)
fuK = 0.0565 # Rye et al. 2024 https://doi.org/10.1021/acs.est.4c04050)

PL_new = fu/fuL
PK_new = fu/fuK

Partition_coefficient = c(PL_old, PK_old, PL_new, PK_new)
Fraction_unbound = c(FreeL, FreeK, fuL, fuK)
Tissue = c("Liver", "Kidney", "Liver", "Kidney")
Version = c("Old", "Old", "New", "New")

df <- cbind(Partition_coefficient, Fraction_unbound, Tissue, Version) %>%  
  as.data.frame() %>% 
  mutate(across(c(Partition_coefficient, Fraction_unbound), as.numeric))

comparison_results <- df %>%
  pivot_longer(cols = c(Partition_coefficient, Fraction_unbound),
               names_to = "What", values_to = "Value") %>%
  pivot_wider(names_from = Version, values_from = Value) %>%
  mutate(Difference = (Old - New)/Old*100)


plot_fu <- ggplot(df, aes(x = Tissue, y = Fraction_unbound, fill = Version)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  labs(title = "fu", y = "Fraction Unbound", x = "Tissue") +
  theme_minimal() +
  scale_fill_manual(values = c("New" = "black", "Old" = "grey")) +
  theme_CP() 
plot_fu

# Plot for Partition Coefficient
plot_pc <- ggplot(df, aes(x = Tissue, y = Partition_coefficient, fill = Version)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  labs(title = "Partition Coefficient", y = "Partition Coefficient", x = "Tissue") +
  theme_minimal() +
  scale_fill_manual(values = c("New" = "black", "Old" = "grey"))+
  theme_CP() 
plot_pc

# Combine the two plots vertically
combined_plot <- plot_fu / plot_pc  # Patchwork syntax for vertical stacking
combined_plot

ggsave(filename=file.path(Saveoutput, "NewOld_Fu_PC.png"),
       dpi = 300)
