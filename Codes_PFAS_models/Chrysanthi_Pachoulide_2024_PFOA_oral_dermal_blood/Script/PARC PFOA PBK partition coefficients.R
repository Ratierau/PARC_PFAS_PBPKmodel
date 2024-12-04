# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Calculating Tissue:Plasma Partition Coefficients
# CP, 29-11-2024
# --------------------------------------------------------------------------- #


rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
INPUT = file.path("Input", Sys.Date(), workingtime)
dir.create(INPUT, recursive = TRUE)
setwd(INPUT)


# Load packages
library(ggplot2)
library(tidyverse)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()
theme_CP <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 18, lineheight = unit(0.5, "lines")), # lineheight is adjusting the space between lines
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.line = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      panel.border = element_blank(), 
      panel.background = element_rect((fill = "grey94")),
      panel.grid = element_line(linewidth = 0.1,5, colour = "grey100"), 
      strip.background = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0, "cm"),
      legend.key.width = unit(0.2, "cm"),  # Make legend key width span the whole plot
      legend.key.height = unit(0.2, "cm"),  # Adjust legend key height
    )
}


# Tissue Composition Physiology

Physio.data <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/TissueComposition.csv")
# Physio.Allend <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/TissueCompositionAsAllendorf.csv")
# Physio.Utsey <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/TissueCompositionAsUtsey.csv")


# PFOA Matrix/Water Distribution Coefficients

logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
logSP <- 1.61 # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954

k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
k_SP <- 10^logSP       # structural protein:water partition coefficient
k_ALB <- 10^logALB     # albumin:water partition coefficient
k_SL <- 10^logSL       # structural lipids:water partition coefficient
k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient

# D_MatrixWater.df <- data.frame(
#   Parameter = c("ML", "SL", "SP", "ALB", "FABP"),
#   Mean = c(3.52, -1.37, 1.61, 4.33, 4.3),
#   STD = c(0.08, 0.01, 0.15, 0.05, 0)
# )
# 
# D_MatrixWater.plot <- ggplot(D_MatrixWater.df, aes(x = Parameter, y = Mean)) +
#   geom_boxplot(aes(lower = Mean - STD, 
#                    middle = Mean, 
#                    upper = Mean + STD, 
#                    ymin = Mean - 2*STD, 
#                    ymax = Mean + 2*STD),
#                stat = "identity", 
#                fill = "white", color = "grey2") +
#   labs(x = "", y = "logDistribution: Matrix/Water") +
#   scale_y_continuous(breaks = seq(-1, max(D_MatrixWater.df$Mean + D_MatrixWater.df$STD), by = 1)) + 
#   scale_x_discrete(labels = c("ML" = "Membrane Lipids",
#                               "SL" = "Storage Lipids",
#                               "SP" = "Structural Proteins", 
#                               "ALB" = "Albumin",
#                               "FABP" = "Fatty Acid Binding Proteins")) +
#   coord_flip()+
#   theme_CP()
# D_MatrixWater.plot
# ggsave("MatrixDistribution.png")

fup <- 0.006 # Fisher 2024

## Tissue:Plasma partitioning ####

### As per Allendorf ####

Physio.dat <- Physio.data %>% filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>% 
  select(-Comment) %>% 
  mutate(across(-Tissue, as.numeric))

Physio.Tissues <- Physio.dat %>% filter(!Tissue %in% c("Blood", "Plasma")) 

Kp <- (Physio.Tissues$f_W +
         (k_ML*Physio.Tissues$f_ML) +
         (k_SL*Physio.Tissues$f_SL) +
         (k_SP*Physio.Tissues$f_SP) +
         (k_ALB*Physio.Tissues$f_ALB) +
         (k_FABP*Physio.Tissues$f_FABP)
) * fup

Names <- Physio.Tissues$Tissue %>% substr(1,2)
Names <- paste("Kp",Names ,sep="")

Kp.df <- as.data.frame(list(Tissue = Names)) %>%
  mutate(Value = Kp)
write.csv(Kp.df, "Kp.csv", row.names = FALSE)


# PlotKps <- Kp.df %>% 
#   ggplot(aes(x = reorder(Tissue, Value), y = Value)) +
#   geom_bar(stat = "identity", fill = "grey75") +
#   scale_y_continuous(breaks = seq(0, 1.5, by = 0.25)) + 
#   scale_x_discrete(labels = c("KpSk" = "Skin",
#                               "KpBr" = "Brain",
#                               "KpAd" = "Adipose",
#                               "KpKi" = "Kidney",
#                               "KpLi" = "Liver",
#                               "KpLu" = "Lung",
#                               "KpSp" = "Spleen",
#                               "KpGu" = "Gut",
#                               "KpGo" = "Gonads",
#                               "KpMu" = "Muscle",
#                               "KpHe" = "Heart",
#                               "KpBo" = "Bone")) +
#   coord_flip() +
#   labs(
#     x = "",
#     y = "Kp"
#   ) +
#   theme_CP()
# ggsave("TissueDistribution.png")


Kp.df <- Kp.df %>% pivot_wider(names_from = Tissue, values_from = Value)

# To calculate Kp for the rest compartment
# The current PFOA model only has: skin, adipose, kidney, liver, gut and plasma compartments, all other compartments are lumped together
MaleVariables_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-12/MaleVariables.csv") %>% 
  as.data.frame() %>% 
  filter(age == 40) 

KpRe = (Kp.df$KpBr * MaleVariables_df$V_brain_M +
          Kp.df$KpHe*MaleVariables_df$V_heart_M +
              Kp.df$KpLu*MaleVariables_df$V_lung_M +
              Kp.df$KpMu*MaleVariables_df$V_muscle_M +
              Kp.df$KpSp*MaleVariables_df$V_spleen_M +
              Kp.df$KpGo*MaleVariables_df$V_repro_M +
              Kp.df$KpBo*MaleVariables_df$V_bone_M) / (
                MaleVariables_df$V_brain_M +
                  MaleVariables_df$V_heart_M +
                  MaleVariables_df$V_lung_M +
                  MaleVariables_df$V_muscle_M +
                  MaleVariables_df$V_spleen_M +
                  MaleVariables_df$V_repro_M +
                  MaleVariables_df$V_bone_M)  # Plasma/rest of the body partition coefficient

Kp.df$KpRe <- KpRe

PFOAPBK.PartCoefs.df <- Kp.df %>% 
  select(KpAd, KpGu, KpKi, KpLi, KpSk, KpRe)
write.csv(PFOAPBK.PartCoefs.df, "PFOAPBK.PartCoefs.csv", row.names = FALSE)

# ### As per Utsey - Schmitt ####
# Physio.dat <- Physio.Utsey %>%
#   mutate(across(-tissue, as.numeric))
# 
# Physio.Tissues <- Physio.dat %>% filter(!tissue %in% c("RBCs", "Plasma")) %>%
#   rename(Tissue = tissue) %>%
#   rename(f_W = f_water) %>%
#   mutate(f_SL = f_lipids - f_pl) %>%
#   mutate(f_ML = f_n_pl) %>%
#   mutate(f_Prot = f_proteins) %>%
#   select(Tissue, f_W, f_SL, f_ML, f_Prot)
# 
# Kp <- (Physio.Tissues$f_W +
#          (k_ML*Physio.Tissues$f_ML) +
#          # (k_SP*Physio.Tissues$f_SP) +
#          (k_prot*Physio.Tissues$f_Prot) +
#          # (k_ALB*Physio.Tissues$f_ALB) +
#          (k_SL*Physio.Tissues$f_SL)
#          # (k_FABP*Physio.Tissues$f_FABP)
#        ) * fup
# 
# Names <- Physio.Tissues$Tissue %>% substr(1,2)
# Names <- paste("Kp",Names ,sep="")
# 
# KpAsUtsey.df <- as.data.frame(list(Tissue = Names)) %>%
#   mutate(Value = Kp)
# write.csv(KpAsUtsey.df, "KpAsUtsey.csv", row.names = FALSE)
# 
# 
# 
# 
# 
# ### As per Allendorf simplified ####
# Physio.dat <- Physio.Allend %>% filter(Tissue != "Comment") %>% # as per Allendorf
#   mutate(across(-Tissue, as.numeric))
# 
# Physio.Tissues <- Physio.dat %>% filter(!Tissue %in% c("Blood", "Plasma")) %>%
#   mutate(f_Prot = f_ALB + f_SP + f_FABP) %>%
#   select(Tissue, f_W, f_SL, f_ML, f_Prot)
# 
# fup <- 0.006 # Fisher 2024
# 
# 
# logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
# logSP <- 1.61 # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
# logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
# logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
# logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
# 
# k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
# k_SP <- 10^logSP       # structural protein:water partition coefficient
# k_ALB <- 10^logALB     # albumin:water partition coefficient
# k_SL <- 10^logSL       # structural lipids:water partition coefficient
# k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
# k_prot <- (k_ALB + k_SP + k_FABP)/3
# 
# 
# Kp <- (Physio.Tissues$f_W +
#          (k_ML*Physio.Tissues$f_ML) +
#          (k_prot*Physio.Tissues$f_Prot) +
#          (k_SL*Physio.Tissues$f_SL)
#        ) * fup
# 
# Names <- Physio.Tissues$Tissue %>% substr(1,2)
# Names <- paste("Kp",Names ,sep="")
# 
# KpAsAllendorfSimpl.df <- as.data.frame(list(Tissue = Names)) %>%
#   mutate(Value = Kp)
# write.csv(KpAsAllendorfSimpl.df, "KpAsAllendorfSimpl.csv", row.names = FALSE)
# 
# 
# 
# 
