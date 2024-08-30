# Lifetime equations check
rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PFAS_PARC"
setwd(HOME)

OUTPUT <- file.path("Output/Data", Sys.Date())
dir.create(OUTPUT, recursive = TRUE)

CP_theme <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 7, lineheight = unit(0.2, "lines")), # lineheight is adjusting the space between lines
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.line = element_line(linewidth = 0.05),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      panel.border = element_blank(), 
      panel.background = element_blank(),
      panel.grid = element_line(linewidth = 0.1), 
      strip.background = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0, "cm"),
      #legend.key.width = unit(0.1, "cm"),  # Make legend key width span the whole plot
      legend.key.height = unit(0.2, "cm"),  # Adjust legend key height
      legend.text = element_text(size = 7)
    )
}

# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)


# Duration of lifetime (0 - 80 years old) 
TSTART <- 0
TSTOP <- 365*80 # years in days
DT <- 1
TIME <- seq(TSTART,TSTOP,by=DT)

Variables_df <- as.data.frame(list(TIME = TIME)) #df column 1 = simulation time, every step is 1 day
Variables_df <- Variables_df %>%
  mutate(age = TIME/365) # add column 2 = age in days

# Body weight
Variables_df <- Variables_df %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))
#?? What is the difference between BW and BDW, we actually only use BDW in the script ?? The difference appears after 19.00274 years, where BW starts to be bigger than BDW

# Using Ratier et al. (2024) model, Fraction of arterial plasma, calculated from Filser 2000 p.43
Fr_art_plasma = 0.0178 / (0.0178 + 0.0533)#fraction of arterial blood (corrected for plasma)

# Hematocrit - male                                                              # From Supp mat of Brochot et al. 2019
Param1_M = 33.455469
Param2_M = 53.206039
Param3_M = 8.277945
Param4_M = 40.492556
Param5_M = 46.899695

b1_M = (Param4_M - Param1_M -(Param2_M - Param1_M) * exp(-Param3_M))/5
a1_M = Param4_M - 6*b1_M
b2_M = (Param5_M - Param4_M)/5
a2_M = Param4_M - 15*b2_M

# Hematocrit - non pragnant female
Param1_F = 32.617402
Param2_F = 53.188459
Param3_F = 7.699418
Param4_F = 37.531463
Param5_F = 40.055284

b1_F = (Param4_F - Param1_F - (Param2_F-Param1_F)*exp(-Param3_F))/2
a1_F = Param4_F - 3*b1_F
b2_F = (Param5_F - Param4_F)/7
a2_F = Param5_F - 10*b2_F


### Changes in body composition over time  ####                                            
Variables_df <- Variables_df %>%
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
  ## Hematocrit
  mutate(Hct_ven_M = if_else(age < 1, (Param1_M +(Param2_M-Param1_M)*exp(-Param3_M*age))*0.01,
                             if_else(age < 6, (a1_M + b1_M*age)*0.01,
                                     if_else(age < 15, Param4_M*0.01,
                                             if_else(age < 20, (a2_M + b2_M*age)*0.01,
                                                     Param5_M*0.01))))) %>%
  mutate(Hct_M = Hct_ven_M*0.91) %>%
  mutate(Hct_ven_F = if_else(age < 1, (Param1_F +(Param2_F-Param1_F)*exp(-Param3_F*age))*0.01,
                             if_else(age < 3, (a1_F + b1_F*age)*0.01,
                                     if_else(age < 10, (a2_F + b2_F*age)*0.01,
                                             Param5_F*0.01)))) %>%
  mutate(Hct_F = Hct_ven_F*0.91) %>%
  
  
  ## Cardiac output
  mutate(CardOut_M = if_else(age < 33.37, (6.642 + (0.6 - 6.642)*exp(-0.1323*age))*(1-Hct_M),
                             (-0.000895*age^2 + 0.0607*age + 5.54)*(1-Hct_M))) %>%
  mutate(CardOut_F = if_else(age < 16.027, (7.734 + (0.6 - 7.734)*exp(-0.09747*age))*(1-Hct_F),
                             (0.000473*age^2 - 0.0782*age + 7.37)*(1-Hct_F))) %>%
  
  ## Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
  mutate(VplasmaFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                     0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(Vplasma_M = VplasmaFraction_M*BDW_M) %>%
  mutate(VplasmaFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                     if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                             0.065))) %>%
  mutate(Vplasma_F = VplasmaFraction_F*BDW_F) %>% 
  
  ## Total volume (sum of all actually included organs, different from BDW as we are missing some) 
  mutate(VtotalFraction_M = VplasmaFraction_M) %>%
  mutate(Vtotal_M = VtotalFraction_M*BDW_M) %>% 
  mutate(VtotalFraction_F = VplasmaFraction_F) %>%  
  mutate(Vtotal_F = VtotalFraction_F*BDW_F) %>% 
  
  ## Rest (Rest = Total - Organs included in the model)
  mutate(VrestFraction_M = 1 - VtotalFraction_M) %>%
  mutate(Vrest_M = VrestFraction_M * BDW_M) %>% 
  mutate(VrestFraction_F = 1 - VtotalFraction_F) %>%
  mutate(Vrest_F = VrestFraction_F * BDW_F) %>% 
  
  ## MassBalance Volumes #better be 0
  mutate(Calculated_M = Vplasma_M + Vrest_M) %>% 
  mutate(Calculated_F = Vplasma_F + Vrest_F) %>% 
  mutate(Vmass_balance_M = BDW_M - Calculated_M) %>%
  mutate(Vmass_balance_F = BDW_F - Calculated_F)

  


Variables_df %>% select(age, BDW_M, BDW_F, Vplasma_M, Vplasma_F, Vtotal_M, Vtotal_F, 
                        Vrest_M, Vrest_F, Calculated_M, Calculated_F) %>% 
  ggplot() +
  geom_line(aes(age, Vplasma_M), colour = "orange")+
  geom_line(aes(age, Vplasma_F), colour = "blue") +
  geom_line(aes(age, Vtotal_M), colour = "purple") +
  geom_line(aes(age, Vtotal_F), color = "darkgreen") +
  geom_line(aes(age, Vrest_M), colour = "red") +
  geom_line(aes(age, Vrest_F), colour = "darkred") +
  geom_line(aes(age, BDW_M), colour = "lightblue") +
  geom_line(aes(age, BDW_F), colour = "pink") +
  # geom_line(aes(age, Calculated_M), colour = "darkblue") +
  # geom_line(aes(age, Calculated_F), colour = "black") + 
  CP_theme()

# Plot mass balance organ volumes
VolumeMassBalance <- Variables_df %>% select(age, Vmass_balance_M, Vmass_balance_F) %>% 
  pivot_longer(names_to = "Gender", values_to = "MB", Vmass_balance_M:Vmass_balance_F) %>% 
  ggplot(aes(age, MB)) +
  geom_path() +
  facet_wrap(~Gender) +
  CP_theme()
VolumeMassBalance
ggsave("VolumeMassBalance0.png", dpi = 300)


  # ## Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  # mutate(VliverFraction_M = 0.0247 + (0.0409 - 0.0247)*exp(-0.218*age)) %>%
  # mutate(Vliver_M = VliverFraction_M*BDW_M) %>%
  # mutate(QliverFraction_M = (VliverFraction_M/0.0247)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  # mutate(VliverFraction_F = 0.0233 + (0.038 - 0.0233)*exp(-0.122*age)) %>%
  # mutate(Vliver_F = VliverFraction_F*BDW_F) %>%
  # mutate(QliverFraction_F = (VliverFraction_F/0.0233)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  
  ## Total volume (sum of all actually included organs, different from BDW as we are missing some) 
  mutate(VtotalFraction_M = 1) %>%
  mutate(VtotalFraction_F = 1) %>%  
  
  
  ## Flow rates
  mutate(QtotalFraction_M = 1) %>%
  mutate(QtotalFraction_F = 1) %>%
  
  # mutate(Qliver_M = QliverFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  # mutate(Qliver_F = QliverFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  # #Qliver is what was used to be Qhepatic
  
  

  mutate(QrestFraction_M = 1) %>% 
  mutate(Qrest_M = QrestFraction_M/QtotalFraction_M*CardOut_M) %>% 
  

  mutate(QrestFraction_F = 1) %>%
  mutate(Qrest_F = QrestFraction_F/QtotalFraction_F*CardOut_F) %>% 
  
  
  ## MassBalance Flow #better be 0
  mutate(Qmass_balance_M = CardOut_M - (Qrest_M)) %>% 
  mutate(Qmass_balance_F = CardOut_F - (Qrest_F)) %>% 
  
  ## MassBalance Volumes #better be 0
  mutate(Vmass_balance_M = BDW_M - (Vplasma_M + Vrest_M)) %>%
  mutate(Vmass_balance_F = BDW_F - (Vplasma_F + Vrest_F))



### Plot mass balance organ flow  ####
FlowMassBalance <- Variables_df %>% select(age, Qmass_balance_M, Qmass_balance_F) %>% 
  pivot_longer(names_to = "Gender", values_to = "MB", Qmass_balance_M:Qmass_balance_F) %>% 
  ggplot(aes(age, MB, colour = Gender)) +
  geom_path() +
  CP_theme()
FlowMassBalance
ggsave("FlowMassBalance0.png", dpi = 300)


# 
# Vtotal_M = 
#   Vblood_M + 
#   Vliver_M + 
#   Vstomach_M + 
#   Vgut_M + 
#   Vkidney_M + 
#   Vurinarytract_M +
#   Vskin_M + 
#   Vadipose_M + 
#   Vadrenal_M + 
#   Vbone_M + 
#   Vbonenonperfused_M + 
#   Vbrain_M + 
#   Vbreast_M + 
#   Vheart_M + 
#   Vmarrow_M + 
#   Vmuscle_M + 
#   Vrepro_M + 
#   Vpancreas_M + 
#   Vspleen_M + 
#   Vthyroid_M + 
#   Vlung_M


