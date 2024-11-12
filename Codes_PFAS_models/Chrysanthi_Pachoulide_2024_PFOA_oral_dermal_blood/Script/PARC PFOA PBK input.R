# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Input file
# CP, 10-11-2024
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set input storage directory
INPUT <- file.path("Input/Data", Sys.Date())
dir.create(INPUT, recursive = TRUE)
setwd(INPUT)


# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)


# PHYSIOLOGICAL PARAMETERS ####
# ------------------------------------------------------ #

## Lifetime equations ####

lifeTSTOP <- 80 # duration of lifetime (0 - 80 years old)  of simulation
TSTART <- 0
TSTOP <- 365*lifeTSTOP # years in days
DT <- 1
TIME <- seq(TSTART,TSTOP,by=DT)


# Creating a dataframe for all the variables
Variables_df <- as.data.frame(list(TIME = TIME)) #df column 1 = simulation time, every step is 1 day
Variables_df <- Variables_df %>%
  mutate(age = TIME/365) # add column 2 = age in days

Variables_df <- Variables_df %>%
  
  ### Fractional Volumes ####
  # Body weight
  # BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
  mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                    -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))
  
# Blood/Plasma/Hematocrit 
# Using Ratier et al. (2024) model, Fraction of arterial plasma, calculated from Filser 2000 p.43
Fr_art_plasma = 0.0178 / (0.0178 + 0.0533) #fraction of arterial blood (corrected for plasma)
  
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

Variables_df <- Variables_df %>%
  select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
  rename(BW_M = BW_M_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BW_F = BW_F_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
  # Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_adrenalFraction_M = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(V_adrenalFraction_F = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  
 # Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_boneFraction_M = (0.313 + (0.506 - 0.313)*exp(-0.0907*age))*0.095) %>%
  mutate(V_bonenonperfusedFraction_M = 0.095 - V_boneFraction_M) %>%
  mutate(V_boneFraction_F = (0.298 + (0.505 - 0.298)*exp(-0.0792*age))*0.085) %>%
  mutate(V_bonenonperfusedFraction_F = 0.085 - V_boneFraction_F) %>%
  
  # Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_brainFraction_M = (1.450 + (0.353 - 1.450) * exp (-0.440*age))/BDW_M) %>%
  mutate(V_brainFraction_F = (1.300 + (0.347 - 1.300) * exp (-0.573*age))/BDW_F) %>%
  
  # Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_breastFraction_M = 3.42E-4*1/(1 + exp(-1.42*age + 20.1))) %>%
  mutate(V_breastFraction_F = 0.00833/(1 + exp(-1.92*age+ 28.6))) %>%
  
  # Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_heartFraction_M = 0.0045) %>%
  mutate(V_heartFraction_F = 0.0042) %>%
  
  # Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_marrowFraction_M = 0.05 + (0.0138 - 0.05)*exp(-0.112*age)) %>%
  mutate(V_marrowFraction_F = 0.045 + (0.0138 - 0.045)*exp(-0.136*age)) %>%
  
  # Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(MuscleAtrophy_M = if_else(age < 24.3, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(V_muscleFraction_M = (0.3973 + (0.201 - 0.3973)*exp(-0.141*age)) * MuscleAtrophy_M) %>%
  mutate(MuscleAtrophy_F = if_else(age < 25.90709, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(V_muscleFraction_F = (0.2917 + (0.207 - 0.2917)*exp(-0.339*age)) * MuscleAtrophy_F) %>%
  
  # Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_reproFraction_M = if_else(age < 20.01, -1.5156E-07*age^3 + 9.3351E-06*age^2 - 1.1177E-04*age + 4.7966E-04,
                                    0.0008)) %>%
  mutate(V_reproFraction_F = if_else(age < 1, -1.064E-3*age + 1.338E-3,
                                    if_else(age < 20, 2.6380E-7*age^3 - 1.7943E-6*age^2 - 5.6465E-6*age + 2.8105E-4,
                                            0.001552))) %>%
  
  # Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_pancreasFraction_M = 0.00192) %>%
  mutate(V_pancreasFraction_F = 0.002) %>%
  
  # Skin; compartment [10] in Ratier 2024
  mutate(V_skinFraction_M = if_else(age < 20.01, -1.1706E-05*age^3 + 5.4130E-04*age^2 - 6.1966E-03*age + 4.6231E-02,
                                   0.0452)) %>%
  mutate(V_skinFraction_F = if_else(age < 19.45, -7.8882E-06*age^3 + 4.0224E-04*age^2 - 5.2146E-03*age + 4.5605E-02,
                                   0.0383)) %>%

  # Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_spleenFraction_M = 0.0021) %>%
  mutate(V_spleenFraction_F = 0.0022) %>%
  
  # Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_thyroidFraction_M = 0.000274) %>%
  mutate(V_thyroidFraction_F = 0.0003) %>%
  
  # Urinary tract (bladder, ureters, urethra); compartment [13] in Ratier 2024
  mutate(V_urinarytractFraction_M = 0.00104) %>%
  mutate(V_urinarytractFraction_F = 0.0010) %>%
  
  # Kidney; compartment [14] in Ratier 2024
  mutate(V_kidneyFraction_M = 0.0042 + (0.00767 - 0.0042)*exp(-0.206*age)) %>%
  mutate(V_kidneyFraction_F = 0.0046 + (0.0071 - 0.0046)*exp(-0.221*age)) %>%
  
  # Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  # Comment Chrysa on 18-10-2024: Shouldn't the lungs take 100% of the blood flow?
  mutate(V_lungFraction_M = 0.0068) %>%
  mutate(V_lungFraction_F = 0.0070) %>%
  
  # Gut; compartment [16] in Ratier 2024
  mutate(V_gutFraction_M = if_else(age < 16, -0.000082562*age^2 + 0.0013523*age + 0.01293,
                                  0.0140)) %>%
  mutate(V_gutFraction_F = if_else(age < 14.453301, -7.42E-5*age^2 + 1.28E-3*age + 1.30E-2,
                                  0.0160)) %>%
  
  # Stomach; compartment [17] in Ratier 2024
  mutate(V_stomachFraction_M = 0.0021) %>%
  mutate(V_stomachFraction_F = 0.0023) %>%
  
  # Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(V_liverFraction_M = 0.0247 + (0.0409 - 0.0247)*exp(-0.218*age)) %>%
  mutate(V_liverFraction_F = 0.0233 + (0.038 - 0.0233)*exp(-0.122*age)) %>%
  
  # Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
  mutate(V_plasmaFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                     0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(V_plasmaFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                     if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                             0.065))) %>%
  
  # Adipose tissue 
  mutate(V_adiposeFraction_M = 0.96 - V_adrenalFraction_M - V_boneFraction_M - V_bonenonperfusedFraction_M - V_brainFraction_M - V_breastFraction_M +
           - V_heartFraction_M - V_marrowFraction_M - V_muscleFraction_M - V_reproFraction_M - V_pancreasFraction_M +
           - V_skinFraction_M - V_spleenFraction_M - V_thyroidFraction_M - V_urinarytractFraction_M - V_kidneyFraction_M +
           - V_lungFraction_M - V_gutFraction_M - V_stomachFraction_M - V_liverFraction_M - V_plasmaFraction_M) %>%
  mutate(AdiposeMass_M = if_else(age < 19.00093277, 0,
                                 (-0.01129273*age^2 + 1.11817056*age + 56.74397436) - (74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))))) %>% # age and not BW dependent 
  mutate(V_adiposeFraction_F = 0.96 - V_adrenalFraction_F - V_boneFraction_F - V_bonenonperfusedFraction_F - V_brainFraction_F - V_breastFraction_F +
           - V_heartFraction_F - V_marrowFraction_F - V_muscleFraction_F - V_reproFraction_F - V_pancreasFraction_F +
           - V_skinFraction_F - V_spleenFraction_F - V_thyroidFraction_F - V_urinarytractFraction_F - V_kidneyFraction_F +
           - V_lungFraction_F - V_gutFraction_F - V_stomachFraction_F - V_liverFraction_F - V_plasmaFraction_F) %>%
  mutate(AdiposeMass_F = if_else(age < 17.9374115, 0,
                                 if_else(((-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))) < 0, 0,
                                         (-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))))) # what does this mean?!!
  
  
### Fractional Blood Flows ####

  Variables_df <- Variables_df %>% 
  
  # Hematocrit
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
  
  # Cardiac output (plasma; L/min*60*24 = L/d)
  mutate(CardOut_M = if_else(age < 33.37, (6.642 + (0.6 - 6.642)*exp(-0.1323*age))*(1-Hct_M)*60*24,
                             (-0.000895*age^2 + 0.0607*age + 5.54)*(1-Hct_M)*60*24)) %>%
  mutate(CardOut_F = if_else(age < 16.027, (7.734 + (0.6 - 7.734)*exp(-0.09747*age))*(1-Hct_F)*60*24,
                             (0.000473*age^2 - 0.0782*age + 7.37)*(1-Hct_F)*60*24)) %>%
  
  # Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_adrenalFraction_M = (V_adrenalFraction_M/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_adrenalFraction_F = (V_adrenalFraction_F/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_boneFraction_M = (V_boneFraction_M/(0.095*0.32))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_boneFraction_F = (V_boneFraction_F/(0.085*0.298))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_brainFraction_M = (V_brainFraction_M/0.01986)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_brainFraction_F = (V_brainFraction_F/0.0217)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_breastFraction_M = (V_breastFraction_M/0.00035)*0.0002) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_breastFraction_F = (V_breastFraction_F/0.0083)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_heartFraction_M = (V_heartFraction_M/0.0045)*0.041) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_heartFraction_F = (V_heartFraction_F/0.0042)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_marrowFraction_M = (V_marrowFraction_M/0.050)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_marrowFraction_F = (V_marrowFraction_F/0.045)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_muscleFraction_M = (V_muscleFraction_M/0.3973)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_muscleFraction_F = (V_muscleFraction_F/0.2917)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_reproFraction_M = (V_reproFraction_M/0.0008)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_reproFraction_F = (V_reproFraction_F/0.0016)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_pancreasFraction_M = (V_pancreasFraction_M/0.00192)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_pancreasFraction_F = (V_pancreasFraction_F/0.002)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Skin; compartment [10] in Ratier 2024
  mutate(Q_skinFraction_M = (V_skinFraction_M/0.0452)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_skinFraction_F = (V_skinFraction_F/0.0383)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_spleenFraction_M = (V_spleenFraction_M/0.0021)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_spleenFraction_F = (V_spleenFraction_F/0.0022)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_thyroidFraction_M = (V_thyroidFraction_M/0.000274)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_thyroidFraction_F = (V_thyroidFraction_F/0.0003)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Urinary tract (bladder, ureters, urethra); compartment [13] in Ratier 2024
  mutate(Q_urinarytractFraction_M = (V_urinarytractFraction_M/0.00104)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_urinarytractFraction_F = (V_urinarytractFraction_F/0.0010)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Kidney; compartment [14] in Ratier 2024
  mutate(Q_kidneyFraction_M = (V_kidneyFraction_M/0.0042)*0.196) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_kidneyFraction_F = (V_kidneyFraction_F/0.0046)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  # Comment Chrysa on 18-10-2024: Shouldn't the lungs take 100% of the blood flow?
  mutate(Q_lungFraction_M = (V_lungFraction_M/0.0068)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_lungFraction_F = (V_lungFraction_F/0.0070)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Gut; compartment [16] in Ratier 2024
  mutate(Q_gutFraction_M = (V_gutFraction_M/0.0140)*0.144) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_gutFraction_F = (V_gutFraction_F/0.0160)*0.165) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Stomach; compartment [17] in Ratier 2024
  mutate(Q_stomachFraction_M = (V_stomachFraction_M/0.0021)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_stomachFraction_F = (V_stomachFraction_F/0.0023)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(Q_liverFraction_M = (V_liverFraction_M/0.0247)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  mutate(Q_liverFraction_F = (V_liverFraction_F/0.0233)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  
  # Adipose tissue 
  mutate(Q_adiposeFraction_M = (V_adiposeFraction_M/0.20)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_adiposeFraction_F = (V_adiposeFraction_F/0.3167)*0.087) # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  

### Check mass balance volumes and flows ####
Variables_M_df <- Variables_df %>%  
  select(.,ends_with("_M")) %>% 
  mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% 
  mutate(VolumesSum = rowSums(select(., starts_with("V_")))) %>% 
  mutate(TIME = TIME) %>% 
  mutate(age = TIME/365) 

Variables_F_df <- Variables_df %>% 
  select(.,ends_with("_F")) %>% 
  mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% 
  mutate(VolumesSum = rowSums(select(., starts_with("V_")))) %>% 
  mutate(TIME = TIME) %>% 
  mutate(age = TIME/365)

PLOT_VolumeTotal <-  
  ggplot() + 
  geom_path(data = Variables_M_df, aes(TIME, VolumesSum), colour = "lavenderblush4") +
  geom_path(data = Variables_F_df, aes(TIME, VolumesSum), colour = "purple")
PLOT_VolumeTotal # Is 1

PLOT_BloodFlowTotal <-  
  ggplot()+ 
  geom_path(data = Variables_M_df, aes(TIME, BloodFlowSum), colour = "lavenderblush4") +
  geom_path(data = Variables_F_df, aes(TIME, BloodFlowSum), colour = "purple")
PLOT_BloodFlowTotal # Not 1; shouldn't it be 1?


### Actual Volumes and Blood Flows ####
MaleVariables_df <- Variables_df %>% 
  select(ends_with('_M')) %>% 
  select(!c(V_boneFraction_M, V_bonenonperfusedFraction_M, V_adiposeFraction_M, V_plasmaFraction_M)) %>% 
  mutate(across(starts_with('V_'), ~ . * BDW_M)) %>% # Final Organ Volume = Fractional Organ Volume * Body Weight (age specific)
  mutate(V_boneFraction_M = Variables_df$V_boneFraction_M*BDW_M/2, # 2 is the bone density
         V_bonenonperfusedFraction_M = Variables_df$V_bonenonperfusedFraction_M*BDW_M/2, # 2 is the bone density
         V_adiposeFraction_M = (Variables_df$V_adiposeFraction_M*BDW_M/0.9) + (AdiposeMass_M/0.9), # 0.9 is the bone density
         V_plasmaFraction_M = Variables_df$V_plasmaFraction_M*BDW_M*(1-Variables_df$Hct_M)) %>%  # (1-Hematocrite) corrects for the plasma volume (if not is total blood) 
  mutate(TotalBloodFlow = Variables_M_df$BloodFlowSum) %>% 
  mutate(across(starts_with("Q_"), ~ . /TotalBloodFlow * CardOut_M))  # Final Blood Flow = Fractional Blood Flow/Total Fractional Blood Flow * Cardiac Output (age specific)
  #select((matches("(Q_|V_)"))) 
colnames(MaleVariables_df) <- str_remove(colnames(MaleVariables_df), "Fraction")

FemaleVariables_df <- Variables_df %>% 
  select(ends_with('_F')) %>% 
  select(!c(V_boneFraction_F, V_bonenonperfusedFraction_F, V_adiposeFraction_F, V_plasmaFraction_F)) %>% 
  mutate(across(starts_with('V_'), ~ . * BDW_F)) %>% # Final Organ Volume = Fractional Organ Volume * Body Weight (age specific)
  mutate(V_boneFraction_F = Variables_df$V_boneFraction_F*BDW_F/2, # 2 is the bone density
         V_bonenonperfusedFraction_F = Variables_df$V_bonenonperfusedFraction_F*BDW_F/2, # 2 is the bone density
         V_adiposeFraction_F = (Variables_df$V_adiposeFraction_F*BDW_F/0.9) + (AdiposeMass_F/0.9), # 0.9 is the bone density
         V_plasmaFraction_F = Variables_df$V_plasmaFraction_F*BDW_F*(1-Variables_df$Hct_F)) %>%  # (1-Hematocrite) corrects for the plasma volume (if not is total blood) 
  mutate(TotalBloodFlow = Variables_F_df$BloodFlowSum) %>% 
  mutate(across(starts_with("Q_"), ~ . /TotalBloodFlow * CardOut_F)) %>%  # Final Blood Flow = Fractional Blood Flow/Total Fractional Blood Flow * Cardiac Output (age specific)
  select((matches("(Q_|V_)"))) 
colnames(FemaleVariables_df) <- str_remove(colnames(FemaleVariables_df), "Fraction")


### PFAS SPECIFIC ####

# Model compartments: Skin, Kidney (K.Plasma, K.Tissue, K.Filtrate), Gut, Liver, Plasma, Adipose

# Rest compartment

 

mutate(Vrest_M = Vtotal_M - Vskin_M - Vkidney_M - Vgut_M - Vliver_M - Vplasma_M - Vadipose_M) %>%
  mutate(Qrest_M = Qtotal_M - Qskin_M - Qkidney_M - Qhepatic_M - Qgut_M - Qadipose_M) %>%
  mutate(Vrest_F = Vtotal_F - Vskin_F - Vkidney_F - Vgut_F - Vliver_F - Vplasma_F - Vadipose_F) %>%
  mutate(Qrest_F = Qtotal_F - Qskin_F - Qkidney_F - Qhepatic_F - Qgut_F - Qadipose_F) %>%
  


  


## Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
mutate(VplasmaFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                   0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(Vplasma_M = VplasmaFraction_M*BDW_M) %>%
  mutate(Vart_M = Vplasma_M*Fr_art_plasma*(1-Hct_M)) %>% #arterial blood (plasma)
  mutate(Vven_M = Vplasma_M*(1-Hct_M) - Vart_M) %>% #venous blood (plasma)
  mutate(VplasmaFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                     if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                             0.065))) %>%
  mutate(Vplasma_F = VplasmaFraction_F*BDW_F) %>%
  mutate(Vart_F = Vplasma_F*Fr_art_plasma*(1-Hct_F)) %>%
  mutate(Vven_F = Vplasma_F*(1-Hct_F) - Vart_F) %>%
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  mutate(Vtotal_M = Vplasma_M + Vadrenal_M + Vbone_M + Vbrain_M + Vbreast_M + 
           Vheart_M + Vmarrow_M + Vmuscle_M + Vrepro_M + Vpancreas_M +
           Vskin_M + Vspleen_M + Vthyroid_M + Vurinarytract_M + Vkidney_M +
           Vlung_M + Vgut_M + Vstomach_M + Vliver_M + Vadipose_M) %>%
  mutate(Vtotal_F = Vplasma_F + Vadrenal_F + Vbone_F + Vbrain_F + Vbreast_F + 
           Vheart_F + Vmarrow_F + Vmuscle_F + Vrepro_F + Vpancreas_F +
           Vskin_F + Vspleen_F + Vthyroid_F + Vurinarytract_F + Vkidney_F +
           Vlung_F + Vgut_F + Vstomach_F + Vliver_F + Vadipose_F) %>%
  
  ## Flow rates
  mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M +
           QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
           QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
           QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M) %>%
  mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F +
           QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
           QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
           QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F) %>%
  
  mutate(TestQF_M = CardOut_M - QtotalFraction_M) %>% 
  mutate(TestQF_F = CardOut_F - QtotalFraction_F) %>% 
  
  
  mutate(Qliver_M = (QliverFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  # Comment Chrysa: Qliver is what was used to be Qhepatic
  mutate(Qstomach_M = (QstomachFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_M = (QgutFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_M = (QkidneyFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_M = (QurinarytractFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_M = (QskinFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_M = (QadrenalFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_M = (QboneFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_M = (QbrainFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_M = (QbreastFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_M = (QheartFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_M = (QmarrowFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_M = (QmuscleFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_M = (QreproFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_M = (QpancreasFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_M = (QspleenFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_M = (QthyroidFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_M = (QlungFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_M = (QadiposeFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  
  mutate(Qliver_F = (QliverFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qstomach_F = (QstomachFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_F = (QgutFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_F = (QkidneyFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_F = (QurinarytractFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_F = (QskinFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_F = (QadrenalFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_F = (QboneFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_F = (QbrainFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_F = (QbreastFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_F = (QheartFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_F = (QmarrowFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_F = (QmuscleFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_F = (QreproFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_F = (QpancreasFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_F = (QspleenFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_F = (QthyroidFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_F = (QlungFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_F = (QadiposeFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  
  ## Total flow
  mutate(Qtotal_M = Qadrenal_M + Qbone_M + Qbrain_M + Qbreast_M + 
           Qheart_M + Qmarrow_M + Qmuscle_M + Qrepro_M + Qpancreas_M +
           Qskin_M + Qspleen_M + Qthyroid_M + Qurinarytract_M + Qkidney_M +
           Qlung_M + Qgut_M + Qstomach_M + Qliver_M + Qadipose_M) %>%
  mutate(Qtotal_F = Qadrenal_F + Qbone_F + Qbrain_F + Qbreast_F + 
           Qheart_F + Qmarrow_F + Qmuscle_F + Qrepro_F + Qpancreas_F +
           Qskin_F + Qspleen_F + Qthyroid_F + Qurinarytract_F + Qkidney_F +
           Qlung_F + Qgut_F + Qstomach_F + Qliver_F + Qadipose_F) %>% 
  
  
# Physicochemical
  mutate(CLdermalabs_M = ((Papp*fSkbarea_M)/1000)*24) %>% # (L/d) ; cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  mutate(CLdermalabs_F = ((Papp*fSkbarea_F)/1000)*24) %>% # (L/d) ; cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  
  ## Rest (Rest = Total - Organs included in the model)
  # mutate(VrestFraction_M = VtotalFraction_M - VskinFraction_M - VurinarytractFraction_M - VkidneyFraction_M - VgutFraction_M - VliverFraction_M - VplasmaFraction_M - VadiposeFraction_M) %>%
  # mutate(Vrest_M = VrestFraction_M*BDW_M) %>%
  # mutate(QrestFraction_M = Qtotal_M - Qskin_M - Qurinarytract_M - Qkidney_M - Qhepatic_M - Qgut_M - Qadipose_M) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  # mutate(Qrest_M = QrestFraction_M/QtotalFraction_M*CardOut_M) %>% 
  # mutate(VrestFraction_F = VtotalFraction_F - VskinFraction_F - VurinarytractFraction_F - VkidneyFraction_F - VgutFraction_F - VliverFraction_F - VplasmaFraction_F - VadiposeFraction_F) %>%
  # mutate(Vrest_F = VrestFraction_F*BDW_F) %>%
  # mutate(QrestFraction_F = Qtotal_F - Qskin_F - Qurinarytract_F - Qkidney_F - Qhepatic_F - Qgut_F - Qadipose_F) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  # mutate(Qrest_F = QrestFraction_F/QtotalFraction_F*CardOut_F) %>% 
  mutate(Vrest_M = Vtotal_M - Vskin_M - Vkidney_M - Vgut_M - Vliver_M - Vplasma_M - Vadipose_M) %>%
  mutate(Qrest_M = Qtotal_M - Qskin_M - Qkidney_M - Qhepatic_M - Qgut_M - Qadipose_M) %>%
  mutate(Vrest_F = Vtotal_F - Vskin_F - Vkidney_F - Vgut_F - Vliver_F - Vplasma_F - Vadipose_F) %>%
  mutate(Qrest_F = Qtotal_F - Qskin_F - Qkidney_F - Qhepatic_F - Qgut_F - Qadipose_F) %>%
  
  ## MassBalance Flow #better be 0
  # !!!!! Issue with the mass balance of the flows
  mutate(Qmass_balance_M = CardOut_M - (Qskin_M + Qkidney_M + Qgut_M + Qhepatic_M + Qadipose_M + Qrest_M)) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  mutate(Qmass_balance_F = CardOut_F - (Qskin_F + Qkidney_F + Qgut_F + Qhepatic_F + Qadipose_F + Qrest_F)) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  
  ## MassBalance Volumes #better be 0
  # !!!!! Issue with the mass balance of the flows
  mutate(Vmass_balance_M = BW_M - (Vplasma_M + Vskin_M + Vkidney_M + Vgut_M + Vliver_M + Vadipose_M + Vrest_M)) %>% 
  mutate(Vmass_balance_F = BW_F - (Vplasma_F + Vskin_F + Vkidney_F + Vgut_F + Vliver_F + Vadipose_F + Vrest_F)) %>% 
  
  ## Biliary clearance
  # mutate(CLbiliaryFraction_M = 0.000109) %>% #from Husoy; 2.62 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day taken from Fujii et al 2015 
  # mutate(CLbiliary_M = CLbiliaryFraction_M*BDW_M^VliverFraction_M) %>% #from Husoy; L/h biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  # mutate(CLbiliaryFraction_F = 0.000109) %>% #from Husoy; 2.62/24/1000 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day taken from Fujii et al 2015 
  # mutate(CLbiliary_F = CLbiliaryFraction_M*BDW_F^VliverFraction_F) %>% #from Husoy; L/h biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  # Comment Chrysa 18-10-2024: Shouldn't we be correcting for the fraction of liver changing during life?
  mutate(CLbiliary_M = CLbiliaryc*(BW_M^0.1)) %>% #from Husoy; L/d biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  mutate(CLbiliary_F = CLbiliaryc*(BW_F^0.1)) %>% #from Husoy; L/d biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  
  ## Fecal clearance
  # mutate(CLfecalFraction_M = 0.00262) %>% #from Husoy; 0.052/24/1000 feces clearance L/h/kg clearance in feces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg
  # mutate(CLfecal_M = CLfecalFraction_M*BDW_M^VgutFraction_M) %>% # L/h faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  # mutate(CLfecalFraction_F = 0.000052) %>% #from Husoy; 0.052/24/1000 feces clearance L/h/kg clearance in feces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg
  # mutate(CLfecal_F = CLfecalFraction_F*BDW_F^VgutFraction_F) %>% # L/h faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  mutate(CLfecal_M = CLfaecesc*(BW_M^0.001)) %>% #from Husoy; L/d faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  mutate(CLfecal_F = CLfaecesc*(BW_F^0.001)) %>% #from Husoy; L/d faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  
  ## Urinary clearance 
  # mutate(CLurine_M = CLurinec*BW_M^(-0.25)) %>% # clearance urine (L/d) # NOT USED
  # mutate(CLurine_F = CLurinec*BW_F^(-0.25)) %>% # clearance urine (L/d) # NOT USED
  # mutate(MPT_M = Vkidney_M*1000*Kcells*Kprotein) %>%	# mass proximal tubule cells in gram based on BW
  # mutate(MPT_F = Vkidney_F*1000*Kcells*Kprotein) %>%	# mass proximal tubule cells in gram based on BW
  # mutate(Vmax_M = Vmaxc * MPT_M * SFOAT4) %>% #ug/d
  # mutate(Vmax_F = Vmaxc * MPT_F * SFOAT4) %>% #ug/d
  mutate(CL_PltPT = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTCPGK * KW_cortex) %>% #L/d plasma to proximal tubule clearance
  mutate(CL_FiltPT = (CL_OAT4*REF_OAT4) * PTCPGK * KW_cortex) #L/d filtrate to proximal tubule clearance 

## New Chrysa 
# Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
mutate(SkbTarea_M = 9.1*((BW_M*1000)^0.666)) %>%  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  mutate(SkbTarea_F = 9.1*((BW_F*1000)^0.666)) %>%  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  mutate(skin_fraction = 0.05) %>% # fraction of the total skin surface exposed; hands are 5% https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf
  # Note Chrysa 18-10-2024: the mean % of total surface area that is hands is 5.2 (in adult male 21+ years) and 4.8 (in adult female), so that would be 5 +/- 0.2 %
  # Note Chrysa 24-20-2024: assumption that the main skin area exposed is the hands
  mutate(fSkbarea_M = skin_fraction*SkbTarea_M) %>% # 1070 (cm^2); exposed skin area = surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
  mutate(fSkbarea_F = skin_fraction*SkbTarea_F) %>% # 1070 (cm^2); exposed skin area = surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
  mutate(Skbthickness_M = 83.1/1000) %>% # (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
  mutate(Skbthickness_F = 83.1/1000) %>% # (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
  mutate(VSkb_M = (fSkbarea_M*Skbthickness_M)/1000) %>%  #(L); Skin barrier volume; as previously coded by Trine
  mutate(VSkb_F = (fSkbarea_F*Skbthickness_F)/1000) %>%  #(L); Skin barrier volume; as previously coded by Trine 
  







# view(Variables_df)

##### Gender-specific data frames ----  
# Variables_M <- Variables_df %>%
#   select(TIME, age, BW_M, Hct_M, CardOut_M, Vplasma_M, Vart_M, Vven_M, Vliver_M, Vstomach_M, Vgut_M, Vkidney_M, Vurinarytract_M,
#          Vskin_M, Vadipose_M, Vadrenal_M, Vbone_M, Vbonenonperfused_M, Vbrain_M, Vbreast_M, Vheart_M, Vmarrow_M, Vmuscle_M,
#          Vrepro_M, Vpancreas_M, Vspleen_M, Vthyroid_M, Vlung_M,
#          Qliver_M, Qstomach_M, Qgut_M, Qkidney_M, Qurinarytract_M,
#          Qskin_M, Qadipose_M, Qadrenal_M, Qbone_M, Qbrain_M, Qbreast_M, Qheart_M, Qmarrow_M, Qmuscle_M,
#          Qrepro_M, Qpancreas_M, Qspleen_M, Qthyroid_M, Qlung_M) %>%
#   mutate(Vtotal_M = Vplasma_M + Vliver_M + Vstomach_M + Vgut_M + Vkidney_M + Vurinarytract_M +
#            Vskin_M + Vadipose_M + Vadrenal_M + Vbone_M + Vbonenonperfused_M + Vbrain_M + Vbreast_M + Vheart_M + Vmarrow_M + Vmuscle_M +
#            Vrepro_M + Vpancreas_M + Vspleen_M + Vthyroid_M + Vlung_M) %>%
#   mutate(Fraction_M = Vtotal_M/BW_M) %>%
#   mutate(Qtotal_M = Qliver_M + Qstomach_M + Qgut_M + Qkidney_M + Qurinarytract_M +
#            Qskin_M + Qadipose_M + Qadrenal_M + Qbone_M + Qbrain_M + Qbreast_M + Qheart_M + Qmarrow_M + Qmuscle_M +
#            Qrepro_M + Qpancreas_M + Qspleen_M + Qthyroid_M + Qlung_M)
# 
# Variables_fractions_M <- Variables_df %>%
#   select(TIME, age, VadrenalFraction_M, VboneFraction_M, VbonenonperfusedFraction_M, VbrainFraction_M, VbreastFraction_M, 
#          VheartFraction_M, VmarrowFraction_M, VmuscleFraction_M, VreproFraction_M, VpancreasFraction_M,
#          VskinFraction_M, VspleenFraction_M, VthyroidFraction_M, VurinarytractFraction_M, VkidneyFraction_M,
#          VlungFraction_M, VgutFraction_M, VstomachFraction_M, VliverFraction_M, VplasmaFraction_M, VadiposeFraction_M,
#          QadrenalFraction_M, QboneFraction_M, QbrainFraction_M, QbreastFraction_M, 
#          QheartFraction_M, QmarrowFraction_M, QmuscleFraction_M, QreproFraction_M, QpancreasFraction_M,
#          QskinFraction_M, QspleenFraction_M, QthyroidFraction_M, QurinarytractFraction_M, QkidneyFraction_M,
#          QlungFraction_M, QgutFraction_M, QstomachFraction_M, QliverFraction_M, QadiposeFraction_M) %>%
#   mutate(VtotalFraction_M = VadrenalFraction_M + VboneFraction_M + VbonenonperfusedFraction_M + VbrainFraction_M + VbreastFraction_M + 
#            VheartFraction_M + VmarrowFraction_M + VmuscleFraction_M + VreproFraction_M + VpancreasFraction_M +
#            VskinFraction_M + VspleenFraction_M + VthyroidFraction_M + VurinarytractFraction_M + VkidneyFraction_M +
#            VlungFraction_M + VgutFraction_M + VstomachFraction_M + VliverFraction_M + VplasmaFraction_M + VadiposeFraction_M) %>%
#   mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M + 
#            QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
#            QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
#            QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M)
# 
# Variables_F <- Variables_df %>%
#   select(TIME, age, BW_F, Hct_F, CardOut_F, Vplasma_F, Vart_F, Vven_F, Vliver_F, Vstomach_F, Vgut_F, Vkidney_F, Vurinarytract_F,
#          Vskin_F, Vadipose_F, Vadrenal_F, Vbone_F, Vbonenonperfused_F, Vbrain_F, Vbreast_F, Vheart_F, Vmarrow_F, Vmuscle_F,
#          Vrepro_F, Vpancreas_F, Vspleen_F, Vthyroid_F, Vlung_F,
#          Qliver_F, Qstomach_F, Qgut_F, Qkidney_F, Qurinarytract_F,
#          Qskin_F, Qadipose_F, Qadrenal_F, Qbone_F, Qbrain_F, Qbreast_F, Qheart_F, Qmarrow_F, Qmuscle_F,
#          Qrepro_F, Qpancreas_F, Qspleen_F, Qthyroid_F, Qlung_F) %>%
#   mutate(Vtotal_F = Vplasma_F + Vliver_F + Vstomach_F + Vgut_F + Vkidney_F + Vurinarytract_F +
#            Vskin_F + Vadipose_F + Vadrenal_F + Vbone_F + Vbonenonperfused_F + Vbrain_F + Vbreast_F + Vheart_F + Vmarrow_F + Vmuscle_F +
#            Vrepro_F + Vpancreas_F + Vspleen_F + Vthyroid_F + Vlung_F) %>%
#   mutate(Fraction_F = Vtotal_F/BW_F) %>%
#   mutate(Qtotal_F = Qliver_F + Qstomach_F + Qgut_F + Qkidney_F + Qurinarytract_F +
#            Qskin_F + Qadipose_F + Qadrenal_F + Qbone_F + Qbrain_F + Qbreast_F + Qheart_F + Qmarrow_F + Qmuscle_F +
#            Qrepro_F + Qpancreas_F + Qspleen_F + Qthyroid_F + Qlung_F)
# 
# Variables_fractions_F <- Variables_df %>%
#   select(TIME, age, VadrenalFraction_F, VboneFraction_F, VbonenonperfusedFraction_F, VbrainFraction_F, VbreastFraction_F, 
#          VheartFraction_F, VmarrowFraction_F, VmuscleFraction_F, VreproFraction_F, VpancreasFraction_F,
#          VskinFraction_F, VspleenFraction_F, VthyroidFraction_F, VurinarytractFraction_F, VkidneyFraction_F,
#          VlungFraction_F, VgutFraction_F, VstomachFraction_F, VliverFraction_F, VplasmaFraction_F, VadiposeFraction_F,
#          QadrenalFraction_F, QboneFraction_F, QbrainFraction_F, QbreastFraction_F, 
#          QheartFraction_F, QmarrowFraction_F, QmuscleFraction_F, QreproFraction_F, QpancreasFraction_F,
#          QskinFraction_F, QspleenFraction_F, QthyroidFraction_F, QurinarytractFraction_F, QkidneyFraction_F,
#          QlungFraction_F, QgutFraction_F, QstomachFraction_F, QliverFraction_F, QadiposeFraction_F) %>%
#   mutate(VtotalFraction_F = VadrenalFraction_F + VboneFraction_F + VbonenonperfusedFraction_F + VbrainFraction_F + VbreastFraction_F + 
#            VheartFraction_F + VmarrowFraction_F + VmuscleFraction_F + VreproFraction_F + VpancreasFraction_F +
#            VskinFraction_F + VspleenFraction_F + VthyroidFraction_F + VurinarytractFraction_F + VkidneyFraction_F +
#            VlungFraction_F + VgutFraction_F + VstomachFraction_F + VliverFraction_F + VplasmaFraction_F + VadiposeFraction_F) %>%
#   mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F + 
#            QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
#            QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
#            QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F)
# 
# Figure male organ volumes
# Figure_M <- ggplot() +
#   geom_line(data = Variables_M, aes(x=age, y=BW_M, color="01_BW")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vplasma_M, color="02_Vplasma")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vliver_M, color="03_Vliver")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vstomach_M, color="04_Vstomach")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vgut_M, color="05_Vgut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vkidney_M, color="06_Vkidney")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vskin_M, color="07_Vskin")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vadipose_M, color="08_Vadipose")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vbrain_M, color="09_Vbrain")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vmuscle_M, color="10_Vmuscle")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vtotal_M, color="11_Vtotal")) +
# 
#   scale_colour_manual(name='',
#                       values=c('01_BW'='black',
#                                '02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey'),
#                       labels=c('01_BW'='Body weight',
#                                '02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues')) +
# 
#   theme_bw() +
#   labs(title="Organ weight changes over time - males",x="\nAge (y)", y="Weight (kg)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_M
# 
# Figure male organ fractions 
# Figure_fractions_M <- ggplot() +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VplasmaFraction_M, color="02_Vplasma")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VliverFraction_M, color="03_Vliver")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VstomachFraction_M, color="04_Vstomach")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VgutFraction_M, color="05_Vgut")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VkidneyFraction_M, color="06_Vkidney")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VskinFraction_M, color="07_Vskin")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VadiposeFraction_M, color="08_Vadipose")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VbrainFraction_M, color="09_Vbrain")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VmuscleFraction_M, color="10_Vmuscle")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VtotalFraction_M, color="11_Vtotal")) +
#   geom_line(data = Variables_M, aes(x=age, y=Fraction_M, color="12_Fraction")) +
# 
#   scale_colour_manual(name='',
#                       values=c('02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey',
#                                '12_Fraction'='black'),
#                       labels=c('02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues',
#                                '12_Fraction'='Sum tissues / BW')) +
# 
#   theme_bw() +
#   labs(title="Organ fractional volume changes over time - males",x="\nAge (y)", y="Fraction of body weight\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_fractions_M
# 
# Figure male blood flows 
# Figure_Q_M <- ggplot() +
#   geom_line(data = Variables_M, aes(x=age, y=CardOut_M, color="01_CardOut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qliver_M, color="03_Qliver")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qstomach_M, color="04_Qstomach")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qgut_M, color="05_Qgut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qkidney_M, color="06_Qkidney")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qskin_M, color="07_Qskin")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qadipose_M, color="08_Qadipose")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qbrain_M, color="09_Qbrain")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qmuscle_M, color="10_Qmuscle")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qtotal_M, color="11_Qtotal")) +
# 
#   scale_colour_manual(name='',
#                       values=c('01_CardOut'='black',
#                                '03_Qliver'='darkred',
#                                '04_Qstomach'='indianred',
#                                '05_Qgut'='lightpink',
#                                '06_Qkidney'='indianred4',
#                                '07_Qskin'='darksalmon',
#                                '08_Qadipose'='khaki',
#                                '09_Qbrain'='mistyrose',
#                                '10_Qmuscle'='rosybrown1',
#                                '11_Qtotal'='grey'),
#                       labels=c('01_CardOut'='Cardiac output',
#                                '03_Qliver'='Liver',
#                                '04_Qstomach'='Stomach',
#                                '05_Qgut'='Gut',
#                                '06_Qkidney'='Kidney',
#                                '07_Qskin'='Skin',
#                                '08_Qadipose'='Adipose tissue',
#                                '09_Qbrain'='Brain',
#                                '10_Qmuscle'='Muscle tissue',
#                                '11_Qtotal'='Sum of all tissues')) +
# 
#   theme_bw() +
#   labs(title="Plasma flow changes over time - males",x="\nAge (y)", y="Plasma flow (L/d)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_Q_M
# 
# Figure female organ volumes 
# Figure_F <- ggplot() +
#   geom_line(data = Variables_F, aes(x=age, y=BW_F, color="01_BW")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vplasma_F, color="02_Vplasma")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vliver_F, color="03_Vliver")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vstomach_F, color="04_Vstomach")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vgut_F, color="05_Vgut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vkidney_F, color="06_Vkidney")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vskin_F, color="07_Vskin")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vadipose_F, color="08_Vadipose")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vbrain_F, color="09_Vbrain")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vmuscle_F, color="10_Vmuscle")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vtotal_F, color="11_Vtotal")) +
#   
#   scale_colour_manual(name='',
#                       values=c('01_BW'='black',
#                                '02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey'),
#                       labels=c('01_BW'='Body weight',
#                                '02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues')) +
#   
#   theme_bw() +
#   labs(title="Organ weight changes over time - females",x="\nAge (y)", y="Weight (kg)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_F
# 
# Figure female organ fractions 
# Figure_fractions_F <- ggplot() +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VplasmaFraction_F, color="02_Vplasma")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VliverFraction_F, color="03_Vliver")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VstomachFraction_F, color="04_Vstomach")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VgutFraction_F, color="05_Vgut")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VkidneyFraction_F, color="06_Vkidney")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VskinFraction_F, color="07_Vskin")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VadiposeFraction_F, color="08_Vadipose")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VbrainFraction_F, color="09_Vbrain")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VmuscleFraction_F, color="10_Vmuscle")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VtotalFraction_F, color="11_Vtotal")) +
#   geom_line(data = Variables_F, aes(x=age, y=Fraction_F, color="12_Fraction")) +
#   
#   scale_colour_manual(name='',
#                       values=c('02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey',
#                                '12_Fraction'='black'),
#                       labels=c('02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues',
#                                '12_Fraction'='Sum tissues / BW')) +
#   
#   theme_bw() +
#   labs(title="Organ fractional volume changes over time - females",x="\nAge (y)", y="Fraction of body weight\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_fractions_F
# 
# Figure female blood flows 
# Figure_Q_F <- ggplot() +
#   geom_line(data = Variables_F, aes(x=age, y=CardOut_F, color="01_CardOut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qliver_F, color="03_Qliver")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qstomach_F, color="04_Qstomach")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qgut_F, color="05_Qgut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qkidney_F, color="06_Qkidney")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qskin_F, color="07_Qskin")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qadipose_F, color="08_Qadipose")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qbrain_F, color="09_Qbrain")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qmuscle_F, color="10_Qmuscle")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qtotal_F, color="11_Qtotal")) +
#   
#   scale_colour_manual(name='',
#                       values=c('01_CardOut'='black',
#                                '03_Qliver'='darkred',
#                                '04_Qstomach'='indianred',
#                                '05_Qgut'='lightpink',
#                                '06_Qkidney'='indianred4',
#                                '07_Qskin'='darksalmon',
#                                '08_Qadipose'='khaki',
#                                '09_Qbrain'='mistyrose',
#                                '10_Qmuscle'='rosybrown1',
#                                '11_Qtotal'='grey'),
#                       labels=c('01_CardOut'='Cardiac output',
#                                '03_Qliver'='Liver',
#                                '04_Qstomach'='Stomach',
#                                '05_Qgut'='Gut',
#                                '06_Qkidney'='Kidney',
#                                '07_Qskin'='Skin',
#                                '08_Qadipose'='Adipose tissue',
#                                '09_Qbrain'='Brain',
#                                '10_Qmuscle'='Muscle tissue',
#                                '11_Qtotal'='Sum of all tissues')) +
#   
#   theme_bw() +
#   labs(title="Plasma flow changes over time - females",x="\nAge (y)", y="Plasma flow (L/d)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_Q_F
# 
# Plot mass balance organ flow  
# 
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# FlowMassBalance <- Variables_df %>% select(age, Qmass_balance_M, Qmass_balance_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "MB", Qmass_balance_M:Qmass_balance_F) %>% 
#   ggplot(aes(age, MB)) +
#   geom_path() +
#   facet_wrap(~Gender)
# FlowMassBalance
# 
# # Plot mass balance organ volumes
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# VolumeMassBalance <- Variables_df %>% select(age, Vmass_balance_M, Vmass_balance_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "MB", Vmass_balance_M:Vmass_balance_F) %>% 
#   ggplot(aes(age, MB)) +
#   geom_path() +
#   facet_wrap(~Gender)
# VolumeMassBalance
# 
# # Plot mass balance organ volumes
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# VolumeTest <- Variables_df %>% select(age, VtotalFraction_M, VtotalFraction_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "Volume", VtotalFraction_M:VtotalFraction_F) %>% 
#   ggplot(aes(age, Volume, color = Gender)) +
#   geom_path() 
# VolumeTest



# mutate(Oraldose_M = Oraldose*BW_M) %>%
# mutate(Oraldose_F = Oraldose*BW_F) %>%
# mutate(Oraldose_M = if_else(age <= exposure_stop,Oraldose_M,0)) %>%
# mutate(Oraldose_F = if_else(age <= exposure_stop,Oraldose_F,0)) %>%
# 
# mutate(Dermaldose_M = AbsPFOA*Dermconc*BW_M) %>% # dermal concentration is expressed as ug/kg BW/day
# mutate(Dermaldose_F = AbsPFOA*Dermconc*BW_F) %>% # dermal concentration is expressed as ug/kg BW/day
# mutate(Dermaldose_M = if_else(age <= exposure_stop,Dermaldose_M,0)) %>%
# mutate(Dermaldose_F = if_else(age <= exposure_stop,Dermaldose_F,0)) %>%


#### Plots ----
#Plot BW changes over time
#Comment Chrysa on 18-10-2024: MassBalance issue: should we then use the BDW term that is also changing during adulthood or are we OK with the massbalance issue in the total body volume?

# PlotBDW <- Variables_df %>% select(c(age, BDW_M_Ratier_2024, BDW_F_Ratier_2024)) %>%
#   rename(Male = BDW_M_Ratier_2024, Female = BDW_F_Ratier_2024) %>%
#   pivot_longer(names_to = "Gender", values_to = "BDW", Male:Female)
# 
# PlotBW <- Variables_df %>% select(c(age, BW_M_Ratier_2024, BW_F_Ratier_2024)) %>%
#   rename(Male = BW_M_Ratier_2024, Female = BW_F_Ratier_2024) %>%
#   pivot_longer(names_to = "Gender", values_to = "BW", Male:Female)
# ggplot() +
#   geom_path(data = PlotBDW, aes(age, BDW, colour = Gender, linetype = "BDW")) +
#   geom_path(data = PlotBW, aes(age, BW, colour = Gender, linetype = "BW")) +
#   labs(x = "Age",
#        y = "Values",
#        colour = "Gender",
#        linetype = "Bodyweight method") +
#   scale_linetype_manual(values = c("BDW" = "dashed", "BW" = "solid"))
# ggsave("BWovertime.png", dpi = 300)