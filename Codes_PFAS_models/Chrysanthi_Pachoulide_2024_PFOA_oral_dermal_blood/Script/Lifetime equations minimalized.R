# Version to be used for the PFAS PARC model

# Lifetime equations
rm(list=ls()) # to clear out the global environment

library(tidyverse)

#### Settings ----
TSTART <- 0
TSTOP <- 365*80 # days
DT <- 1
TIME <- seq(TSTART,TSTOP,by=DT)

Variables_df <- as.data.frame(list(TIME = TIME))
Variables_df <- Variables_df %>%
  mutate(age = TIME/365) 

Variables_df <- Variables_df %>%
  # BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
  mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                    -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))

#### Organ volumes ----
# Using Ratier et al. (2024) model
# Fraction of arterial blood, calculated from Filser 2000 p.43
Fr_art_blood = 0.0178 / (0.0178 + 0.0533)

# Hematocrit ### IN MAN ### From Supp mat of Brochot et al. 2019
# ?????????????????? not sure what is this ?????????????????????
Param1_M = 33.455469
Param2_M = 53.206039
Param3_M = 8.277945
Param4_M = 40.492556
Param5_M = 46.899695

b1_M = (Param4_M - Param1_M -(Param2_M - Param1_M) * exp(-Param3_M))/5
a1_M = Param4_M - 6*b1_M
b2_M = (Param5_M - Param4_M)/5
a2_M = Param4_M - 15*b2_M

# Hematocrit ### IN NON-PREGNANT WOMAN ### 
Param1_F = 32.617402
Param2_F = 53.188459
Param3_F = 7.699418
Param4_F = 37.531463
Param5_F = 40.055284

b1_F = (Param4_F - Param1_F - (Param2_F-Param1_F)*exp(-Param3_F))/2
a1_F = Param4_F - 3*b1_F
b2_F = (Param5_F - Param4_F)/7
a2_F = Param5_F - 10*b2_F


# Changes in body composition over time
Variables_df <- Variables_df %>%
  select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
  rename(BW_M = BW_M_Ratier_2024) %>%
  rename(BW_F = BW_F_Ratier_2024) %>%
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
  # Hematocrit
  # ?????????? Is this the hematocrit change over lifetime ????????????
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
  
  # Cardiac output
  mutate(CardOut_M = if_else(age < 33.37, (6.642 + (0.6 - 6.642)*exp(-0.1323*age))*(1-Hct_M),
                             (-0.000895*age^2 + 0.0607*age + 5.54)*(1-Hct_M))) %>%
  mutate(CardOut_F = if_else(age < 16.027, (7.734 + (0.6 - 7.734)*exp(-0.09747*age))*(1-Hct_F),
                             (0.000473*age^2 - 0.0782*age + 7.37)*(1-Hct_F))) %>%
  
  # Blood volume; compartment [22] in Ratier 2024
  mutate(VbloodFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                    0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(Vblood_M = VbloodFraction_M*BDW_M) %>%
  mutate(Vart_M = Vblood_M*Fr_art_blood*(1-Hct_M)) %>%
  mutate(Vven_M = Vblood_M*(1-Hct_M) - Vart_M) %>%
  mutate(VbloodFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                    if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                            0.065))) %>%
  mutate(Vblood_F = VbloodFraction_F*BDW_F) %>%
  mutate(Vart_F = Vblood_F*Fr_art_blood*(1-Hct_F)) %>%
  mutate(Vven_F = Vblood_F*(1-Hct_F) - Vart_F) %>%
  
  # Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(VliverFraction_M = 0.0247 + (0.0409 - 0.0247)*exp(-0.218*age)) %>%
  mutate(Vliver_M = VliverFraction_M*BDW_M) %>%
  mutate(QliverFraction_M = (VliverFraction_M/0.0247)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  mutate(VliverFraction_F = 0.0233 + (0.038 - 0.0233)*exp(-0.122*age)) %>%
  mutate(Vliver_F = VliverFraction_F*BDW_F) %>%
  mutate(QliverFraction_F = (VliverFraction_F/0.0233)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  
  # Stomach; compartment [17] in Ratier 2024
  mutate(VstomachFraction_M = 0.0021) %>%
  mutate(Vstomach_M = VstomachFraction_M*BDW_M) %>%
  mutate(QstomachFraction_M = (VstomachFraction_M/0.0021)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VstomachFraction_F = 0.0023) %>%
  mutate(Vstomach_F = VstomachFraction_F*BDW_F) %>%
  mutate(QstomachFraction_F = (VstomachFraction_F/0.0023)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Gut; compartment [16] in Ratier 2024
  mutate(VgutFraction_M = if_else(age < 16, -0.000082562*age^2 + 0.0013523*age + 0.01293,
                                  0.0140)) %>%
  mutate(Vgut_M = VgutFraction_M*BDW_M) %>%
  mutate(QgutFraction_M = (VgutFraction_M/0.0140)*0.144) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VgutFraction_F = if_else(age < 14.453301, -7.42E-5*age^2 + 1.28E-3*age + 1.30E-2,
                                  0.0160)) %>%
  mutate(Vgut_F = VgutFraction_F*BDW_F) %>%
  mutate(QgutFraction_F = (VgutFraction_F/0.0160)*0.165) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Kidney; compartment [14] in Ratier 2024
  mutate(VkidneyFraction_M = 0.0042 + (0.00767 - 0.0042)*exp(-0.206*age)) %>%
  mutate(Vkidney_M = VkidneyFraction_M*BDW_M) %>%
  mutate(QkidneyFraction_M = (VkidneyFraction_M/0.0042)*0.196) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VkidneyFraction_F = 0.0046 + (0.0071 - 0.0046)*exp(-0.221*age)) %>%
  mutate(Vkidney_F = VkidneyFraction_F*BDW_F) %>%
  mutate(QkidneyFraction_F = (VkidneyFraction_F/0.0046)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Filtrate; compartment [13] (urinary tract) in Ratier 2024
  mutate(VurinarytractFraction_M = 0.00104) %>%
  mutate(Vurinarytract_M = VurinarytractFraction_M*BDW_M) %>%
  mutate(QurinarytractFraction_M = (VurinarytractFraction_M/0.00104)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VurinarytractFraction_F = 0.0010) %>%
  mutate(Vurinarytract_F = VurinarytractFraction_F*BDW_F) %>%
  mutate(QurinarytractFraction_F = (VurinarytractFraction_F/0.0010)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Skin; compartment [10] in Ratier 2024
  mutate(VskinFraction_M = if_else(age < 20.01, -1.1706E-05*age^3 + 5.4130E-04*age^2 - 6.1966E-03*age + 4.6231E-02,
                                   0.0452)) %>%
  mutate(Vskin_M = VskinFraction_M*BDW_M) %>%
  mutate(QskinFraction_M = (VskinFraction_M/0.0452)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VskinFraction_F = if_else(age < 19.45, -7.8882E-06*age^3 + 4.0224E-04*age^2 - 5.2146E-03*age + 4.5605E-02,
                                   0.0383)) %>%
  mutate(Vskin_F = VskinFraction_F*BDW_F) %>%
  mutate(QskinFraction_F = (VskinFraction_F/0.0383)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VadrenalFraction_M = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(Vadrenal_M = VadrenalFraction_M*BDW_M) %>%
  mutate(QadrenalFraction_M = (VadrenalFraction_M/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VadrenalFraction_F = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(Vadrenal_F = VadrenalFraction_F*BDW_F) %>%
  mutate(QadrenalFraction_F = (VadrenalFraction_F/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VboneFraction_M = (0.313 + (0.506 - 0.313)*exp(-0.0907*age))*0.095) %>%
  mutate(VbonenonperfusedFraction_M = 0.095 - VboneFraction_M) %>%
  mutate(Vbone_M = VboneFraction_M*BDW_M/2) %>% # 2 is bone density
  mutate(Vbonenonperfused_M = VbonenonperfusedFraction_M*BDW_M/2) %>% # 2 is bone density
  mutate(QboneFraction_M = (VboneFraction_M/(0.095*0.32))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VboneFraction_F = (0.298 + (0.505 - 0.298)*exp(-0.0792*age))*0.085) %>%
  mutate(VbonenonperfusedFraction_F = 0.085 - VboneFraction_F) %>%
  mutate(Vbone_F = VboneFraction_F*BDW_F/2) %>% # 2 is bone density
  mutate(Vbonenonperfused_F = VbonenonperfusedFraction_F*BDW_F/2) %>% # 2 is bone density
  mutate(QboneFraction_F = (VboneFraction_F/(0.085*0.298))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VbrainFraction_M = (1.450 + (0.353 - 1.450) * exp (-0.440*age))/BDW_M) %>%
  mutate(Vbrain_M = VbrainFraction_M*BDW_M) %>%
  mutate(QbrainFraction_M = (VbrainFraction_M/0.01986)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VbrainFraction_F = (1.300 + (0.347 - 1.300) * exp (-0.573*age))/BDW_F) %>%
  mutate(Vbrain_F = VbrainFraction_F*BDW_F) %>%
  mutate(QbrainFraction_F = (VbrainFraction_F/0.0217)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VbreastFraction_M = 3.42E-4*1/(1 + exp(-1.42*age + 20.1))) %>%
  mutate(Vbreast_M = VbreastFraction_M*BDW_M) %>%
  mutate(QbreastFraction_M = (VbreastFraction_M/0.00035)*0.0002) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VbreastFraction_F = 0.00833/(1 + exp(-1.92*age+ 28.6))) %>%
  mutate(Vbreast_F = VbreastFraction_F*BDW_F) %>%
  mutate(QbreastFraction_F = (VbreastFraction_F/0.0083)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VheartFraction_M = 0.0045) %>%
  mutate(Vheart_M = VheartFraction_M*BDW_M) %>%
  mutate(QheartFraction_M = (VheartFraction_M/0.0045)*0.041) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VheartFraction_F = 0.0042) %>%
  mutate(Vheart_F = VheartFraction_F*BDW_F) %>%
  mutate(QheartFraction_F = (VheartFraction_F/0.0042)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VmarrowFraction_M = 0.05 + (0.0138 - 0.05)*exp(-0.112*age)) %>%
  mutate(Vmarrow_M = VmarrowFraction_M*BDW_M) %>%
  mutate(QmarrowFraction_M = (VmarrowFraction_M/0.050)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VmarrowFraction_F = 0.045 + (0.0138 - 0.045)*exp(-0.136*age)) %>%
  mutate(Vmarrow_F = VmarrowFraction_F*BDW_F) %>%
  mutate(QmarrowFraction_F = (VmarrowFraction_F/0.045)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(MuscleAtrophy_M = if_else(age < 24.3, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(VmuscleFraction_M = (0.3973 + (0.201 - 0.3973)*exp(-0.141*age)) * MuscleAtrophy_M) %>%
  mutate(Vmuscle_M = VmuscleFraction_M*BDW_M) %>%
  mutate(QmuscleFraction_M = (VmuscleFraction_M/0.3973)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(MuscleAtrophy_F = if_else(age < 25.90709, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(VmuscleFraction_F = (0.2917 + (0.207 - 0.2917)*exp(-0.339*age)) * MuscleAtrophy_F) %>%
  mutate(Vmuscle_F = VmuscleFraction_F*BDW_F) %>%
  mutate(QmuscleFraction_F = (VmuscleFraction_F/0.2917)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VreproFraction_M = if_else(age < 20.01, -1.5156E-07*age^3 + 9.3351E-06*age^2 - 1.1177E-04*age + 4.7966E-04,
                                    0.0008)) %>%
  mutate(Vrepro_M = VreproFraction_M*BDW_M) %>%
  mutate(QreproFraction_M = (VreproFraction_M/0.0008)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VreproFraction_F = if_else(age < 1, -1.064E-3*age + 1.338E-3,
                                    if_else(age < 20, 2.6380E-7*age^3 - 1.7943E-6*age^2 - 5.6465E-6*age + 2.8105E-4,
                                            0.001552))) %>%
  mutate(Vrepro_F = VreproFraction_F*BDW_F) %>%
  mutate(QreproFraction_F = (VreproFraction_F/0.0016)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VpancreasFraction_M = 0.00192) %>%
  mutate(Vpancreas_M = VpancreasFraction_M*BDW_M) %>%
  mutate(QpancreasFraction_M = (VpancreasFraction_M/0.00192)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VpancreasFraction_F = 0.002) %>%
  mutate(Vpancreas_F = VpancreasFraction_F*BDW_F) %>%
  mutate(QpancreasFraction_F = (VpancreasFraction_F/0.002)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VspleenFraction_M = 0.0021) %>%
  mutate(Vspleen_M = VspleenFraction_M*BDW_M) %>%
  mutate(QspleenFraction_M = (VspleenFraction_M/0.0021)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VspleenFraction_F = 0.0022) %>%
  mutate(Vspleen_F = VspleenFraction_F*BDW_F) %>%
  mutate(QspleenFraction_F = (VspleenFraction_F/0.0022)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VthyroidFraction_M = 0.000274) %>%
  mutate(Vthyroid_M = VthyroidFraction_M*BDW_M) %>%
  mutate(QthyroidFraction_M = (VthyroidFraction_M/0.000274)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VthyroidFraction_F = 0.0003) %>%
  mutate(Vthyroid_F = VthyroidFraction_F*BDW_F) %>%
  mutate(QthyroidFraction_F = (VthyroidFraction_F/0.0003)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VlungFraction_M = 0.0068) %>%
  mutate(Vlung_M = VlungFraction_M*BDW_M) %>%
  mutate(QlungFraction_M = (VlungFraction_M/0.0068)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VlungFraction_F = 0.0070) %>%
  mutate(Vlung_F = VlungFraction_F*BDW_F) %>%
  mutate(QlungFraction_F = (VlungFraction_F/0.0070)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Adipose tissue
  mutate(VadiposeFraction_M = 0.96 - VadrenalFraction_M - VboneFraction_M - VbonenonperfusedFraction_M - VbrainFraction_M - VbreastFraction_M - 
           VheartFraction_M - VmarrowFraction_M - VmuscleFraction_M - VreproFraction_M - VpancreasFraction_M -
           VskinFraction_M - VspleenFraction_M - VthyroidFraction_M - VurinarytractFraction_M - VkidneyFraction_M -
           VlungFraction_M - VgutFraction_M - VstomachFraction_M - VliverFraction_M - VbloodFraction_M) %>%
  mutate(AdiposeMass_M = if_else(age < 19.00093277, 0,
                                 (-0.01129273*age^2 + 1.11817056*age + 56.74397436)-(74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))))) %>%
  mutate(Vadipose_M = (AdiposeMass_M/0.9) + VadiposeFraction_M*BDW_M/0.9) %>% # 0.9 is adipose tissue density
  mutate(QadiposeFraction_M = (VadiposeFraction_M/0.20)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VadiposeFraction_F = 0.96 - VadrenalFraction_F - VboneFraction_F - VbonenonperfusedFraction_F - VbrainFraction_F - VbreastFraction_F - 
           VheartFraction_F - VmarrowFraction_F - VmuscleFraction_F - VreproFraction_F - VpancreasFraction_F -
           VskinFraction_F - VspleenFraction_F - VthyroidFraction_F - VurinarytractFraction_F - VkidneyFraction_F -
           VlungFraction_F - VgutFraction_F - VstomachFraction_F - VliverFraction_F - VbloodFraction_F) %>%
  mutate(AdiposeMass_F = if_else(age < 17.9374115, 0,
                                 if_else(((-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))) < 0, 0,
                                         (-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))))) %>%
  mutate(Vadipose_F = (AdiposeMass_F/0.9) + VadiposeFraction_F*BDW_F/0.9) %>% # 0.9 is adipose tissue density
  mutate(QadiposeFraction_F = (VadiposeFraction_F/0.3167)*0.087) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M + 
           QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
           QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
           QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M) %>%
  mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F + 
           QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
           QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
           QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F) %>%
  
  mutate(Qliver_M = QliverFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qstomach_M = QstomachFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_M = QgutFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_M = QkidneyFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_M = QurinarytractFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_M = QskinFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_M = QadrenalFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_M = QboneFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_M = QbrainFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_M = QbreastFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_M = QheartFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_M = QmarrowFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_M = QmuscleFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_M = QreproFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_M = QpancreasFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_M = QspleenFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_M = QthyroidFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_M = QlungFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_M = QadiposeFraction_M/QtotalFraction_M*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  
  mutate(Qliver_F = QliverFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qstomach_F = QstomachFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_F = QgutFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_F = QkidneyFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_F = QurinarytractFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_F = QskinFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_F = QadrenalFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_F = QboneFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_F = QbrainFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_F = QbreastFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_F = QheartFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_F = QmarrowFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_F = QmuscleFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_F = QreproFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_F = QpancreasFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_F = QspleenFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_F = QthyroidFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_F = QlungFraction_F/QtotalFraction_F*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_F = QadiposeFraction_F/QtotalFraction_F*CardOut_F) # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;

#### Functions for inclusion in model code ----
#### Is this where the BW, V and Qs are changed over time?, 
# Volumes - Male
varBW_M <- approxfun(Variables_df$TIME, Variables_df$BW_M, rule = 2)
# so varX <- changed(in time(dependent), based on the body weight(independent), Constant?! )
varVblood_M <- approxfun(Variables_df$TIME, Variables_df$Vblood_M, rule = 2)
varVliver_M <- approxfun(Variables_df$TIME, Variables_df$Vliver_M, rule = 2)
varVstomach_M <- approxfun(Variables_df$TIME, Variables_df$Vstomach_M, rule = 2)
varVgut_M <- approxfun(Variables_df$TIME, Variables_df$Vgut_M, rule = 2)
varVkidney_M <- approxfun(Variables_df$TIME, Variables_df$Vkidney_M, rule = 2)
varVskin_M <- approxfun(Variables_df$TIME, Variables_df$Vskin_M, rule = 2)
varVadipose_M <- approxfun(Variables_df$TIME, Variables_df$Vadipose_M, rule = 2)
varVbrain_M <- approxfun(Variables_df$TIME, Variables_df$Vbrain_M, rule = 2)
varVmuscle_M <- approxfun(Variables_df$TIME, Variables_df$Vmuscle_M, rule = 2)
varVtotal_M <- approxfun(Variables_df$TIME, Variables_df$Vtotal_M, rule = 2)
varVurinarytract_M <- approxfun(Variables_df$TIME, Variables_df$Vurinarytract_M, rule = 2)
varVart_M <- approxfun(Variables_df$TIME, Variables_df$Vart_M, rule = 2)
varVven_M <- approxfun(Variables_df$TIME, Variables_df$Vven_M, rule = 2)
varVadrenal_M <- approxfun(Variables_df$TIME, Variables_df$Vadrenal_M, rule = 2)
varVbone_M <- approxfun(Variables_df$TIME, Variables_df$Vbone_M, rule = 2)
varVbonenonperfused_M <- approxfun(Variables_df$TIME, Variables_df$Vbonenonperfused_M, rule = 2)
varVbreast_M <- approxfun(Variables_df$TIME, Variables_df$Vbreast_M, rule = 2)
varVheart_M <- approxfun(Variables_df$TIME, Variables_df$Vheart_M, rule = 2)
varVmarrow_M <- approxfun(Variables_df$TIME, Variables_df$Vmarrow_M, rule = 2)
varVrepro_M <- approxfun(Variables_df$TIME, Variables_df$Vrepro_M, rule = 2)
varVpancreas_M <- approxfun(Variables_df$TIME, Variables_df$Vpancreas_M, rule = 2)
varVspleen_M <- approxfun(Variables_df$TIME, Variables_df$Vspleen_M, rule = 2)
varVthyroid_M <- approxfun(Variables_df$TIME, Variables_df$Vthyroid_M, rule = 2)
varVlung_M <- approxfun(Variables_df$TIME, Variables_df$Vlung_M, rule = 2)
# Blood flows - Male
varCardOut_M <- approxfun(Variables_df$TIME, Variables_df$CardOut_M, rule = 2)
varQliver_M <- approxfun(Variables_df$TIME, Variables_df$Qliver_M, rule = 2)
varQstomach_M <- approxfun(Variables_df$TIME, Variables_df$Qstomach_M, rule = 2)
varQgut_M <- approxfun(Variables_df$TIME, Variables_df$Qgut_M, rule = 2)
varQkidney_M <- approxfun(Variables_df$TIME, Variables_df$Qkidney_M, rule = 2)
varQskin_M <- approxfun(Variables_df$TIME, Variables_df$Qskin_M, rule = 2)
varQadipose_M <- approxfun(Variables_df$TIME, Variables_df$Qadipose_M, rule = 2)
varQbrain_M <- approxfun(Variables_df$TIME, Variables_df$Qbrain_M, rule = 2)
varQmuscle_M <- approxfun(Variables_df$TIME, Variables_df$Qmuscle_M, rule = 2)
varQtotal_M <- approxfun(Variables_df$TIME, Variables_df$Qtotal_M, rule = 2)
varQurinarytract_M <- approxfun(Variables_df$TIME, Variables_df$Qurinarytract_M, rule = 2)
varQadrenal_M <- approxfun(Variables_df$TIME, Variables_df$Qadrenal_M, rule = 2)
varQbone_M <- approxfun(Variables_df$TIME, Variables_df$Qbone_M, rule = 2)
varQbreast_M <- approxfun(Variables_df$TIME, Variables_df$Qbreast_M, rule = 2)
varQheart_M <- approxfun(Variables_df$TIME, Variables_df$Qheart_M, rule = 2)
varQmarrow_M <- approxfun(Variables_df$TIME, Variables_df$Qmarrow_M, rule = 2)
varQrepro_M <- approxfun(Variables_df$TIME, Variables_df$Qrepro_M, rule = 2)
varQpancreas_M <- approxfun(Variables_df$TIME, Variables_df$Qpancreas_M, rule = 2)
varQspleen_M <- approxfun(Variables_df$TIME, Variables_df$Qspleen_M, rule = 2)
varQthyroid_M <- approxfun(Variables_df$TIME, Variables_df$Qthyroid_M, rule = 2)
varQlung_M <- approxfun(Variables_df$TIME, Variables_df$Qlung_M, rule = 2)

# Volumes - Female
varBW_F <- approxfun(Variables_df$TIME, Variables_df$BW_F, rule = 2)
varVblood_F <- approxfun(Variables_df$TIME, Variables_df$Vblood_F, rule = 2)
varVliver_F <- approxfun(Variables_df$TIME, Variables_df$Vliver_F, rule = 2)
varVstomach_F <- approxfun(Variables_df$TIME, Variables_df$Vstomach_F, rule = 2)
varVgut_F <- approxfun(Variables_df$TIME, Variables_df$Vgut_F, rule = 2)
varVkidney_F <- approxfun(Variables_df$TIME, Variables_df$Vkidney_F, rule = 2)
varVskin_F <- approxfun(Variables_df$TIME, Variables_df$Vskin_F, rule = 2)
varVadipose_F <- approxfun(Variables_df$TIME, Variables_df$Vadipose_F, rule = 2)
varVbrain_F <- approxfun(Variables_df$TIME, Variables_df$Vbrain_F, rule = 2)
varVmuscle_F <- approxfun(Variables_df$TIME, Variables_df$Vmuscle_F, rule = 2)
varVtotal_F <- approxfun(Variables_df$TIME, Variables_df$Vtotal_F, rule = 2)
varVurinarytract_F <- approxfun(Variables_df$TIME, Variables_df$Vurinarytract_F, rule = 2)
varVart_F <- approxfun(Variables_df$TIME, Variables_df$Vart_F, rule = 2)
varVven_F <- approxfun(Variables_df$TIME, Variables_df$Vven_F, rule = 2)
varVadrenal_F <- approxfun(Variables_df$TIME, Variables_df$Vadrenal_F, rule = 2)
varVbone_F <- approxfun(Variables_df$TIME, Variables_df$Vbone_F, rule = 2)
varVbonenonperfused_F <- approxfun(Variables_df$TIME, Variables_df$Vbonenonperfused_F, rule = 2)
varVbreast_F <- approxfun(Variables_df$TIME, Variables_df$Vbreast_F, rule = 2)
varVheart_F <- approxfun(Variables_df$TIME, Variables_df$Vheart_F, rule = 2)
varVmarrow_F <- approxfun(Variables_df$TIME, Variables_df$Vmarrow_F, rule = 2)
varVrepro_F <- approxfun(Variables_df$TIME, Variables_df$Vrepro_F, rule = 2)
varVpancreas_F <- approxfun(Variables_df$TIME, Variables_df$Vpancreas_F, rule = 2)
varVspleen_F <- approxfun(Variables_df$TIME, Variables_df$Vspleen_F, rule = 2)
varVthyroid_F <- approxfun(Variables_df$TIME, Variables_df$Vthyroid_F, rule = 2)
varVlung_F <- approxfun(Variables_df$TIME, Variables_df$Vlung_F, rule = 2)
# Blood flows - Female
varCardOut_F <- approxfun(Variables_df$TIME, Variables_df$CardOut_F, rule = 2)
varQliver_F <- approxfun(Variables_df$TIME, Variables_df$Qliver_F, rule = 2)
varQstomach_F <- approxfun(Variables_df$TIME, Variables_df$Qstomach_F, rule = 2)
varQgut_F <- approxfun(Variables_df$TIME, Variables_df$Qgut_F, rule = 2)
varQkidney_F <- approxfun(Variables_df$TIME, Variables_df$Qkidney_F, rule = 2)
varQskin_F <- approxfun(Variables_df$TIME, Variables_df$Qskin_F, rule = 2)
varQadipose_F <- approxfun(Variables_df$TIME, Variables_df$Qadipose_F, rule = 2)
varQbrain_F <- approxfun(Variables_df$TIME, Variables_df$Qbrain_F, rule = 2)
varQmuscle_F <- approxfun(Variables_df$TIME, Variables_df$Qmuscle_F, rule = 2)
varQtotal_F <- approxfun(Variables_df$TIME, Variables_df$Qtotal_F, rule = 2)
varQurinarytract_F <- approxfun(Variables_df$TIME, Variables_df$Qurinarytract_F, rule = 2)
varQadrenal_F <- approxfun(Variables_df$TIME, Variables_df$Qadrenal_F, rule = 2)
varQbone_F <- approxfun(Variables_df$TIME, Variables_df$Qbone_F, rule = 2)
varQbreast_F <- approxfun(Variables_df$TIME, Variables_df$Qbreast_F, rule = 2)
varQheart_F <- approxfun(Variables_df$TIME, Variables_df$Qheart_F, rule = 2)
varQmarrow_F <- approxfun(Variables_df$TIME, Variables_df$Qmarrow_F, rule = 2)
varQrepro_F <- approxfun(Variables_df$TIME, Variables_df$Qrepro_F, rule = 2)
varQpancreas_F <- approxfun(Variables_df$TIME, Variables_df$Qpancreas_F, rule = 2)
varQspleen_F <- approxfun(Variables_df$TIME, Variables_df$Qspleen_F, rule = 2)
varQthyroid_F <- approxfun(Variables_df$TIME, Variables_df$Qthyroid_F, rule = 2)
varQlung_F <- approxfun(Variables_df$TIME, Variables_df$Qlung_F, rule = 2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Model for PFOA
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Comment JW: The first part of the PBK model description should include the time-dependent variables
# Please make sure that all relevant compartments are represented

PBPKmodPFOA <- function(t,state,parameters){
  with(as.list(c(state,parameters)), {
    QCP <- varCardOut_M(t) #finds the value at TIME == t? 
    QG <- varQgut_M(t) # needs to check what G includes (is this only intestine?)
    QL <- varQliver_M(t)
    QF <- varQadipose_M(t)
    QK <- varQkidney_M(t)
    Qfil <- varQurinarytract_M(t)
    QSk <- varQskin_M(t)
    #QR <- varQrest(t) # needs to be calculated based on the lumped compartments: 
    #Qp <- varQp(t) # needs to be incorporated when inhalation exposure is included
    Vart <- varVart_M(t)
    Vven <- varVven_M(t)
    VG <- varVgut_M(t)
    VL <- varVliver_M(t)
    VF <- varVadipose_M(t)
    VK <- varVkidney_M(t)
    Vfil <- varVurinarytract_M(t)
    VSk <- varVskin_M(t)
    
    if(t<tchng){DoseOn=1} else{DoseOn=0}
    
    Input1 <- Oraldose/Tinput*(t %% tinterval<Tinput) 
    Input2 <- Dermdose/Tinput*(t %% tinterval<Tinput)
    
