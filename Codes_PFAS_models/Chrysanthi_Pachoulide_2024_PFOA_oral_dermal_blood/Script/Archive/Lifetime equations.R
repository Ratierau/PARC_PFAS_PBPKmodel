# Lifetime equations
rm(list=ls()) # to clear out the global environment

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PFAS_PARC"
setwd(HOME)

library(deSolve)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#### Settings ----
TSTART <- 0
TSTOP <- 365*80 # days
DT <- 1
TIME <- seq(TSTART,TSTOP,by=DT)

# Data from EFSA 2012 [EFSA Journal 2012;10(3):2579]
# list.files(path=getwd())
# file <- file.choose()
file <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PFAS_PARC/Input/EFSA 2012 BW dataset.csv"
EFSA_2012 <- read.csv2(file,na.strings=c(""," ","NA"),stringsAsFactors = FALSE)
EFSA_2012[] <- lapply(EFSA_2012, as.numeric)

Variables_df <- as.data.frame(list(TIME = TIME))
Variables_df <- Variables_df %>%
  mutate(age = TIME/365) %>%
  mutate(age_m = TIME/(29200/960)) # %>%
# mutate(day = TIME+1) %>%
# mutate(month = ceiling(day/(365/12))-(((ceiling(day/365))-1)*12)) %>%
# mutate(year = ceiling(day/(365)))

Variables_df <- Variables_df %>%
  # BW_M & BW_F = Manual polynomial fit by JW of EFSA 2012 data, forced through a BWbirth of 3.5 kg
  mutate(BW_M = if_else(age <= 16, -0.00336521*age^4 + 0.116*age^3 - 1.17642062*age^2 + 6.44121320*age + 3.5,
                        0.0000223806*age^3 - 0.01881935*age^2 + 1.8946701*age + 34.41136163)) %>%
  mutate(BW_F = if_else(age <= 16, -0.00336521*age^4 + 0.116*age^3 - 1.17642062*age^2 + 6.44121320*age + 3.5,
                        -0.00024183*age^3 + 0.02946050*age^2 -0.81218974*age + 66.44367325)) %>%
  # BW_M_WG & BW_F_WG = Polynomial equation as used in BodyParameters file from WG parameters
  # mutate(BW_M_WG = 8.917 - 3.499e-02 * age_m + 4.152e-03 * (age_m)^2 - 1.917e-05 * (age_m)^3 +
  #          3.621e-08 * (age_m)^4 - 3.126e-11 * (age_m)^5 + 1.019e-14 * (age_m)^6) %>%
  # mutate(BW_F_WG = 6.564 + 1.193e-01 * age_m + 2.412e-03 * (age_m)^2 - 1.357e-05 * (age_m)^3 +
  #          2.832e-08 * (age_m)^4 - 2.630e-11 * (age_m)^5 + 9.097e-15 * (age_m)^6) %>%
  mutate(BW_M_WG = 3.938425 + 7.518199e-01* age_m - 2.023793e-02 * (age_m)^2 + 2.921682e-04 * (age_m)^3 - 
           2.06762e-06 * (age_m)^4 + 8.469e-09 * (age_m)^5 - 2.188427e-11 * (age_m)^6 + 
           3.699776e-14 * (age_m)^7 - 4.099077e-17 * (age_m)^8 + 2.874804e-20 * (age_m)^9 - 
           1.159732e-23 * (age_m)^10 + 2.052602e-27 * (age_m)^11) %>%
  mutate(BW_F_WG = 3.932403 + 6.866462e-01 * age_m - 1.949911e-02 * (age_m)^2 + 3.1311e-04 * (age_m)^3 - 
           2.466654e-06 * (age_m)^4 + 1.113217e-08 * (age_m)^5 - 3.131402e-11 * (age_m)^6 + 
           5.693737e-14 * (age_m)^7 - 6.706947e-17 * (age_m)^8 + 4.947858e-20 * (age_m)^9 - 
           2.079251e-23 * (age_m)^10 + 3.800367e-27 * (age_m)^11) %>%
  
  # BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
  mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                    -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))

BW_df <- Variables_df
#### Figure body weights ----
Figure_BW <- ggplot() +
  # geom_line(data=BW_df, aes(x=age, y=BW_M, color="01_BW_M", lty="01_BW_M")) +
  # geom_line(data=BW_df, aes(x=age, y=BW_F, color="02_BW_F", lty="02_BW_F")) +
  geom_line(data=BW_df, aes(x=age, y=BW_M_WG, color="03_BW_M_WG", lty="03_BW_M_WG")) +
  geom_line(data=BW_df, aes(x=age, y=BW_F_WG, color="04_BW_F_WG", lty="04_BW_F_WG")) +
  geom_line(data=BW_df, aes(x=age, y=BW_M_Ratier_2024, color="05_BW_M_Ratier", lty="05_BW_M_Ratier")) +
  geom_line(data=BW_df, aes(x=age, y=BW_F_Ratier_2024, color="06_BW_F_Ratier", lty="06_BW_F_Ratier")) +
  # geom_line(data=BW_df, aes(x=age, y=BDW_M_Ratier_2024, color="07_BDW_M_Ratier", lty="07_BDW_M_Ratier")) +
  # geom_line(data=BW_df, aes(x=age, y=BDW_F_Ratier_2024, color="08_BDW_F_Ratier", lty="08_BDW_F_Ratier")) +
  geom_point(data=EFSA_2012, aes(x=Age,y=BW_M_EFSA, color="09_BW_M_EFSA",pch="09_BW_M_EFSA")) +
  geom_point(data=EFSA_2012, aes(x=Age,y=BW_F_EFSA, color="10_BW_F_EFSA",pch="10_BW_F_EFSA")) +
  geom_ribbon(data=EFSA_2012, aes(ymin = BW_M_EFSA_P5, ymax = BW_M_EFSA_P95, x = Age), fill = "cadetblue1", alpha = 0.3) +
  geom_ribbon(data=EFSA_2012, aes(ymin = BW_F_EFSA_P5, ymax = BW_F_EFSA_P95, x = Age), fill = "lightpink1", alpha = 0.3) +
  
  scale_colour_manual(name='',
                      values=c(#'01_BW_M'='cadetblue1',
                        #'02_BW_F'='lightpink1',
                        '03_BW_M_WG'='dodgerblue',
                        '04_BW_F_WG'='violet',
                        '05_BW_M_Ratier'='dodgerblue4',
                        '06_BW_F_Ratier'='violetred4',
                        # '07_BDW_M_Ratier'='blue',
                        # '08_BDW_F_Ratier'='red',
                        '09_BW_M_EFSA'='cadetblue1',
                        '10_BW_F_EFSA'='lightpink1'),
                      labels=c('01_BW_M'='Body weight EFSA males',
                               '02_BW_F'='Body weight EFSA females',
                               '03_BW_M_WG'='Body weight males WG parameters equation',
                               '04_BW_F_WG'='Body weight females WG parameters equation',
                               '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                               '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                               # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                               # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                               '09_BW_M_EFSA'='Body weight data EFSA males',
                               '10_BW_F_EFSA'='Body weight data EFSA females'),
                      guide="none") +
  
  scale_linetype_manual(name='Predicted',
                        values=c(#'01_BW_M'='solid',
                          #'02_BW_F'='solid',
                          '03_BW_M_WG'='solid',
                          '04_BW_F_WG'='solid',
                          '05_BW_M_Ratier'='solid',
                          '06_BW_F_Ratier'='solid',
                          # '07_BDW_M_Ratier'='solid',
                          # '08_BDW_F_Ratier'='solid',
                          '09_BW_M_EFSA'=NA,
                          '10_BW_F_EFSA'=NA),
                        labels=c(#'01_BW_M'='Body weight EFSA males',
                          #'02_BW_F'='Body weight EFSA females',
                          '03_BW_M_WG'='Body weight males WG parameters equation',
                          '04_BW_F_WG'='Body weight females WG parameters equation',
                          '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                          '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                          # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                          # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                          '09_BW_M_EFSA'='Body weight data EFSA males',
                          '10_BW_F_EFSA'='Body weight data EFSA females')) +
  guides(linetype = guide_legend(override.aes = list(color = c("dodgerblue","violet","dodgerblue4","violetred4")))) + #"cadetblue1","lightpink1", ... ,"blue","red")))) +
  
  scale_shape_manual(name='Observed',
                     values=c(#'01_BW_M'=NA,
                       #'02_BW_F'=NA,
                       '03_BW_M_WG'=NA,
                       '04_BW_F_WG'=NA,
                       '05_BW_M_Ratier'=NA,
                       '06_BW_F_Ratier'=NA,
                       # '07_BDW_M_Ratier'=NA,
                       # '08_BDW_F_Ratier'=NA,
                       '09_BW_M_EFSA'=15,
                       '10_BW_F_EFSA'=16),
                     labels=c(#'01_BW_M'='Body weight EFSA males',
                       #'02_BW_F'='Body weight EFSA females',
                       '03_BW_M_WG'='Body weight males WG parameters equation',
                       '04_BW_F_WG'='Body weight females WG parameters equation',
                       '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                       '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                       # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                       # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                       '09_BW_M_EFSA'='Body weight data EFSA males',
                       '10_BW_F_EFSA'='Body weight data EFSA females')) +
  guides(shape = guide_legend(override.aes = list(color = c("cadetblue1","lightpink1")))) + 
  
  theme_bw() +
  labs(title="Body weight changes over time",x="\nAge (y)", y="Body weight (kg)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_BW

BW_df_zoom <- BW_df %>%
  filter(age <= 2)
EFSA_2012_zoom <- EFSA_2012 %>%
  filter(Age <= 2)

Figure_BW_zoom <- ggplot() +
  # geom_line(data=BW_df_zoom, aes(x=age, y=BW_M, color="01_BW_M", lty="01_BW_M")) +
  # geom_line(data=BW_df_zoom, aes(x=age, y=BW_F, color="02_BW_F", lty="02_BW_F")) +
  geom_line(data=BW_df_zoom, aes(x=age, y=BW_M_WG, color="03_BW_M_WG", lty="03_BW_M_WG")) +
  geom_line(data=BW_df_zoom, aes(x=age, y=BW_F_WG, color="04_BW_F_WG", lty="04_BW_F_WG")) +
  geom_line(data=BW_df_zoom, aes(x=age, y=BW_M_Ratier_2024, color="05_BW_M_Ratier", lty="05_BW_M_Ratier")) +
  geom_line(data=BW_df_zoom, aes(x=age, y=BW_F_Ratier_2024, color="06_BW_F_Ratier", lty="06_BW_F_Ratier")) +
  # geom_line(data=BW_df_zoom, aes(x=age, y=BDW_M_Ratier_2024, color="07_BDW_M_Ratier", lty="07_BDW_M_Ratier")) +
  # geom_line(data=BW_df_zoom, aes(x=age, y=BDW_F_Ratier_2024, color="08_BDW_F_Ratier", lty="08_BDW_F_Ratier")) +
  geom_point(data=EFSA_2012_zoom, aes(x=Age,y=BW_M_EFSA, color="09_BW_M_EFSA",pch="09_BW_M_EFSA")) +
  geom_point(data=EFSA_2012_zoom, aes(x=Age,y=BW_F_EFSA, color="10_BW_F_EFSA",pch="10_BW_F_EFSA")) +
  geom_ribbon(data=EFSA_2012_zoom, aes(ymin = BW_M_EFSA_P5, ymax = BW_M_EFSA_P95, x = Age), fill = "cadetblue1", alpha = 0.3) +
  geom_ribbon(data=EFSA_2012_zoom, aes(ymin = BW_F_EFSA_P5, ymax = BW_F_EFSA_P95, x = Age), fill = "lightpink1", alpha = 0.3) +
  
  scale_colour_manual(name='',
                      values=c(#'01_BW_M'='cadetblue1',
                        #'02_BW_F'='lightpink1',
                        '03_BW_M_WG'='dodgerblue',
                        '04_BW_F_WG'='violet',
                        '05_BW_M_Ratier'='dodgerblue4',
                        '06_BW_F_Ratier'='violetred4',
                        # '07_BDW_M_Ratier'='blue',
                        # '08_BDW_F_Ratier'='red',
                        '09_BW_M_EFSA'='cadetblue1',
                        '10_BW_F_EFSA'='lightpink1'),
                      labels=c('01_BW_M'='Body weight EFSA males',
                               '02_BW_F'='Body weight EFSA females',
                               '03_BW_M_WG'='Body weight males WG parameters equation',
                               '04_BW_F_WG'='Body weight females WG parameters equation',
                               '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                               '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                               # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                               # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                               '09_BW_M_EFSA'='Body weight data EFSA males',
                               '10_BW_F_EFSA'='Body weight data EFSA females'),
                      guide="none") +
  
  scale_linetype_manual(name='Predicted',
                        values=c(#'01_BW_M'='solid',
                          #'02_BW_F'='solid',
                          '03_BW_M_WG'='solid',
                          '04_BW_F_WG'='solid',
                          '05_BW_M_Ratier'='solid',
                          '06_BW_F_Ratier'='solid',
                          # '07_BDW_M_Ratier'='solid',
                          # '08_BDW_F_Ratier'='solid',
                          '09_BW_M_EFSA'=NA,
                          '10_BW_F_EFSA'=NA),
                        labels=c(#'01_BW_M'='Body weight EFSA males',
                          #'02_BW_F'='Body weight EFSA females',
                          '03_BW_M_WG'='Body weight males WG parameters equation',
                          '04_BW_F_WG'='Body weight females WG parameters equation',
                          '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                          '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                          # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                          # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                          '09_BW_M_EFSA'='Body weight data EFSA males',
                          '10_BW_F_EFSA'='Body weight data EFSA females')) +
  guides(linetype = guide_legend(override.aes = list(color = c("dodgerblue","violet","dodgerblue4","violetred4")))) + #"cadetblue1","lightpink1", ... ,"blue","red")))) +
  
  scale_shape_manual(name='Observed',
                     values=c(#'01_BW_M'=NA,
                       #'02_BW_F'=NA,
                       '03_BW_M_WG'=NA,
                       '04_BW_F_WG'=NA,
                       '05_BW_M_Ratier'=NA,
                       '06_BW_F_Ratier'=NA,
                       # '07_BDW_M_Ratier'=NA,
                       # '08_BDW_F_Ratier'=NA,
                       '09_BW_M_EFSA'=15,
                       '10_BW_F_EFSA'=16),
                     labels=c(#'01_BW_M'='Body weight EFSA males',
                       #'02_BW_F'='Body weight EFSA females',
                       '03_BW_M_WG'='Body weight males WG parameters equation',
                       '04_BW_F_WG'='Body weight females WG parameters equation',
                       '05_BW_M_Ratier'='Ratier 2024 BodyWeight males',
                       '06_BW_F_Ratier'='Ratier 2024 BodyWeight females',
                       # '07_BDW_M_Ratier'='Ratier 2024 BDW males',
                       # '08_BDW_F_Ratier'='Ratier 2024 BDW females',
                       '09_BW_M_EFSA'='Body weight data EFSA males',
                       '10_BW_F_EFSA'='Body weight data EFSA females')) +
  guides(shape = guide_legend(override.aes = list(color = c("cadetblue1","lightpink1")))) + 
  
  theme_bw() +
  labs(title="Body weight changes over time",x="\nAge (y)", y="Body weight (kg)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_BW_zoom

Variables_df_backup <- Variables_df

#### Organ volumes ----
# Using Ratier et al. (2024) model
# Fraction of arterial blood, calculated from Filser 2000 p.43
Fr_art_blood = 0.0178 / (0.0178 + 0.0533);

# Hematocrit ### IN MAN ### From Supp mat of Brochot et al. 2019
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

Variables_df <- Variables_df_backup

Variables_df <- Variables_df %>%
  select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
  rename(BW_M = BW_M_Ratier_2024) %>%
  rename(BW_F = BW_F_Ratier_2024) %>%
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
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
  
  #here all QtissueFraction, QtotalFraction and CardiacOutput are changing over time?
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

#### Gender-specific data frames ----  
Variables_M <- Variables_df %>%
  select(TIME, age, BW_M, Hct_M, CardOut_M, Vblood_M, Vart_M, Vven_M, Vliver_M, Vstomach_M, Vgut_M, Vkidney_M, Vurinarytract_M,
         Vskin_M, Vadipose_M, Vadrenal_M, Vbone_M, Vbonenonperfused_M, Vbrain_M, Vbreast_M, Vheart_M, Vmarrow_M, Vmuscle_M,
         Vrepro_M, Vpancreas_M, Vspleen_M, Vthyroid_M, Vlung_M,
         Qliver_M, Qstomach_M, Qgut_M, Qkidney_M, Qurinarytract_M,
         Qskin_M, Qadipose_M, Qadrenal_M, Qbone_M, Qbrain_M, Qbreast_M, Qheart_M, Qmarrow_M, Qmuscle_M,
         Qrepro_M, Qpancreas_M, Qspleen_M, Qthyroid_M, Qlung_M) %>%
  mutate(Vtotal_M = Vblood_M + Vliver_M + Vstomach_M + Vgut_M + Vkidney_M + Vurinarytract_M +
           Vskin_M + Vadipose_M + Vadrenal_M + Vbone_M + Vbonenonperfused_M + Vbrain_M + Vbreast_M + Vheart_M + Vmarrow_M + Vmuscle_M +
           Vrepro_M + Vpancreas_M + Vspleen_M + Vthyroid_M + Vlung_M) %>%
  mutate(Fraction_M = Vtotal_M/BW_M) %>%
  mutate(Qtotal_M = Qliver_M + Qstomach_M + Qgut_M + Qkidney_M + Qurinarytract_M +
           Qskin_M + Qadipose_M + Qadrenal_M + Qbone_M + Qbrain_M + Qbreast_M + Qheart_M + Qmarrow_M + Qmuscle_M +
           Qrepro_M + Qpancreas_M + Qspleen_M + Qthyroid_M + Qlung_M)

Variables_fractions_M <- Variables_df %>%
  select(TIME, age, VadrenalFraction_M, VboneFraction_M, VbonenonperfusedFraction_M, VbrainFraction_M, VbreastFraction_M, 
         VheartFraction_M, VmarrowFraction_M, VmuscleFraction_M, VreproFraction_M, VpancreasFraction_M,
         VskinFraction_M, VspleenFraction_M, VthyroidFraction_M, VurinarytractFraction_M, VkidneyFraction_M,
         VlungFraction_M, VgutFraction_M, VstomachFraction_M, VliverFraction_M, VbloodFraction_M, VadiposeFraction_M,
         QadrenalFraction_M, QboneFraction_M, QbrainFraction_M, QbreastFraction_M, 
         QheartFraction_M, QmarrowFraction_M, QmuscleFraction_M, QreproFraction_M, QpancreasFraction_M,
         QskinFraction_M, QspleenFraction_M, QthyroidFraction_M, QurinarytractFraction_M, QkidneyFraction_M,
         QlungFraction_M, QgutFraction_M, QstomachFraction_M, QliverFraction_M, QadiposeFraction_M) %>%
  mutate(VtotalFraction_M = VadrenalFraction_M + VboneFraction_M + VbonenonperfusedFraction_M + VbrainFraction_M + VbreastFraction_M + 
           VheartFraction_M + VmarrowFraction_M + VmuscleFraction_M + VreproFraction_M + VpancreasFraction_M +
           VskinFraction_M + VspleenFraction_M + VthyroidFraction_M + VurinarytractFraction_M + VkidneyFraction_M +
           VlungFraction_M + VgutFraction_M + VstomachFraction_M + VliverFraction_M + VbloodFraction_M + VadiposeFraction_M) %>%
  mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M + 
           QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
           QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
           QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M) 

Variables_F <- Variables_df %>%
  select(TIME, age, BW_F, Hct_F, CardOut_F, Vblood_F, Vart_F, Vven_F, Vliver_F, Vstomach_F, Vgut_F, Vkidney_F, Vurinarytract_F,
         Vskin_F, Vadipose_F, Vadrenal_F, Vbone_F, Vbonenonperfused_F, Vbrain_F, Vbreast_F, Vheart_F, Vmarrow_F, Vmuscle_F,
         Vrepro_F, Vpancreas_F, Vspleen_F, Vthyroid_F, Vlung_F,
         Qliver_F, Qstomach_F, Qgut_F, Qkidney_F, Qurinarytract_F,
         Qskin_F, Qadipose_F, Qadrenal_F, Qbone_F, Qbrain_F, Qbreast_F, Qheart_F, Qmarrow_F, Qmuscle_F,
         Qrepro_F, Qpancreas_F, Qspleen_F, Qthyroid_F, Qlung_F) %>%
  mutate(Vtotal_F = Vblood_F + Vliver_F + Vstomach_F + Vgut_F + Vkidney_F + Vurinarytract_F +
           Vskin_F + Vadipose_F + Vadrenal_F + Vbone_F + Vbonenonperfused_F + Vbrain_F + Vbreast_F + Vheart_F + Vmarrow_F + Vmuscle_F +
           Vrepro_F + Vpancreas_F + Vspleen_F + Vthyroid_F + Vlung_F) %>%
  mutate(Fraction_F = Vtotal_F/BW_F) %>%
  mutate(Qtotal_F = Qliver_F + Qstomach_F + Qgut_F + Qkidney_F + Qurinarytract_F +
           Qskin_F + Qadipose_F + Qadrenal_F + Qbone_F + Qbrain_F + Qbreast_F + Qheart_F + Qmarrow_F + Qmuscle_F +
           Qrepro_F + Qpancreas_F + Qspleen_F + Qthyroid_F + Qlung_F) #so mass balance should be 1

Variables_fractions_F <- Variables_df %>%
  select(TIME, age, VadrenalFraction_F, VboneFraction_F, VbonenonperfusedFraction_F, VbrainFraction_F, VbreastFraction_F, 
         VheartFraction_F, VmarrowFraction_F, VmuscleFraction_F, VreproFraction_F, VpancreasFraction_F,
         VskinFraction_F, VspleenFraction_F, VthyroidFraction_F, VurinarytractFraction_F, VkidneyFraction_F,
         VlungFraction_F, VgutFraction_F, VstomachFraction_F, VliverFraction_F, VbloodFraction_F, VadiposeFraction_F,
         QadrenalFraction_F, QboneFraction_F, QbrainFraction_F, QbreastFraction_F, 
         QheartFraction_F, QmarrowFraction_F, QmuscleFraction_F, QreproFraction_F, QpancreasFraction_F,
         QskinFraction_F, QspleenFraction_F, QthyroidFraction_F, QurinarytractFraction_F, QkidneyFraction_F,
         QlungFraction_F, QgutFraction_F, QstomachFraction_F, QliverFraction_F, QadiposeFraction_F) %>%
  mutate(VtotalFraction_F = VadrenalFraction_F + VboneFraction_F + VbonenonperfusedFraction_F + VbrainFraction_F + VbreastFraction_F + 
           VheartFraction_F + VmarrowFraction_F + VmuscleFraction_F + VreproFraction_F + VpancreasFraction_F +
           VskinFraction_F + VspleenFraction_F + VthyroidFraction_F + VurinarytractFraction_F + VkidneyFraction_F +
           VlungFraction_F + VgutFraction_F + VstomachFraction_F + VliverFraction_F + VbloodFraction_F + VadiposeFraction_F) %>%
  mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F + 
           QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
           QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
           QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F)

#### Figure male organ volumes ----
Figure_M <- ggplot() +
  geom_line(data = Variables_M, aes(x=age, y=BW_M, color="01_BW")) +
  geom_line(data = Variables_M, aes(x=age, y=Vblood_M, color="02_Vblood")) +
  geom_line(data = Variables_M, aes(x=age, y=Vliver_M, color="03_Vliver")) +
  geom_line(data = Variables_M, aes(x=age, y=Vstomach_M, color="04_Vstomach")) +
  geom_line(data = Variables_M, aes(x=age, y=Vgut_M, color="05_Vgut")) +
  geom_line(data = Variables_M, aes(x=age, y=Vkidney_M, color="06_Vkidney")) +
  geom_line(data = Variables_M, aes(x=age, y=Vskin_M, color="07_Vskin")) +
  geom_line(data = Variables_M, aes(x=age, y=Vadipose_M, color="08_Vadipose")) +
  geom_line(data = Variables_M, aes(x=age, y=Vbrain_M, color="09_Vbrain")) +
  geom_line(data = Variables_M, aes(x=age, y=Vmuscle_M, color="10_Vmuscle")) +
  geom_line(data = Variables_M, aes(x=age, y=Vtotal_M, color="11_Vtotal")) +
  
  scale_colour_manual(name='',
                      values=c('01_BW'='black',
                               '02_Vblood'='red',
                               '03_Vliver'='darkred',
                               '04_Vstomach'='indianred',
                               '05_Vgut'='lightpink',
                               '06_Vkidney'='indianred4',
                               '07_Vskin'='darksalmon',
                               '08_Vadipose'='khaki',
                               '09_Vbrain'='mistyrose',
                               '10_Vmuscle'='rosybrown1',
                               '11_Vtotal'='grey'),
                      labels=c('01_BW'='Body weight',
                               '02_Vblood'='Blood',
                               '03_Vliver'='Liver',
                               '04_Vstomach'='Stomach',
                               '05_Vgut'='Gut',
                               '06_Vkidney'='Kidney',
                               '07_Vskin'='Skin',
                               '08_Vadipose'='Adipose tissue',
                               '09_Vbrain'='Brain',
                               '10_Vmuscle'='Muscle tissue',
                               '11_Vtotal'='Sum of all tissues')) +
  
  theme_bw() +
  labs(title="Organ weight changes over time - males",x="\nAge (y)", y="Weight (kg)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_M

#### Figure male organ fractions ----
Figure_fractions_M <- ggplot() +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VbloodFraction_M, color="02_Vblood")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VliverFraction_M, color="03_Vliver")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VstomachFraction_M, color="04_Vstomach")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VgutFraction_M, color="05_Vgut")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VkidneyFraction_M, color="06_Vkidney")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VskinFraction_M, color="07_Vskin")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VadiposeFraction_M, color="08_Vadipose")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VbrainFraction_M, color="09_Vbrain")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VmuscleFraction_M, color="10_Vmuscle")) +
  geom_line(data = Variables_fractions_M, aes(x=age, y=VtotalFraction_M, color="11_Vtotal")) +
  geom_line(data = Variables_M, aes(x=age, y=Fraction_M, color="12_Fraction")) +
  
  scale_colour_manual(name='',
                      values=c('02_Vblood'='red',
                               '03_Vliver'='darkred',
                               '04_Vstomach'='indianred',
                               '05_Vgut'='lightpink',
                               '06_Vkidney'='indianred4',
                               '07_Vskin'='darksalmon',
                               '08_Vadipose'='khaki',
                               '09_Vbrain'='mistyrose',
                               '10_Vmuscle'='rosybrown1',
                               '11_Vtotal'='grey',
                               '12_Fraction'='black'),
                      labels=c('02_Vblood'='Blood',
                               '03_Vliver'='Liver',
                               '04_Vstomach'='Stomach',
                               '05_Vgut'='Gut',
                               '06_Vkidney'='Kidney',
                               '07_Vskin'='Skin',
                               '08_Vadipose'='Adipose tissue',
                               '09_Vbrain'='Brain',
                               '10_Vmuscle'='Muscle tissue',
                               '11_Vtotal'='Sum of all tissues',
                               '12_Fraction'='Sum tissues / BW')) +
  
  theme_bw() +
  labs(title="Organ fractional volume changes over time - males",x="\nAge (y)", y="Fraction of body weight\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_fractions_M

#### Figure male blood flows ----
Figure_Q_M <- ggplot() +
  geom_line(data = Variables_M, aes(x=age, y=CardOut_M, color="01_CardOut")) +
  geom_line(data = Variables_M, aes(x=age, y=Qliver_M, color="03_Qliver")) +
  geom_line(data = Variables_M, aes(x=age, y=Qstomach_M, color="04_Qstomach")) +
  geom_line(data = Variables_M, aes(x=age, y=Qgut_M, color="05_Qgut")) +
  geom_line(data = Variables_M, aes(x=age, y=Qkidney_M, color="06_Qkidney")) +
  geom_line(data = Variables_M, aes(x=age, y=Qskin_M, color="07_Qskin")) +
  geom_line(data = Variables_M, aes(x=age, y=Qadipose_M, color="08_Qadipose")) +
  geom_line(data = Variables_M, aes(x=age, y=Qbrain_M, color="09_Qbrain")) +
  geom_line(data = Variables_M, aes(x=age, y=Qmuscle_M, color="10_Qmuscle")) +
  geom_line(data = Variables_M, aes(x=age, y=Qtotal_M, color="11_Qtotal")) +
  
  scale_colour_manual(name='',
                      values=c('01_CardOut'='black',
                               '03_Qliver'='darkred',
                               '04_Qstomach'='indianred',
                               '05_Qgut'='lightpink',
                               '06_Qkidney'='indianred4',
                               '07_Qskin'='darksalmon',
                               '08_Qadipose'='khaki',
                               '09_Qbrain'='mistyrose',
                               '10_Qmuscle'='rosybrown1',
                               '11_Qtotal'='grey'),
                      labels=c('01_CardOut'='Cardiac output',
                               '03_Qliver'='Liver',
                               '04_Qstomach'='Stomach',
                               '05_Qgut'='Gut',
                               '06_Qkidney'='Kidney',
                               '07_Qskin'='Skin',
                               '08_Qadipose'='Adipose tissue',
                               '09_Qbrain'='Brain',
                               '10_Qmuscle'='Muscle tissue',
                               '11_Qtotal'='Sum of all tissues')) +
  
  theme_bw() +
  labs(title="Blood flow changes over time - males",x="\nAge (y)", y="Blood flow (L/min)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_Q_M

#### Figure female organ volumes ----
Figure_F <- ggplot() +
  geom_line(data = Variables_F, aes(x=age, y=BW_F, color="01_BW")) +
  geom_line(data = Variables_F, aes(x=age, y=Vblood_F, color="02_Vblood")) +
  geom_line(data = Variables_F, aes(x=age, y=Vliver_F, color="03_Vliver")) +
  geom_line(data = Variables_F, aes(x=age, y=Vstomach_F, color="04_Vstomach")) +
  geom_line(data = Variables_F, aes(x=age, y=Vgut_F, color="05_Vgut")) +
  geom_line(data = Variables_F, aes(x=age, y=Vkidney_F, color="06_Vkidney")) +
  geom_line(data = Variables_F, aes(x=age, y=Vskin_F, color="07_Vskin")) +
  geom_line(data = Variables_F, aes(x=age, y=Vadipose_F, color="08_Vadipose")) +
  geom_line(data = Variables_F, aes(x=age, y=Vbrain_F, color="09_Vbrain")) +
  geom_line(data = Variables_F, aes(x=age, y=Vmuscle_F, color="10_Vmuscle")) +
  geom_line(data = Variables_F, aes(x=age, y=Vtotal_F, color="11_Vtotal")) +
  
  scale_colour_manual(name='',
                      values=c('01_BW'='black',
                               '02_Vblood'='red',
                               '03_Vliver'='darkred',
                               '04_Vstomach'='indianred',
                               '05_Vgut'='lightpink',
                               '06_Vkidney'='indianred4',
                               '07_Vskin'='darksalmon',
                               '08_Vadipose'='khaki',
                               '09_Vbrain'='mistyrose',
                               '10_Vmuscle'='rosybrown1',
                               '11_Vtotal'='grey'),
                      labels=c('01_BW'='Body weight',
                               '02_Vblood'='Blood',
                               '03_Vliver'='Liver',
                               '04_Vstomach'='Stomach',
                               '05_Vgut'='Gut',
                               '06_Vkidney'='Kidney',
                               '07_Vskin'='Skin',
                               '08_Vadipose'='Adipose tissue',
                               '09_Vbrain'='Brain',
                               '10_Vmuscle'='Muscle tissue',
                               '11_Vtotal'='Sum of all tissues')) +
  
  theme_bw() +
  labs(title="Organ weight changes over time - females",x="\nAge (y)", y="Weight (kg)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_F

#### Figure female organ fractions ----
Figure_fractions_F <- ggplot() +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VbloodFraction_F, color="02_Vblood")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VliverFraction_F, color="03_Vliver")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VstomachFraction_F, color="04_Vstomach")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VgutFraction_F, color="05_Vgut")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VkidneyFraction_F, color="06_Vkidney")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VskinFraction_F, color="07_Vskin")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VadiposeFraction_F, color="08_Vadipose")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VbrainFraction_F, color="09_Vbrain")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VmuscleFraction_F, color="10_Vmuscle")) +
  geom_line(data = Variables_fractions_F, aes(x=age, y=VtotalFraction_F, color="11_Vtotal")) +
  geom_line(data = Variables_F, aes(x=age, y=Fraction_F, color="12_Fraction")) +
  
  scale_colour_manual(name='',
                      values=c('02_Vblood'='red',
                               '03_Vliver'='darkred',
                               '04_Vstomach'='indianred',
                               '05_Vgut'='lightpink',
                               '06_Vkidney'='indianred4',
                               '07_Vskin'='darksalmon',
                               '08_Vadipose'='khaki',
                               '09_Vbrain'='mistyrose',
                               '10_Vmuscle'='rosybrown1',
                               '11_Vtotal'='grey',
                               '12_Fraction'='black'),
                      labels=c('02_Vblood'='Blood',
                               '03_Vliver'='Liver',
                               '04_Vstomach'='Stomach',
                               '05_Vgut'='Gut',
                               '06_Vkidney'='Kidney',
                               '07_Vskin'='Skin',
                               '08_Vadipose'='Adipose tissue',
                               '09_Vbrain'='Brain',
                               '10_Vmuscle'='Muscle tissue',
                               '11_Vtotal'='Sum of all tissues',
                               '12_Fraction'='Sum tissues / BW')) +
  
  theme_bw() +
  labs(title="Organ fractional volume changes over time - females",x="\nAge (y)", y="Fraction of body weight\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_fractions_F

#### Figure female blood flows ----
Figure_Q_F <- ggplot() +
  geom_line(data = Variables_F, aes(x=age, y=CardOut_F, color="01_CardOut")) +
  geom_line(data = Variables_F, aes(x=age, y=Qliver_F, color="03_Qliver")) +
  geom_line(data = Variables_F, aes(x=age, y=Qstomach_F, color="04_Qstomach")) +
  geom_line(data = Variables_F, aes(x=age, y=Qgut_F, color="05_Qgut")) +
  geom_line(data = Variables_F, aes(x=age, y=Qkidney_F, color="06_Qkidney")) +
  geom_line(data = Variables_F, aes(x=age, y=Qskin_F, color="07_Qskin")) +
  geom_line(data = Variables_F, aes(x=age, y=Qadipose_F, color="08_Qadipose")) +
  geom_line(data = Variables_F, aes(x=age, y=Qbrain_F, color="09_Qbrain")) +
  geom_line(data = Variables_F, aes(x=age, y=Qmuscle_F, color="10_Qmuscle")) +
  geom_line(data = Variables_F, aes(x=age, y=Qtotal_F, color="11_Qtotal")) +
  
  scale_colour_manual(name='',
                      values=c('01_CardOut'='black',
                               '03_Qliver'='darkred',
                               '04_Qstomach'='indianred',
                               '05_Qgut'='lightpink',
                               '06_Qkidney'='indianred4',
                               '07_Qskin'='darksalmon',
                               '08_Qadipose'='khaki',
                               '09_Qbrain'='mistyrose',
                               '10_Qmuscle'='rosybrown1',
                               '11_Qtotal'='grey'),
                      labels=c('01_CardOut'='Cardiac output',
                               '03_Qliver'='Liver',
                               '04_Qstomach'='Stomach',
                               '05_Qgut'='Gut',
                               '06_Qkidney'='Kidney',
                               '07_Qskin'='Skin',
                               '08_Qadipose'='Adipose tissue',
                               '09_Qbrain'='Brain',
                               '10_Qmuscle'='Muscle tissue',
                               '11_Qtotal'='Sum of all tissues')) +
  
  theme_bw() +
  labs(title="Blood flow changes over time - females",x="\nAge (y)", y="Blood flow (L/min)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_Q_F

#### Functions for inclusion in model code ----
# Volumes - Male
varBW_M <- approxfun(Variables_M$TIME, Variables_M$BW_M, rule = 2)
varVblood_M <- approxfun(Variables_M$TIME, Variables_M$Vblood_M, rule = 2)
varVliver_M <- approxfun(Variables_M$TIME, Variables_M$Vliver_M, rule = 2)
varVstomach_M <- approxfun(Variables_M$TIME, Variables_M$Vstomach_M, rule = 2)
varVgut_M <- approxfun(Variables_M$TIME, Variables_M$Vgut_M, rule = 2)
varVkidney_M <- approxfun(Variables_M$TIME, Variables_M$Vkidney_M, rule = 2)
varVskin_M <- approxfun(Variables_M$TIME, Variables_M$Vskin_M, rule = 2)
varVadipose_M <- approxfun(Variables_M$TIME, Variables_M$Vadipose_M, rule = 2)
varVbrain_M <- approxfun(Variables_M$TIME, Variables_M$Vbrain_M, rule = 2)
varVmuscle_M <- approxfun(Variables_M$TIME, Variables_M$Vmuscle_M, rule = 2)
varVtotal_M <- approxfun(Variables_M$TIME, Variables_M$Vtotal_M, rule = 2)
varVurinarytract_M <- approxfun(Variables_M$TIME, Variables_M$Vurinarytract_M, rule = 2)
varVart_M <- approxfun(Variables_M$TIME, Variables_M$Vart_M, rule = 2)
varVven_M <- approxfun(Variables_M$TIME, Variables_M$Vven_M, rule = 2)
varVadrenal_M <- approxfun(Variables_M$TIME, Variables_M$Vadrenal_M, rule = 2)
varVbone_M <- approxfun(Variables_M$TIME, Variables_M$Vbone_M, rule = 2)
varVbonenonperfused_M <- approxfun(Variables_M$TIME, Variables_M$Vbonenonperfused_M, rule = 2)
varVbreast_M <- approxfun(Variables_M$TIME, Variables_M$Vbreast_M, rule = 2)
varVheart_M <- approxfun(Variables_M$TIME, Variables_M$Vheart_M, rule = 2)
varVmarrow_M <- approxfun(Variables_M$TIME, Variables_M$Vmarrow_M, rule = 2)
varVrepro_M <- approxfun(Variables_M$TIME, Variables_M$Vrepro_M, rule = 2)
varVpancreas_M <- approxfun(Variables_M$TIME, Variables_M$Vpancreas_M, rule = 2)
varVspleen_M <- approxfun(Variables_M$TIME, Variables_M$Vspleen_M, rule = 2)
varVthyroid_M <- approxfun(Variables_M$TIME, Variables_M$Vthyroid_M, rule = 2)
varVlung_M <- approxfun(Variables_M$TIME, Variables_M$Vlung_M, rule = 2)
# Blood flows - Male
varCardOut_M <- approxfun(Variables_M$TIME, Variables_M$CardOut_M, rule = 2)
varQliver_M <- approxfun(Variables_M$TIME, Variables_M$Qliver_M, rule = 2)
varQstomach_M <- approxfun(Variables_M$TIME, Variables_M$Qstomach_M, rule = 2)
varQgut_M <- approxfun(Variables_M$TIME, Variables_M$Qgut_M, rule = 2)
varQkidney_M <- approxfun(Variables_M$TIME, Variables_M$Qkidney_M, rule = 2)
varQskin_M <- approxfun(Variables_M$TIME, Variables_M$Qskin_M, rule = 2)
varQadipose_M <- approxfun(Variables_M$TIME, Variables_M$Qadipose_M, rule = 2)
varQbrain_M <- approxfun(Variables_M$TIME, Variables_M$Qbrain_M, rule = 2)
varQmuscle_M <- approxfun(Variables_M$TIME, Variables_M$Qmuscle_M, rule = 2)
varQtotal_M <- approxfun(Variables_M$TIME, Variables_M$Qtotal_M, rule = 2)
varQurinarytract_M <- approxfun(Variables_M$TIME, Variables_M$Qurinarytract_M, rule = 2)
varQadrenal_M <- approxfun(Variables_M$TIME, Variables_M$Qadrenal_M, rule = 2)
varQbone_M <- approxfun(Variables_M$TIME, Variables_M$Qbone_M, rule = 2)
varQbreast_M <- approxfun(Variables_M$TIME, Variables_M$Qbreast_M, rule = 2)
varQheart_M <- approxfun(Variables_M$TIME, Variables_M$Qheart_M, rule = 2)
varQmarrow_M <- approxfun(Variables_M$TIME, Variables_M$Qmarrow_M, rule = 2)
varQrepro_M <- approxfun(Variables_M$TIME, Variables_M$Qrepro_M, rule = 2)
varQpancreas_M <- approxfun(Variables_M$TIME, Variables_M$Qpancreas_M, rule = 2)
varQspleen_M <- approxfun(Variables_M$TIME, Variables_M$Qspleen_M, rule = 2)
varQthyroid_M <- approxfun(Variables_M$TIME, Variables_M$Qthyroid_M, rule = 2)
varQlung_M <- approxfun(Variables_M$TIME, Variables_M$Qlung_M, rule = 2)

# Volumes - Female
varBW_F <- approxfun(Variables_F$TIME, Variables_F$BW_F, rule = 2)
varVblood_F <- approxfun(Variables_F$TIME, Variables_F$Vblood_F, rule = 2)
varVliver_F <- approxfun(Variables_F$TIME, Variables_F$Vliver_F, rule = 2)
varVstomach_F <- approxfun(Variables_F$TIME, Variables_F$Vstomach_F, rule = 2)
varVgut_F <- approxfun(Variables_F$TIME, Variables_F$Vgut_F, rule = 2)
varVkidney_F <- approxfun(Variables_F$TIME, Variables_F$Vkidney_F, rule = 2)
varVskin_F <- approxfun(Variables_F$TIME, Variables_F$Vskin_F, rule = 2)
varVadipose_F <- approxfun(Variables_F$TIME, Variables_F$Vadipose_F, rule = 2)
varVbrain_F <- approxfun(Variables_F$TIME, Variables_F$Vbrain_F, rule = 2)
varVmuscle_F <- approxfun(Variables_F$TIME, Variables_F$Vmuscle_F, rule = 2)
varVtotal_F <- approxfun(Variables_F$TIME, Variables_F$Vtotal_F, rule = 2)
varVurinarytract_F <- approxfun(Variables_F$TIME, Variables_F$Vurinarytract_F, rule = 2)
varVart_F <- approxfun(Variables_F$TIME, Variables_F$Vart_F, rule = 2)
varVven_F <- approxfun(Variables_F$TIME, Variables_F$Vven_F, rule = 2)
varVadrenal_F <- approxfun(Variables_F$TIME, Variables_F$Vadrenal_F, rule = 2)
varVbone_F <- approxfun(Variables_F$TIME, Variables_F$Vbone_F, rule = 2)
varVbonenonperfused_F <- approxfun(Variables_F$TIME, Variables_F$Vbonenonperfused_F, rule = 2)
varVbreast_F <- approxfun(Variables_F$TIME, Variables_F$Vbreast_F, rule = 2)
varVheart_F <- approxfun(Variables_F$TIME, Variables_F$Vheart_F, rule = 2)
varVmarrow_F <- approxfun(Variables_F$TIME, Variables_F$Vmarrow_F, rule = 2)
varVrepro_F <- approxfun(Variables_F$TIME, Variables_F$Vrepro_F, rule = 2)
varVpancreas_F <- approxfun(Variables_F$TIME, Variables_F$Vpancreas_F, rule = 2)
varVspleen_F <- approxfun(Variables_F$TIME, Variables_F$Vspleen_F, rule = 2)
varVthyroid_F <- approxfun(Variables_F$TIME, Variables_F$Vthyroid_F, rule = 2)
varVlung_F <- approxfun(Variables_F$TIME, Variables_F$Vlung_F, rule = 2)
# Blood flows - Female
varCardOut_F <- approxfun(Variables_F$TIME, Variables_F$CardOut_F, rule = 2)
varQliver_F <- approxfun(Variables_F$TIME, Variables_F$Qliver_F, rule = 2)
varQstomach_F <- approxfun(Variables_F$TIME, Variables_F$Qstomach_F, rule = 2)
varQgut_F <- approxfun(Variables_F$TIME, Variables_F$Qgut_F, rule = 2)
varQkidney_F <- approxfun(Variables_F$TIME, Variables_F$Qkidney_F, rule = 2)
varQskin_F <- approxfun(Variables_F$TIME, Variables_F$Qskin_F, rule = 2)
varQadipose_F <- approxfun(Variables_F$TIME, Variables_F$Qadipose_F, rule = 2)
varQbrain_F <- approxfun(Variables_F$TIME, Variables_F$Qbrain_F, rule = 2)
varQmuscle_F <- approxfun(Variables_F$TIME, Variables_F$Qmuscle_F, rule = 2)
varQtotal_F <- approxfun(Variables_F$TIME, Variables_F$Qtotal_F, rule = 2)
varQurinarytract_F <- approxfun(Variables_F$TIME, Variables_F$Qurinarytract_F, rule = 2)
varQadrenal_F <- approxfun(Variables_F$TIME, Variables_F$Qadrenal_F, rule = 2)
varQbone_F <- approxfun(Variables_F$TIME, Variables_F$Qbone_F, rule = 2)
varQbreast_F <- approxfun(Variables_F$TIME, Variables_F$Qbreast_F, rule = 2)
varQheart_F <- approxfun(Variables_F$TIME, Variables_F$Qheart_F, rule = 2)
varQmarrow_F <- approxfun(Variables_F$TIME, Variables_F$Qmarrow_F, rule = 2)
varQrepro_F <- approxfun(Variables_F$TIME, Variables_F$Qrepro_F, rule = 2)
varQpancreas_F <- approxfun(Variables_F$TIME, Variables_F$Qpancreas_F, rule = 2)
varQspleen_F <- approxfun(Variables_F$TIME, Variables_F$Qspleen_F, rule = 2)
varQthyroid_F <- approxfun(Variables_F$TIME, Variables_F$Qthyroid_F, rule = 2)
varQlung_F <- approxfun(Variables_F$TIME, Variables_F$Qlung_F, rule = 2)

BDWseventy <- Variables_df %>% filter(BW_M >= 70.3052 & BW_M < 70.31 )

fiftyYO$VskinFraction_M
