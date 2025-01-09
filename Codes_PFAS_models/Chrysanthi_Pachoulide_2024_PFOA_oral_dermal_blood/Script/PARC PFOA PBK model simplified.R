# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Model File
# CP, 10-11-2024
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodelCodes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(ggplot2)
library(deSolve)
library(tidyverse)
library(PKNCA)
library(pracma)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()
theme_CP <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 22, lineheight = unit(0.5, "lines")), # lineheight is adjusting the space between lines
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      axis.line = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      panel.border = element_blank(), 
      panel.background = element_rect((fill = "grey94")),
      panel.grid = element_line(linewidth = 0.2,5, colour = "grey100"), 
      strip.background = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0, "cm"),
      legend.key.width = unit(0.2, "cm"),  # Make legend key width span the whole plot
      legend.key.height = unit(0.2, "cm"),  # Adjust legend key height
    )
}


# INPUT ####
# ------------------------------------------------------ #

# Input variables
# Final_variables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-18/Final_variables_M.csv") %>% as.data.frame()
Final_variables_M_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-18/Final_variables_MVolunteer.csv")



# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

exposure_stop <- 1 #50*365         # days
sim_stop <- 80*365 # 80*365   # days

TSTART <- 0
TSTOP <- sim_stop               # years in days
DT <- 1                         # days
TIME <- seq(TSTART,TSTOP,by=DT)


Oralconc <- 0.000187 # ug/kg/day [EFSA 2020; page 143] 3.96ug
Dermconc <- 0.000542 # ug/kg/day; mean of #as.numeric(SumExpPFOA_LB_val[i,14])


Oraldose <- 3.96 # ug https://doi.org/10.1016/j.envint.2024.109047 # Oralconc*BDW # ug/day


# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #

### Constants ####
#### Physiological  ####
Physio_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/PhysioVariables.csv")

# 50 years old adult, male
Physio_params <- Physio_params %>% 
  filter(age == 50) %>% 
  select(ends_with('_M')) %>% 
  mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935; total blood flow as the sum of the fractional blood flows of all organs on which we have data
  mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96; total volume as the sum of the fractional organ volumes of all organs on which we have data

BW <- Physio_params$BDW_M 
QC <- Physio_params$CardOut_M # This is corrected for hematocrit already so it's plasma


# Fractional organ volumes
VAc <- Physio_params$V_adiposeFraction_M   
MassAc <- Physio_params$AdiposeMass_M      #adipose mass
VGc <- Physio_params$V_gutFraction_M
VKc <- Physio_params$V_kidneyFraction_M
VLc <- Physio_params$V_liverFraction_M
VPc <- Physio_params$V_plasmaFraction_M
Hct <- Physio_params$Hct_M
VSkc <- Physio_params$V_skinFraction_M
VTotc <- Physio_params$VolumesSum 

# Fractional organ blood flows
QAc <- Physio_params$Q_adiposeFraction_M/Physio_params$BloodFlowSum
QGc <- Physio_params$Q_gutFraction_M/Physio_params$BloodFlowSum 
QKc <- Physio_params$Q_kidneyFraction_M/Physio_params$BloodFlowSum 
QLc <- Physio_params$Q_liverFraction_M/Physio_params$BloodFlowSum 
QSkc <- Physio_params$Q_skinFraction_M/Physio_params$BloodFlowSum 


#### Chemical Specific ####
PFOA_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/PFOAParams.csv")

MW <- PFOA_params$MW 

fup <- PFOA_params$fup 

# Partition coefficients
# Calculate Plasma/Rest of the body partition coefficient
# KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
KpRe <- (PFOA_params$KpBr*Physio_params$V_brainFraction_M +
           PFOA_params$KpHe*Physio_params$V_heartFraction_M + 
           PFOA_params$KpLu*Physio_params$V_lungFraction_M +
           PFOA_params$KpMu*Physio_params$V_muscleFraction_M +
           PFOA_params$KpSp*Physio_params$V_spleenFraction_M +
           PFOA_params$KpGo*Physio_params$V_reproFraction_M +
           PFOA_params$KpBo*Physio_params$V_boneFraction_M) / (
             Physio_params$V_brainFraction_M +
               Physio_params$V_heartFraction_M +
               Physio_params$V_lungFraction_M +
               Physio_params$V_muscleFraction_M +
               Physio_params$V_spleenFraction_M +
               Physio_params$V_reproFraction_M +
               Physio_params$V_boneFraction_M)  
PFOA_params$KpRe <- KpRe

# Choose between calculated or initial ones (rat)
# In PFOA_params, PF, PG, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
PAc <- PFOA_params$PF #KpAd #PF # Adipose
PGc <- PFOA_params$PG #KpGu #PG # Gut
PKc <- PFOA_params$PK #KpKi #PK # Kidney
PLc <- PFOA_params$PL #KpLi #PL # Liver
PSkc <- PFOA_params$PSk #KpSk #PSk # Skin
PRc <- PFOA_params$PR #KpRe #PR #Rest


# Renal Clearance
CL_OAT4 <- PFOA_params$CL_OAT4 # L/d/kg protein; in vitro PFOA clearance in HEK239 cells overexpressing OAT4, Louisse et al. 2024
REF_OAT4 <- PFOA_params$REF_OAT4 # 1; as we don't have data on the in vitro expression of OAT4

# Biliary and Fecal Clearances
CLbiliaryc <- PFOA_params$CLbiliaryc #L/d/kg Fujii et al 2015 
CLfecalc <- PFOA_params$CLfaecesc # L/d/kg Fujii et al 2015 

### Final Parameters ####

#### Physiological  ####

# Body volumes (L)
VA <- VAc * BW/0.9 + MassAc/0.9
VG <- VGc * BW
VK <- VKc * BW
VL <- VLc * BW
VP <- VPc * (1-Hct) * BW
VSk <- VSkc * BW
VR <- (VTotc*BW) - (VA + VG + VK + VL + VP + VSk)

VFil <- VK * 0.05 # Kidney filtrate compartment, corresponds to the volume of the collecting system in [ICRP 89 page 149] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089

# Body flows (L/d)
QA <- QAc * QC
QG <- QGc * QC
QK <- QKc * QC
QL <- QLc * QC
QSk <- QSkc * QC
QR <- QC - (QA + QG + QK + QL + QSk)

QUr = 0.022 * BW # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
GFR = 0.18 * QK # L/d, Glomerular filtration rate 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089

#### Compound specific ####

# Partition coefficients
PA <- PAc * fup # Adipose
PG <- PGc * fup # Gut
PK <- PKc * fup # Kidney
PL <- PLc * fup # Liver
PSk <- PSkc * fup # Skin
PR <- PRc * fup #Rest

# Renal clearance
CL_FiltPT <- CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * VK # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0

CLbiliary <- CLbiliaryc * BW #L/d  
CLfecal <- CLfecalc * BW # L/d 


parms <- unlist(c(data.frame(BW, #1
                             QC, #2
                             VA, #3
                             VG, #4
                             VK, #5
                             VL, #6
                             VP, #7
                             VSk, #8
                             VR, #9
                             VFil, #10
                             QA, #11
                             QG, #12 
                             QK, #13
                             QL, #14
                             QSk, #15
                             QR, #16
                             QUr, #17
                             GFR, #18
                             PA, #19
                             PG, #20
                             PK, #21
                             PL, #22
                             PSk, #23
                             PR, #24
                             CL_FiltPT, #25
                             CLbiliary, #26
                             CLfecal, #27
                             fup, #28
                             VAc, #29
                             VGc, #30
                             VKc, #31
                             VLc, #32
                             VPc, #33
                             VSkc, #34
                             QAc, #35
                             QGc, #36
                             QKc, #37
                             QLc, #38
                             QSkc, #39
                             PAc, #40
                             PGc, #41
                             PKc, #42
                             PLc, #43
                             PSkc, #44
                             PRc, #45
                             CL_OAT4, #46 
                             REF_OAT4, #47
                             CLbiliaryc, #48
                             CLfecalc #49
)))

parms


# PBK MODEL ####
# ---------------------------------------------------------------------------- #

PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    ## Dose ----
    # Oraldose <- if_else(t <= exposure_stop, Oralconc * BDW, 0)
    # Dermaldose <- if_else(t <= exposure_stop, AbsPFOA * Dermconc * BDW, 0)
    # 
    ## Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # Concentration in plasma (ug/L)
    
    CG <- AG/VG  # Concentration in gut (ug/L)
    CL <- AL/VL  # Concentration in liver (ug/L)
    CA <- AA/VA  # Concentration in adipose (ug/L)
    CR <- AR/VR  # Concentration in rest the body (ug/L)
    
    CSk <- ASk/VSk  # Concentration in skin (ug/L)
    
    CK <- AK/VK  # Concentration in kidney (ug/L)
    # CKP <- AKP/VKP # Concentration in kidney plasma (ug/L) 
    # CKT <- AKT/VKT # Concentration in kidney tissue (ug/l) 
    CFil <- AFil/VFil # Concentration in filtrate compartment
    
    CVG <- CG/PG # Concentration leaving gut (ug/L)
    CVL <- CL/PL  # Concentration leaving liver (ug/L)
    CVA <- CA/PA # Concentration leaving adipose (ug/L)
    CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
    CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
    
    CVK <- CK/PK # Concentration leaving the kidney (ug/L)
    # CVKP <- CKP/PK # Concentration leaving kidney tissue (ug/l)
    
    
    ## Differential equations ----
    
    # Adipose compartment
    dAA <- QA*(CP-CVA) # (ug/h)
    
    # Rest compartment
    dAR <- QR*(CP-CVR) # (ug/h)
    
    # Gut compartment: 
    dAG <- + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #Oraldose
    
    # Excretion fecal: cumulative
    dAEx_feces <- CLfecal*CG 
    
    # Liver compartment
    dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fup 
    
    # Kidney compartment
    dAK <- QK*(CP -CVK) - GFR*fup*CK + CL_FiltPT*fup*CFil
    # dAKP <- QK*(CP - CVKP) - CL_PltPT*fup*CKP + CL_FiltPT*fup*CKT - CL_PltPT*fup*CKP
    # dAKT <- CL_PltPT*fup*CKP + CL_FiltPT*CFil - CL_FiltPT*fup*CKT
    # Filtrate compartment
    dAFil <- GFR*fup*CK - QUr*CFil - CL_FiltPT*fup*CFil 
    # dAFil <- CL_PltPT*fup*CKP - QUr*CFil - CL_FiltPT*CFil 
    
    
    # Urine compartment
    dAUr <- QUr*CFil
    
    # Skin compartment
    dASk <- QSk*(CP-CVSk) #+ Dermaldose 
    
    # Plasma compartment
    dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
      QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR 
    
    
    # Mass Balance
    Atot <- AP + ASk + AG + AL + AA + AK + AFil + AR + AEx_feces + AUr #AKP + AKT
    # dInput <- Oraldose + Dermaldose
    # MB = Input - Atot
    MB = Oraldose - Atot
    
    # End
    
    list(c(dAP, 
           dASk, 
           dAG, 
           dAL, 
           dAA,
           dAK,
           # dAKT,
           # dAKP,
           dAFil,
           dAUr, 
           dAR, 
           dAEx_feces #,
           # dInput
           ), 
         c(CP = CP, 
           CSk = CSk, 
           CG = CG, CVG = CVG, 
           CL = CL, CVL = CVL, 
           CA = CA, CVA = CVA,
           CK = CK,
           # CKP = CKP, CKT = CKT, 
           CFil = CFil, 
           CVK = CVK,
           # CVKP = CVKP, 
           CR = CR, CVR = CVR,
           Atot = Atot,MB = MB)
         )
  })
}

## Initials ####

A_init <- c(AP = 0, 
            # ASkb = 0, 
            ASk = 0, 
            AG = Oraldose, 
            AL = 0, 
            AA = 0, 
            AK = 0,
            # AKP = 0,
            # AKT = 0,
            AFil = 0,
            AUr = 0, 
            AR = 0, 
            AEx_feces = 0 #,
            #Input = 0
            )


## Solving the model ####
output_PFOA <- lsoda(y = A_init, 
                     times = TIME, 
                     func = PBPKmodPFOA_M, 
                     parms = parms)
output.PFOA.df <- as.data.frame(output_PFOA) 
# output.df <- Final_variables_M_df %>%
#   rename(time = TIME) %>%
#   left_join(output.PFOA.df) %>% 
#   mutate(time = time/365) %>% 
#   rename(Years = time) #%>% 
#   filter(Years <= 80) #according to sim_stop

# # RESULTS ####
# # ---------------------------------------------------------------------------- #
# 
# # PFOA in plasma of one individual
# Plot_PFOA_Plasma <- ggplot()+
#   geom_path(data = output.df, aes(x = Years, y = CP))+
#   theme_CP()+
#   # theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   ylab("Plasma (ng/ml)")
# 
# Plot_PFOA_Plasma
# ggsave("PlasmaConcentration.png", dpi = 300)
# 
# 
# # PFOA in Kidney of one individual
# Plot_PFOA_Kidney <- ggplot()+
#   geom_path(data = output.df, aes(x = Years, y = CK, color="Kidney"))+
#   # geom_path(data = output.df, aes(x = Years, y = CKP, color="Kidney plasma"))+
#   # geom_path(data = output.df, aes(x = Years, y = CKT, color="Kidney tissue"))+
#   geom_path(data = output.df, aes(x = Years, y = CFil, color="Filtrate"))+
#   theme_CP()+
#   # theme(axis.text.x = element_text(size = 7),
#   #       axis.text.y = element_text(size = 7),
#   #       axis.title = element_text(size = 8)) +
#   ylab("Kidney Compartments (ng/ml)") +
#   scale_color_manual(values = c("Kidney" = "lightblue",
#                                 # "Kidney plasma" = "lightblue",
#                                 # "Kidney tissue" = "blue",
#                                 "Filtrate" = "darkblue"))
# 
# Plot_PFOA_Kidney
# ggsave("KidneyConcentrations.png", dpi = 300)
# 
# 
# # PFOA in all organs of one individual
# Plot_PFOA_All <- output.df %>% 
#   # mutate(CK = CKP + CKT) %>% 
#   select(Years, CK, CSk, CL, CG, CA, CR, CP) %>% 
#   rename(Kidney = CK, Skin = CSk, Liver = CL, Gut = CG, Adipose = CA, Rest = CR, Plasma = CP) %>% 
#   pivot_longer(names_to = "Organ", values_to = "Concentration", Kidney:Plasma) %>% 
#   ggplot()+
#   geom_path(aes(x = Years, y = Concentration, color = Organ)) +
#   facet_wrap(~ Organ)+
#   theme_CP()+
#   theme(legend.position = "none")+
#   # theme(axis.text.x = element_text(size = 7),
#   #       axis.text.y = element_text(size = 7), 
#   #       axis.title = element_text(size = 8)) +
#   ylab("Organ Concentration (ng/ml)")
# 
# Plot_PFOA_All
# ggsave("OrganConcentrations.png", dpi = 300)
# 
# 


# RESULTS FOR STUDY OF 1 INDIVIDUAL####
# ---------------------------------------------------------------------------- #

output.PFOA.df <- output.PFOA.df %>% 
  # mutate(time = time/365) %>% 
  rename(Days = time)


# PFOA in plasma of one individual
Plot_PFOA_Plasma <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CP))+
  theme_CP()+
  # theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")

Plot_PFOA_Plasma
ggsave("PlasmaConcentration.png", dpi = 300)


# PFOA in Kidney of one individual
Plot_PFOA_Kidney <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CK, color="Kidney"))+
  # geom_path(data = output.PFOA.df, aes(x = Days, y = CKP, color="Kidney plasma"))+
  # geom_path(data = output.PFOA.df, aes(x = Days, y = CKT, color="Kidney tissue"))+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CFil, color="Filtrate"))+
  theme_CP()+
  # theme(axis.text.x = element_text(size = 7),
  #       axis.text.y = element_text(size = 7),
  #       axis.title = element_text(size = 8)) +
  ylab("Kidney Compartments (ng/ml)") +
  scale_color_manual(values = c("Kidney" = "lightblue",
                                # "Kidney plasma" = "lightblue",
                                # "Kidney tissue" = "blue",
                                "Filtrate" = "darkblue"))

Plot_PFOA_Kidney
ggsave("KidneyConcentrations.png", dpi = 300)


# PFOA in all organs of one individual
Plot_PFOA_All <- output.PFOA.df %>% 
  # mutate(CK = CKP + CKT) %>% 
  select(Days, CK, CSk, CL, CG, CA, CR, CP) %>% 
  rename(Kidney = CK, Skin = CSk, Liver = CL, Gut = CG, Adipose = CA, Rest = CR, Plasma = CP) %>% 
  pivot_longer(names_to = "Organ", values_to = "Concentration", Kidney:Plasma) %>% 
  ggplot()+
  geom_path(aes(x = Days, y = Concentration, color = Organ)) +
  facet_wrap(~ Organ)+
  theme_CP()+
  theme(legend.position = "none")+
  # theme(axis.text.x = element_text(size = 7),
  #       axis.text.y = element_text(size = 7), 
  #       axis.title = element_text(size = 8)) +
  ylab("Organ Concentration (ng/ml)")

Plot_PFOA_All
ggsave("OrganConcentrations.png", dpi = 300)





# Calculate AUC and Half life ####
# ---------------------------------------------------------------------------- #

AUC <- trapz(output_PFOA[ , "time"], output_PFOA[ , "CP"]) #ug*day/L
  
# Calculate predicted half-life
time <- output_PFOA[ , "time"] #in days
conc <- output_PFOA[ , "CP"] #(ug/L) or (ng/ml) 
Cmax <- max(conc)
Tmax <- time[which.max(conc)]
tlast <- max(time[conc > 0])
 
half_life <- pk.calc.half.life(
    conc,
    time,
    Tmax,
    tlast
  )

HalfLife <- half_life$half.life/365 #halflife in years

# Experimental
ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/HalfLifes.csv")

Exp_HalfLifes <- ExpData %>%
  filter(species == "human",
         chemical == "pfoa",
         parameter== "HalfLife") %>%
  select(value_average) %>%
  rename(HalfLife = value_average)
Exp_HalfLifes$HalfLife <- as.numeric(Exp_HalfLifes$HalfLife)

Predicted.df <- data.frame(
  HalfLife = HalfLife,
  Origin = "Predicted",
  value = 1)
Experimental.df <- data.frame(
  HalfLife = Exp_HalfLifes$HalfLife,
  Origin = "Experimental",
  value = 1)

HalfLifes <- rbind(Predicted.df, Experimental.df)

Plot_HalfLifes <- HalfLifes %>%
  ggplot()+
  geom_violin(data = Experimental.df, aes(value, HalfLife),
              color = "transparent",
              fill = "grey")+
  geom_point(data = Predicted.df, aes(value, HalfLife),
             color = "slateblue3", size = 5, shape = 18) +
  ylab("Half life (years)") +
  theme_CP() +
  theme(axis.text.x=element_blank(), #remove x axis labels
      axis.ticks.x=element_blank(), #remove x axis ticks
      axis.title.x = element_blank()
      )
Plot_HalfLifes
ggsave("ExpVsSimHalfLife.png", dpi = 300)


## Plasma concentrations ####
ExpPlasma <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/Experimental.Plasma.PFOA.csv")

ExpPlasma <- ExpPlasma %>% 
  # mutate(Years = Time_days/365) %>% 
  # select(!Time_days) %>% 
  rename(Days = Time_days) %>% 
  rename(CP = Plasma_concentration_ngPerml) %>% 
  filter(Days <=TSTOP)


Plot_Plasma <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CP), color = "aquamarine")+
  geom_point(data = ExpPlasma, aes(x = Days, y = CP), color = "black")+
  theme_CP()+
  # scale_color_manual(values = c("Simulated" = "aquamarine",
  #                               "Experimental" = "black")) +
  # # theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")
Plot_Plasma
ggsave("PlasmaExpVsPredicted.png", dpi = 300)
