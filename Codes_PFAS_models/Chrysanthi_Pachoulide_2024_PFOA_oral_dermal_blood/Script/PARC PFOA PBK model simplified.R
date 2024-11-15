# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Model File
# CP, 10-11-2024
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)

# INPUT ####
# ------------------------------------------------------ #

# Input variables
Final_variables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-15/Final_variables_M.csv") %>% as.data.frame()


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

exposure_stop <- 50*365         # days
sim_stop <- 80*365 # 80*365     # days

TSTART <- 0
TSTOP <- sim_stop               # years in days
DT <- 1                         # days
TIME <- seq(TSTART,TSTOP,by=DT)


Oralconc <- 0.000187 # ug/kg/day [EFSA 2020; page 143]
Dermconc <- 0.000542 # ug/kg/day; mean of #as.numeric(SumExpPFOA_LB_val[i,14])

## Dosing for NOT lifestage simulation ####
# Choose age accordingly
BDW <- Final_variables_M_df %>%
  filter(age == 20) %>% pull(BDW_M) 

AbsPFOA <- Final_variables_M_df %>%
  filter(age == 20) %>% pull(AbsPFOA) 

Oraldose <- Oralconc*BDW # ug/day
Dermaldose <- AbsPFOA*Dermconc*BDW # ug/day

 
# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #

## For NOT lifestage simulations ####

Indep_parms <- Final_variables_M_df %>% select(kAbl, kAap, fup, PL, PF, PK, PSk, PR, PG) %>% distinct()

Age_parms <- Final_variables_M_df %>%
  filter(age == 20) %>% # Select age accordingly
  select(
    V_skin_M,
    V_skinBarrier_M,
    # V_skinPlasma_M,
    # V_skinTissue_M,
    V_kidney_M,
    V_kidneyPlasma_M,
    V_kidneyTissue_M,
    V_filtrate_M,
    V_gut_M,
    V_liver_M,
    V_adipose_M,
    V_rest_M,
    V_plasma_M,
    Q_skin_M,
    Q_kidney_M,
    Q_gut_M,
    Q_liver_M,
    Q_hepatic_M,
    Q_adipose_M,
    Q_rest_M,
    # Q_lungs_M, # Needs to be incorporated when inhalation exposure is included
    QUr_M,
    GFR_M,
    CLdermalabs_M,
    CL_PltPT_M,
    CL_FiltPT_M,
    CL_FiltPT_Prot_M,
    CLbiliary_M,
    CLfecal_M
  ) 

parms <- c(
  kAbl = Indep_parms$kAbl, # affinity constant basolateral; this is about OAT1 and OAT3 which have affinity for the uptake (plasma to cells)
  kAap = Indep_parms$kAap,  # affinity constant apical; this is about OAT4 which has affinity for the re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
  fup = Indep_parms$fup,   # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
  PL = Indep_parms$PL,     # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
  PA = Indep_parms$PF,     # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
  PK = Indep_parms$PK,     # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
  PSk = Indep_parms$PSk,   # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
  PR = Indep_parms$PR,     # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
  PG = Indep_parms$PG,     # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
  VSk = Age_parms$V_skin_M,          
  VSkb = Age_parms$V_skinBarrier_M,
  # VSkP = Age_parms$V_skinPlasma_M, For future use when updating skin again
  # VSkT = Age_parms$V_skinTissue_M, For future use when updating skin again
  VK = Age_parms$V_kidney_M,
  VKP = Age_parms$V_kidneyPlasma_M,
  VKT = Age_parms$V_kidneyPlasma_M,
  VFil = Age_parms$V_filtrate_M,
  VG = Age_parms$V_gut_M,
  VL = Age_parms$V_liver_M,
  VA = Age_parms$V_adipose_M,        
  VR = Age_parms$V_rest_M,  
  VP = Age_parms$V_plasma_M, 
  QSk = Age_parms$Q_skin_M, 
  QK = Age_parms$Q_kidney_M, 
  QG = Age_parms$Q_gut_M, 
  QL = Age_parms$Q_liver_M, # Liver artery
  QH = Age_parms$Q_hepatic_M, # Includes Portal (gut) and Liver arteries 
  QA = Age_parms$Q_adipose_M, 
  QR = Age_parms$Q_rest_M,
  # QP = Age_parms$Q_lungs_M, # Needs to be incorporated when inhalation exposure is included
  QUr = Age_parms$QUr_M, 
  GFR = Age_parms$GFR_M, # Used to be called CFil
  CLdermalabs = Age_parms$CLdermalabs_M,       
  CL_PltPT = Age_parms$CL_PltPT_M,          
  CL_FiltPT = Age_parms$CL_FiltPT_M, #CL_FiltPT_Prot_M
  CLbiliary = Age_parms$CLbiliary_M, 
  CLfecal = Age_parms$CLfecal_M 
  )


# PBK MODEL ####
# ---------------------------------------------------------------------------- #

PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    ## Dose ----
    Oraldose <- if_else(t <= exposure_stop, Oralconc * BDW, 0)
    Dermaldose <- if_else(t <= exposure_stop, AbsPFOA * Dermconc * BDW, 0)
    
    ## Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # Concentration in plasma (ug/L)
    CG <- AG/VG  # Concentration in gut (ug/L)
    CL <- AL/VL  # Concentration in liver (ug/L)
    CA <- AA/VA  # Concentration in adipose (ug/L)
    CR <- AR/VR  # Concentration in rest the body (ug/L)
    
    CSk <- ASk/VSk  # Concentration in skin (ug/L)
    
    CK <- AK/VK  # Concentration in kidney (ug/L)
    CFil <- AFil/VFil # Concentration in filtrate compartment
    
    CVG <- CG/PG # Concentration leaving gut (ug/L)
    CVL <- CL/PL  # Concentration leaving liver (ug/L)
    CVA <- CA/PA # Concentration leaving adipose (ug/L)
    CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
    CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
    
    CVK <- CK/PK # Concentration leaving the kidney (ug/L)
    
    
    ## Differential equations ----
    
    # Adipose compartment
    dAA <- QA*(CP-CVA) # (ug/h)
    
    # Rest compartment
    dAR <- QR*(CP-CVR) # (ug/h)
    
    # Gut compartment: 
    dAG <- Oraldose + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #
    
    # Excretion fecal: cumulative
    dAEx_feces <- CLfecal*CG 
    
    # Liver compartment
    dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fup 

    # Kidney compartment
    dAK <- QK*(CP -CVK) - GFR*fup*CK + CL_FiltPT*CFil
    
    # Filtrate compartment
    dAFil <- GFR*fup*CK - QUr*CFil - CL_FiltPT*CFil 
    
    # Urine compartment
    dAUr <- QUr*CFil
    
    # Skin compartment
    dASk <- QSk*(CP-CVSk) + Dermaldose 
    
    # Plasma compartment
    dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
      QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR 
    
    
    # Mass Balance
    Atot <- AP + ASk + AG + AL + AA + AK + AFil + AR + AEx_feces + AUr
    dInput <- Oraldose + Dermaldose
    MB = Input - Atot
    
    # End
    
    list(c(dAP, 
           dASk, 
           dAG, 
           dAL, 
           dAA, 
           dAK,
           dAFil,
           dAUr, 
           dAR, 
           dAEx_feces, 
           dInput), 
         c(CP = CP, 
           CSk = CSk, 
           CG = CG, CVG = CVG, 
           CL = CL, CVL = CVL, 
           CA = CA, CVA = CVA,
           CK = CK, 
           CFil = CFil, CVK = CVK, 
           CR = CR, CVR = CVR,
           Atot = Atot,MB = MB)
         )
  })
}

## Initials ####

A_init <- c(AP = 0, 
            # ASkb = 0, 
            ASk = 0, 
            AG = 0, 
            AL = 0, 
            AA = 0, 
            AK = 0,
            AFil = 0,
            AUr = 0, 
            AR = 0, 
            AEx_feces = 0, 
            Input = 0)


## Solving the model ####
output_PFOA <- lsoda(y = A_init, 
                     times = TIME, 
                     func = PBPKmodPFOA_M, 
                     parms = parms)
output.PFOA.df <- as.data.frame(output_PFOA) 
output.df <- Final_variables_M_df %>%
  rename(time = TIME) %>%
  left_join(output.PFOA.df) %>% 
  mutate(time = time/365) %>% 
  rename(Years = time) %>% 
  filter(Years <= 80) #according to sim_stop

# RESULTS ####
# ---------------------------------------------------------------------------- #

# PFOA in plasma of one individual
Plot_PFOA_Plasma <- ggplot()+
  geom_path(data = output.df, aes(x = Years, y = CP))+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")

Plot_PFOA_Plasma
ggsave("PlasmaConcentration.png", dpi = 300)


# PFOA in Kidney of one individual
Plot_PFOA_Kidney <- ggplot()+
  geom_path(data = output.df, aes(x = Years, y = CK, color="Kidney"))+
  # geom_path(data = output.df, aes(x = Years, y = CKT, color="Kidney tissue"))+
  geom_path(data = output.df, aes(x = Years, y = CFil, color="Filtrate"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8)) +
  ylab("Kidney Compartments (ng/ml)") +
  scale_color_manual(values = c("Kidney" = "lightblue",
                                # "Kidney tissue" = "blue",
                                "Filtrate" = "darkblue"))

Plot_PFOA_Kidney
ggsave("KidneyConcentrations.png", dpi = 300)


# PFOA in all organs of one individual
Plot_PFOA_All <- output.df %>% 
  # mutate(CK = CKP + CKT) %>% 
  select(Years, CK, CSk, CL, CG, CA, CR, CP) %>% 
  rename(Kidney = CK, Skin = CSk, Liver = CL, Gut = CG, Adipose = CA, Rest = CR, Plasma = CP) %>% 
  pivot_longer(names_to = "Organ", values_to = "Concentration", Kidney:Plasma) %>% 
  ggplot()+
  geom_path(aes(x = Years, y = Concentration, color = Organ)) +
  facet_wrap(~ Organ)+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 8)) +
  ylab("Organ Concentration (ng/ml)")

Plot_PFOA_All
ggsave("OrganConcentrations.png", dpi = 300)




# # Plot MB
# Plot_PFOA_MB <- ggplot()+
#   geom_path(data = output.df, aes(x = Years, y = MB))+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   ylab("MB")
# Plot_PFOA_MB


# Plot MB
# Plot_PFOA_InputOutput <- ggplot()+
#   # geom_path(data = output.df, aes(x = Years, y = Input, color = "blue"))+
#   geom_path(data = output.df, aes(x = Years, y = Atot, color = "green"))+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   ylab("InputOutput")
# Plot_PFOA_InputOutput

## Plotting per simulation output.df ####
# Plot_PFOA_doses <- ggplot()+
#   geom_line(data = Variables_df, aes(x = output.df, y = Oraldose_M),col="red")+
#   geom_line(data = Variables_df, aes(x = output.df, y = Dermaldose_M),col="blue")+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   #scale_colour_hue()+
#   ylab("Dose (ug)")
# 
# Plot_PFOA_doses

# ## Plotting per age ####
# # Plot_PFOA_doses <- ggplot()+
# #   geom_line(data = Variables_df, aes(x = age, y = Oraldose_M),col="red")+
# #   geom_line(data = Variables_df, aes(x = age, y = Dermaldose_M),col="blue")+
# #   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
# #   #scale_colour_hue()+
# #   ylab("Dose (ug)")
# # 
# # Plot_PFOA_doses
# 
# 
# # PFOA in plasma of one individual
# Plot_PFOA_Plasma <- ggplot()+
#   geom_path(data = output.df, aes(x = age, y = CP))+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   ylab("Plasma (ng/ml)")
# 
# Plot_PFOA_Plasma
# ggsave("PlasmaConcentration.1.png", dpi = 300)
# 
# 
# # PFOA in Kidney of one individual
# Plot_PFOA_Kidney <- ggplot()+
#   geom_path(data = output.df, aes(x = age, y = CKP, color="Kidney blood"))+
#   geom_path(data = output.df, aes(x = age, y = CKT, color="Kidney tissue"))+
#   geom_path(data = output.df, aes(x = age, y = CFil, color="Filtrate"))+
#   theme(axis.text.x = element_text(size = 7),
#         axis.text.y = element_text(size = 7), 
#         axis.title = element_text(size = 8)) +
#   ylab("Kidney Compartments (ng/ml)") +
#   scale_color_manual(values = c("Kidney blood" = "lightblue", 
#                                 "Kidney tissue" = "blue", 
#                                 "Filtrate" = "darkblue"))
# 
# Plot_PFOA_Kidney
# ggsave("KidneyConcentration.2.png", dpi = 300)



