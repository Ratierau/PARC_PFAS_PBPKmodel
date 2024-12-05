# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Local Sensitivity Analysis
# CP, 05-12-2024
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

library(ggplot2)
library(deSolve)
library(tidyverse)
library(PKNCA)

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
Final_variables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-27/Final_variables_M.csv") %>% as.data.frame()


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

exposure_stop <- 50*365         # days
sim_stop <- 80*365 # 80*365     # days

TSTART <- 0
TSTOP <- sim_stop               # years in days
DT <- 0.1                         # days
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
    CLbiliary_M,
    CLfecal_M
  ) 

parms <- c(
  # kAbl = Indep_parms$kAbl, # affinity constant basolateral; this is about OAT1 and OAT3 which have affinity for the uptake (plasma to cells)
  # kAap = Indep_parms$kAap, # affinity constant apical; this is about OAT4 which has affinity for the re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
  fup = Indep_parms$fup,   # Unbound fraction in plasma Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
  PL = Indep_parms$PL,     # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
  PA = Indep_parms$PF,     # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
  PK = Indep_parms$PK,     # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
  PSk = Indep_parms$PSk,   # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
  PR = Indep_parms$PR,     # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
  PG = Indep_parms$PG,     # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
  VSk = Age_parms$V_skin_M,          
  # VSkb = Age_parms$V_skinBarrier_M,
  # VSkP = Age_parms$V_skinPlasma_M, For future use when updating skin again
  # VSkT = Age_parms$V_skinTissue_M, For future use when updating skin again
  VK = Age_parms$V_kidney_M,
  # VKP = Age_parms$V_kidneyPlasma_M,
  # VKT = Age_parms$V_kidneyPlasma_M,
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
  # QH = Age_parms$Q_hepatic_M, # Includes Portal (gut) and Liver arteries 
  QA = Age_parms$Q_adipose_M, 
  QR = Age_parms$Q_rest_M,
  # QP = Age_parms$Q_lungs_M, # Needs to be incorporated when inhalation exposure is included
  QUr = Age_parms$QUr_M, 
  GFR = Age_parms$GFR_M, # Used to be called CFil
  # CLdermalabs = Age_parms$CLdermalabs_M,       
  # CL_PltPT = Age_parms$CL_PltPT_M,          
  CL_FiltPT = Age_parms$CL_FiltPT_M, 
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
    dAG <- Oraldose + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #
    
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
    dASk <- QSk*(CP-CVSk) + Dermaldose 
    
    # Plasma compartment
    dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
      QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR 
    
    
    # Mass Balance
    Atot <- AP + ASk + AG + AL + AA + AK + AFil + AR + AEx_feces + AUr #AKP + AKT
    dInput <- Oraldose + Dermaldose
    MB = Input - Atot
    
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
           dAEx_feces, 
           dInput), 
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
            AG = 0, 
            AL = 0, 
            AA = 0, 
            AK = 0,
            # AKP = 0,
            # AKT = 0,
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



## LOCAL SENSITIVITY ANALYSIS ####

pars <- parms %>%  as.list()

SENSI_SOLVE <- function(pars) {
  
  PBPKmodPFOA_M
  
  state   <- A_init
  tout    <- TIME
  
  out <- ode(y = state, parms = pars, times = TIME, func = PBPKmodPFOA_M)
  
  return(as.data.frame(out))
}

SENSI_OUT <- SENSI_SOLVE(pars) #saving the results of the model used for Sensitivity Analysis

#All parameters
SENSI_all <- sensFun(
  func = SENSI_SOLVE,
  parms = pars, 
  tiny = 0.1)
head(SENSI_all)

## Analysing the sensitivity on all the variables ####
## ------------------------------------------------------

summary(SENSI_all)

# value : the nominal value of the parameter
# scale : the scale of the parameter
# L1 : the L1 norm of each parameter's sensitivity function
# ---> it's the average of the absolute values of the sensitivities; represents the sensitivity of the x parameter in the y timepoint
# L2 : the L2 norm of each parameter's sensitivity function
# ---> it's the square root of the average of the squares of the sensitivities; represents the sensitivity of the x parameter at the y timepoint
# Mean : the mean value of the sensitivity function
# Min : the minimum value of the sensitivity function
# Max : the maximum value of the sensitivity function
# N : the number of points used in the sensitivity analysis

# L1 and L2 norms provide a measure of the overall sensitivity of the model output to changes in each parameter

DF_summary_SENSI_all <- summary(SENSI_all) %>% as.data.frame() 
DF_summary_SENSI_all <- DF_summary_SENSI_all %>%  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_all))

DF_summary_SENSI_all %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )
ggsave(here(TodaysOuput, "SENSI_all.png"), dpi = 300)


## Analysing the sensitivity on only one compartment ####

SENSI_CBl <- SENSI_all %>% filter(var == "CBl")
summary(SENSI_CBl)

DF_summary_SENSI_CBl <- summary(SENSI_CBl) %>% as.data.frame() 
DF_summary_SENSI_CBl <- DF_summary_SENSI_CBl %>% 
  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl))

DF_summary_SENSI_CBl %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  ) 
ggsave(here(TodaysOuput, "SENSI_CBl.png"), dpi = 300)



## Analysing the sensitivity on only the Blood concentration of Zearalenone glucuronide ####
## ------------------------------------------------------

SENSI_CBl.glu <- SENSI_all %>% filter(var == "CBl.glu")
summary(SENSI_CBl.glu)

DF_summary_SENSI_CBl.glu <- summary(SENSI_CBl.glu) %>% as.data.frame() %>% 
  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl.glu))

DF_summary_SENSI_CBl.glu %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )
ggsave(here(TodaysOuput, "SENSI_CBl.gluc.png"), dpi = 300)
