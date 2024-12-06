# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Global Sensitivity Analysis
# CP, 30-12-2024
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
library(data.table)
library(openxlsx)
library(writexl)
library(sensitivity)
library(pksensi)
library(PKNCA)



# INPUT ####
# ------------------------------------------------------ #

# Input variables
Final_variables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-27/Final_variables_M.csv") %>% as.data.frame()


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

# exposure_stop <- 50  #50*365         # days
# sim_stop <- 365 #0.5* 365  #80*365 # 80*365     # days

## Commment Chrysa 05-12-2024: I'm changing these for the sensitivity analysis. Does this affect the SA result?
TSTART <- 0.01
TSTOP <- 24.01         # years in days
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
                V_kidney_M,
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
                QUr_M,
                GFR_M,
                CL_FiltPT_M,
                CLbiliary_M,
                CLfecal_M
        ) 

parms <- c(
        fup = Indep_parms$fup,   # Unbound fraction in plasma Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
        PL = Indep_parms$PL,     # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
        PA = Indep_parms$PF,     # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
        PK = Indep_parms$PK,     # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
        PSk = Indep_parms$PSk,   # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
        PR = Indep_parms$PR,     # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
        PG = Indep_parms$PG,     # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
        VSk = Age_parms$V_skin_M,          
        VK = Age_parms$V_kidney_M,
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
        QA = Age_parms$Q_adipose_M,
        QR = Age_parms$Q_rest_M,
        QUr = Age_parms$QUr_M,
        GFR = Age_parms$GFR_M, # Used to be called CFil
        CL_FiltPT = Age_parms$CL_FiltPT_M,
        CLbiliary = Age_parms$CLbiliary_M,
        CLfecal = Age_parms$CLfecal_M
)


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
                dAK <- QK*(CP -CVK) - GFR*fup*CK + CL_FiltPT*fup*CFil
                # Filtrate compartment
                dAFil <- GFR*fup*CK - QUr*CFil - CL_FiltPT*fup*CFil


                # Urine compartment
                dAUr <- QUr*CFil
                
                # Skin compartment
                dASk <- QSk*(CP-CVSk) + Dermaldose 
                
                # Plasma compartment
                dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
                        QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR 
                
                
                # # Mass Balance
                # Atot <- AP + ASk + AG + AL + AA + AK + AFil + AR + AEx_feces + AUr #AKP + AKT
                # dInput <- Oraldose + Dermaldose
                # MB = Input - Atot
                
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
                       dAEx_feces #,
                       # dInput
                       ), 
                     c(CP = CP, 
                       CSk = CSk, 
                       CG = CG, CVG = CVG,
                       CL = CL, CVL = CVL,
                       CA = CA, CVA = CVA,
                       CK = CK,
                       CFil = CFil,
                       CVK = CVK,
                       CR = CR, CVR = CVR #,
                       # Atot = Atot,MB = MB
                       )
                     
                )
        })
}

outputs = c("CP") # Variable in test uncertainty/sensitivity

## Initials ####

A_init <- c(AP = 0, 
            ASk = Dermaldose, 
            AG = Oraldose,
            AL = 0,
            AA = 0,
            AK = 0,
            AFil = 0,
            AUr = 0,
            AR = 0,
            AEx_feces = 0 #,
            # Input = 0
            )


## Solving the model ####
output_PFOA <- ode(y = A_init,
                   times = TIME,
                   func = PBPKmodPFOA_M, 
                   parms = parms,
                   method="lsoda"
                   )
output.PFOA.df <- as.data.frame(output_PFOA) 


## Sensitivity analysis using pksensi ####


## Define the distribution of the parameters that you will analyse in the sensitivity test 
q <- c( "qunif", #fup = Indep_parms$fup,   # Unbound fraction in plasma Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
        "qunif", #PL = Indep_parms$PL,     # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #PA = Indep_parms$PF,     # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #PK = Indep_parms$PK,     # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #PSk = Indep_parms$PSk,   # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #PR = Indep_parms$PR,     # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #PG = Indep_parms$PG,     # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
        "qunif", #VSk = Age_parms$V_skin_M,
        "qunif", #VK = Age_parms$V_kidney_M,
        "qunif", #VFil = Age_parms$V_filtrate_M,
        "qunif", #VG = Age_parms$V_gut_M,
        "qunif", #VL = Age_parms$V_liver_M,
        "qunif", #VA = Age_parms$V_adipose_M,
        "qunif", #VR = Age_parms$V_rest_M,
        "qunif", #VP = Age_parms$V_plasma_M,
        "qunif", #QSk = Age_parms$Q_skin_M,
        "qunif", #QK = Age_parms$Q_kidney_M,
        "qunif", #QG = Age_parms$Q_gut_M,
        "qunif", #QL = Age_parms$Q_liver_M, # Liver artery
        "qunif", #QA = Age_parms$Q_adipose_M,
        "qunif", #QR = Age_parms$Q_rest_M,
        "qunif", #QUr = Age_parms$QUr_M,
        "qunif", #GFR = Age_parms$GFR_M, # Used to be called CFil
        "qunif", #CL_FiltPT = Age_parms$CL_FiltPT_M,
        "qunif", #CLbiliary = Age_parms$CLbiliary_M,
        "qunif" #CLfecal = Age_parms$CLfecal_M
        ) 


## Set parameter distribution ##
# we use 10% change in all parameters

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit


q.arg <- list(list(min = parms["fup"]*LL, max= parms["fup"]*UL),
              list(min = parms["PL"]*LL, max = parms["PL"]*UL),
              list(min = parms["PA"]*LL, max = parms["PA"]*UL),
              list(min = parms["PK"]*LL, max = parms["PK"]*UL),
              list(min = parms["PSk"]*LL, max = parms["PSk"]*UL),
              list(min = parms["PR"]*LL, max = parms["PR"]*UL),
              list(min = parms["PG"]*LL, max = parms["PG"]*UL),
              list(min = parms["VSk"]*LL, max = parms["VSk"]*UL),
              list(min = parms["VK"]*LL, max = parms["VK"]*UL),
              list(min = parms["VFil"]*LL, max = parms["VFil"]*UL),
              list(min = parms["VG"]*LL, max = parms["VG"]*UL),
              list(min = parms["VL"]*LL, max = parms["VL"]*UL),
              list(min = parms["VA"]*LL, max = parms["VA"]*UL),
              list(min = parms["VR"]*LL, max = parms["VR"]*UL),
              list(min = parms["VP"]*LL, max = parms["VP"]*UL),
              list(min = parms["QSk"]*LL, max = parms["QSk"]*UL),
              list(min = parms["QK"]*LL, max = parms["QK"]*UL),
              list(min = parms["QG"]*LL, max = parms["QG"]*UL),
              list(min = parms["QL"]*LL, max = parms["QL"]*UL),
              list(min = parms["QA"]*LL, max = parms["QA"]*UL),
              list(min = parms["QR"]*LL, max = parms["QR"]*UL),
              list(min = parms["QUr"]*LL, max = parms["QUr"]*UL),
              list(min = parms["GFR"]*LL, max = parms["GFR"]*UL),
              list(min = parms["CL_FiltPT"]*LL, max = parms["CL_FiltPT"]*UL),
              list(min = parms["CLbiliary"]*LL, max = parms["CLbiliary"]*UL),
              list(min = parms["CLfecal"]*LL, max = parms["CLfecal"]*UL)
)



## Create parameter matrix ##  
set.seed(1234)

params <- c(
        "fup",
        "PL",
        "PA",
        "PK",
        "PSk",
        "PR",
        "PG",
        "VSk",
        "VK",
        "VFil",
        "VG",
        "VL",
        "VA",
        "VR",
        "VP",
        "QSk",
        "QK",
        "QG",
        "QL",
        "QA",
        "QR",
        "QUr",
        "GFR",
        "CL_FiltPT",
        "CLbiliary",
        "CLfecal"
        )

length(params) == length(q)

# Create the sequences for each parameter by eFAST
x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, replicate = 1)

dim(x$a) # the array of c(model evaluation, replication, parameters)

## Conduct simulation ##
## Solve ODE through R deSolve package
out <- solve_fun(x, 
                 time = TIME, 
                 func = PBPKmodPFOA_M, 
                 initState = A_init, 
                 outnames = outputs)

saveRDS(object=out, file="out_scaled.rds")
out <- readRDS("out_scaled.rds")

## Output of the Uncertainty analysis ##
pdf("out_scaled.pdf")
pksim(out) # PK plot of the outputs based on the given parameter (Uncertainty analysis)
dev.off()

## Output from the sensitivity analysis ##
pdf("out_scaled.pdf")
plot(out) # plot Time-dependent sensitivity (with 95 % CI)


dev.off()

output_PFOASI <- as.data.frame(print(out["tSI"]))
output_PFOASI$Times <- rownames(output_PFOASI)
output_PFOASI$Times <- as.numeric(as.character(output_PFOASI$Times))
output_PFOASI$Year <- output_PFOASI$Times/(365*24)

write.csv(output_PFOASI, "output_PFOASI.csv", col.names = TRUE)

check(out) # Check sensitivity measurement for parameter fixing
pdf("heat_check_CI_scaled.pdf")
heat_check(out, index = "CI") # heat_check Create heatmap to overview the result of GSA
dev.off()

pdf("heat_check_all.pdf")
heat_check(out, show.all = TRUE)
dev.off()


