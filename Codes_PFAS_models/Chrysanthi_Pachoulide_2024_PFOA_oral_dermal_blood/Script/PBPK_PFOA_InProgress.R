####################################################################################################
# Updating the Husoy model, in the context of the PARC project
# The goal is to have a PFAS generic PBPK model to be used for risk assessment purposes by EFSA. 
# The model should be as simple as possible and as mechanistic as necessary to capture PFAS specific purposes, where in vitro data are available for different PFAS's
# The first step is to update the PFOA model *Husoy_2023 :
#   - the most updated dermal and oral exposure information
#   - a more mechanistic dermal absorption (as of now the model is using the results of a study on 1 human)
#   - simplify the kidney compartment, change Tm and kt which are fitted values with Vmax and Km from the most recent studies
#   - remove the "Free" term which is a fitted value from by historical PFAS models
#   - consider changing the partition coefficients with either Human biopsy concentrations/ or mechanistic coefficients
# The second step is to change the physiological parameters with human life-stage equations to reflect the ages from 0 to 80 years old
#####################################################################################################

# Last modified by: Chrysanthi Pachoulide
# Notes: Starting from initial Trine's model again 
#####################################################################################################

rm(list=ls()) # to clear out the global environment



# set work directory
HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Script"
setwd(HOME)

# creating a new folder to store the results 
WhatAmITesting <- "CHANGE_KIDNEY"
Saveoutput <- file.path('C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Results', Sys.Date(), WhatAmITesting)
dir.create(Saveoutput, WhatAmITesting, recursive = TRUE)

# bring my plotting function
source("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Script/Theme_CP_function.R")

# load packages
library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)



## Read in data ##

# PFOA_LB_dummy <-  read.delim("./Data/PFOA_LB_dummy.csv", sep = ";")

nPeople <- 1 #as.numeric(nrow(PFOA_LB_dummy))

## empty databases to store the results ##

PFOASerum <- matrix(0, ncol = nPeople, nrow = 438001) # nrow = 50 (years)*365 (days)*24(hours) +1
PFOAFat <- matrix(0, ncol = nPeople, nrow = 438001)
PFOAUrine <- matrix(0, ncol = nPeople, nrow = 438001)
PFOAKidney <- matrix(0, ncol = nPeople, nrow = 438001)
PFOALiver <- matrix(0, ncol = nPeople, nrow = 438001)

# PFOASerumTestTot # ??? what is this ???

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### PBPK model PFOA ####
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for (i in 1:nPeople) { 
  
  ## Time parameters ##
  
  starttime = 0 # in hours
  stoptime = 50*365*24  # end of simmulation (h) - (50 years) run it with 0.5 år
  dtout <- 1 # resolution of the output time (h)
  times <- seq(starttime, stoptime, by=dtout) # not in the original code 
  Tolerance <- 0.01 # default tolerance
  dtmax = 10.0
  dtmin = 0.000001
  year = times/(24*365) # A bit uncertain if I need to define the time
  
  # Do I need to define the time like below
  # Times <- seq(seq(starttime, stoptime, dtout)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Physiological parameters (from Brown, et al 1997)
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  QCC = 12.5  # Cardiac blood output (L/h/kg^0.75)
  QFC = 0.052 # Fraction cardiac output going to fat
  QLC = 0.069 # Fraction cardiac output going to liver suggested by EFSA is 0.069, since they corrected for output to the gut
  # Luccisano used 0.25, but Brown actually report 0.23
  QKC = 0.175 # Fraction cardiac output going to kidney
  QSkC = 0.058 # Fraction cardiac output going to skin
  QGC = 0.181 # Fraction of cardiac output going to gut and the liver via portal arthery
  
  BW = 70.308 # #mean of as.numeric(SumExpPFOA_LB_val[i,3])      # Body weight from the EuroMix study
  
  ## fractional tissue columes ##
  
  VLC = 0.026  # Fraction liver volume
  VFC = 0.214  # Fraction fat volume
  VKC = 0.004  # Fraction kidney volume
  VfilC = 0.0004  # Fraction filtrate compartment volume (10% of kidney volume)
  VGC = 0.0171  # Fraction gut volume
  VPlasC = 0.0428  # Fraction plasma volume (58% of blood)
  Htc = 0.44  # hematocrit
  
  Skinarea = 972  # Exposed area on skin (cm^2); Changes from 5 since 2018 EFSA
  SkinTarea = 9.1*(BW*1000)^0.666  # Total area of the skin (cm^2)
  Skinthickness = 0.1 # Skin thickness (cm)
  
  ## UPDATED: Skin barrier or epidermis (stratum corneum & viable epidermis)
  SkbTarea = 9.1*BW*1000^0.666  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  fSkbarea = 0.107 # (cm^2); exposed skin area = surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
  # checked the reference, 0.972m^2 is the 95th percentile of the lower extremities; I think that was a mistake when copying the information as the hand data is just above. For hands the 95th percentile is 0.131 (mean value is 0.107m^2)
  Skbthickness = 83.1/1000 # (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
  VSkb = fSkbarea*Skbthickness/1000  # (L); Skin barrier volume; as previously coded by Trine

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Scaled cardiac output and blood flows##
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  QC = QCC*BW^0.75  # Cardiac output adjusted for BW (L/h)
  QCP = QC*(1-Htc)  # Cardiac output adjusted for plasma flow (L/h)
  QL = QLC*QCP  # Plasma flow to liver (L/h)
  QF = QFC*QCP  # Plasma flow to fat (L/h)
  QK = QKC*QCP  # Plasma flow to kidney (L/h)
  
  QG = QGC*QCP  # Plasma flow to gut (L/h)
  
  QSk <- QSkC*QCP*(Skinarea/SkinTarea) #ifelse(Dermconc>0.0,QSkC*QCP*(Skinarea/SkinTarea),0.0) # plasma flow to the skin
  
  QR = QCP-QL-QF-QK-QG-QSk # Plasma flow to the rest of the body.
  
  
  ## Scaled tissue volumes ##
  
  VL = VLC*BW  # Liver volume (L)
  VF = VFC*BW  # Fat volume (L)
  VK = VKC*BW  # Kidney volume (L)
  Vfil = VfilC*BW  # Filtrate compartment volume
  VG = VGC*BW  # Gut volume (L)
  VPlas = VPlasC*BW  # Plasma volume(L)
  
  VSk = (Skinarea*Skinthickness)/1000  # Skin volume (L)
  VR = 0.84*BW-VL-VF-VK-Vfil-VG-VPlas-VSk  # Rest of the body volume (L). Need to know where the number 0.84 comes from???
  
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Chemical specific parameters (PFOA)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Free = 0.02 #initial free fraction; 0.061 # from Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415; OLD FREE = 0.02 # Free fraction of PFOA in plasma
  # fut <- 1/(1+((1-Free)/Free)*0.5) #fraction unbound in the tissues
  # Free/fut
  
  PL =2.2  # Plasma/liver partition coefficient
  PF = 0.04  # Plasma/fat partition coefficient
  PK = 1.05  # Plasma/kidney partition coefficient
  PSk = 0.1  # Plasma/skin partition coefficient
  PR = 0.12  # Plasma/rest of the body partition coefficient
  PG = 0.05  # Plasma/gut partition coefficient
  

  AbsPFOA =  0.016 #0.00048     # Changed to the absorption measured by Abraham and Monien 2022, of 1.6% of applied dose from sunscreen. 
  #NO LONGR USING Dermal absorption fraction take from Fasano et al 2008. This is in line with the absorption reported from Franko et al 2012 divided on 1000 since PFOA is ionized 
  
  # kidney and urine #
  
  Tmc = 5842.308 # average of as.numeric(SumExpPFOA_LB_val[i,5]) #5000  # ug/h/kg^0.75 Maximum resorption rate, changed from 6 in the original Loccisano 2011 model (ug)
  # representing a half-life of 2.3 years
  
  Tm = Tmc*BW^0.75 # transporter maximum
  
  Kt = 55  # Resorption affinity, changed from 0.055 in the original Loccisano 2011 model (ug)
  # Expressed in ug to be consistent with the other parameters
  
  ## UPDATED TRANSPORTER KINETICS ####
  ### Kidney scaling factors
  MK = VKC*BW*1000	 # total kidney mass (g)
  Kcells = 6E7	# number of cells per gram kidney
  Kprotein = 2.0e-9	# gram protein/proximal tubule cell
  MPT = MK*Kcells*Kprotein	#mass proximal tubule in gram based on BW
  SFOAT4 = 1 # for now it's 1; could be 1 pmole/g tissue from https://doi.org/10.1124/dmd.119.086579, but this is used for scaling between animals; as they mention in the article: "the data obtained from the absolute peptide approach should not be considered as absolute molar protein abundance data because complete trypsin digestion may not be confirmed" 
  
  MW = 414.07 #PFOA MW
  Vmaxc = 4.5*MW/1000*60 #nmol/min/mg protein *MW/1000*60 ug/h ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-6
  Vmax = Vmaxc * MPT * SFOAT4 #ug/h 
  Km = 47*MW #ug/L ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-6
  
  
  CLurinec = 0.00000183  # 0.044/24/1000 urinary clearance L/h/kg calculated from the urinary clearance of 0.044 ml/h/kg  taken from Fujii et al 2015 
  # 0.0003 urinary clearance (/h/kg^-0.25); estimated from Harada et al 2005 NOT USED
  CLurine = CLurinec*BW^(-0.25) # clearance urine (L/h), 
 

  QUr = 0.0625 #L/h (total urine loss per day is 1.5L, assuming that this also reflects the formation of urine)
  kUr = 0.03125 #/h (total volume capacity of bladder is 500ml, assuming a urine production of 1.5L/24h, the bladder needs to be emptied 3 times per day; the rate constant of emptying is the filling rate (=QUr)*max volume (0.5L))
  
  CLfil = 0.1*QK #L/h GFR
  # Changing the CLfil, as from what I've read GFR is 10% of total kidney blood flow
  # CLfil = 0.2*QK  # (used to be kfil) Clearance from the kidney to the filtrate compartment (L/h); 20% of bloodstream to QK is cleared for 
  
  
  # Not used anymore
  # CLurinec = 0.00000183  # (used to be CLurinec) 0.044/24/1000 urinary clearance L/h/kg calculated from the urinary clearance of 0.044 ml/h/kg  taken from Fujii et al 2015 
  # # 0.0003 urinary clearance (/h/kg^-0.25); estimated from Harada et al 2005 NOT USED
  # CLurine = CLurinec*BW^(-0.25) # (used to be CLurine) clearance urine (L/h), 
  
  
  # Clearance parameters# 
  
  CLbiliaryc = 0.000109 # 2.62/24/1000 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day  taken from Fujii et al 2015 
  CLfaecesc = 0.00000217 # 0.052/24/1000 faeces clearance L/h/kg clearance in faeces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg
  
  CLbiliary = CLbiliaryc*BW^0.1# L/h biliary clearence, BW adjusted to the volume from liver
  CLfaeces = CLfaecesc*BW^0.001  # L/h faeces clearance, BW adjusted to the volume of GI tract
  
  #kfil = 0.2*QK  # Clearance from the kidney to the filtrate compartment (L/h); 20% of bloodstream to QK is cleared for 
  Qfil = 0.2&QK
  
  ## Free fraction of chemical in tissues ##
  
  FreeL = Free/PL  # liver
  FreeF = Free/PF  # fat
  FreeK = Free/PK  # kidney
  FreeSk = Free/PSk  # skin
  FreeR = Free/PR  # rest of the body
  FreeG = Free/PG  # gut
  
  ## THIS IS ASSUMING THAT THE PARTITION COEFFICIENTS ARE CALCULATED USING THE FREE FRACTION ONLY; HOWEVER BASED ON THE PT RESOURCE, THE CONCENTRATIONS WERE TOTAL TISSUE CONCENTRATIONS 
  # TO BE REPLACED BY THE DATA FROM THE NEW FISCHER STUDY AND INTEGRATED IN THE TISSUE PARTITION COEFFICIENT: https://pubs.acs.org/doi/10.1021/acs.est.4c04050#
  # PL = Free/FreeP # liver
  # PF = Free/FreeF # fat
  # PK = Free/FreeK  # kidney
  # PSk = Free/FreeSk  # skin
  # PR = Free/FreePR  # rest of the body
  # PG = Free/FreePG  # gut
  
  # # UPDATED Skin absorption
  # Papp = 3.82 * 10^-3 #cm/h ref: https://doi.org/10.1016/j.envint.2024.108772
  # kperm = Papp/Skbthickness  #/h

  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Dosing parameters
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  tchng = 50*365*24  # Duration of exposure (h); 50 years; turn dose on and off
  
  ## Dermal exposure ##
  
  Dermconc = 0.000542 #mean of #as.numeric(SumExpPFOA_LB_val[i,14])
  Dermdose = Dermconc*BW*AbsPFOA     # Internal dose from dermal absorption (Ug/day)
  #Dermdose = Dermconc*BW #(ug/day)
  
  
  ## Oral exposure ##
  
  Oralconc = 0.000103 # mean of as.numeric(SumExpPFOA_LB_val[i,7])  # Oral uptake /ug/kg bw/day
  Oraldose = Oralconc*BW  # (ug/day)
  
  
  
  Tinput = 24  # durration of dose (h)
  
  tinterval = 24 # the interval the dose should be repeated (h). Do not have to be similar to Tinput
  # Not originally in the EFSA model since in Bercley Madonna this is included in an predefined function
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Compile parameters 
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  para <- unlist(c(data.frame(
    Htc,
    Tmc,
    Kt,
    Vmax,
    Km, 
    kUr,
    Free,
    BW,
    CLurinec,
    PL,
    PF,
    PK,
    PSk,
    PR,
    PG,
    CLurine,
    CLbiliaryc,
    CLfaecesc,
    CLbiliary, 
    CLfaeces,
    FreeL,
    FreeF,
    FreeK,
    FreeSk,
    FreeR,
    FreeG,
    tchng,
    QC,
    QCP,
    QL,
    QF,
    QK,
    kfil,
    Qfil,
    QG,
    QSk,
    QR,
    VL,
    VF,
    VK,
    Vfil,
    VG,
    VPlas,
    VSk,
    VR,
    Tm,
    SkinTarea,
    AbsPFOA
    #kperm
  )))
  
  para
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Initial conditions
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  yini <- unlist(c(data.frame(
    Input1 = 0.0,  # amount of PFOA entering the body via oral dosing
    Input2 = 0.0, # amount of PFOA entering the body via dermal exposure
    APlas = 0.0,  # amount of PFOA in Plasma
    AG = 0.0,  # amount of PFOA in gut
    AL = 0.0,  # amount of PFOA in liver
    AF = 0.0,  # amount of PFOA in fat
    AK = 0.0,  # amount of PFOA in kidney
    # Afil = 0.0,  # amount of PFOA in kidney filtrate
    Afil = 0.0,  # amount of PFOA in kidney filtrate
    Adelay = 0.0,  # amount of PFOA in storage compartment of urine
    Aurine = 0.0,  # amount of PFOA in urine
    # AUr = 0.0,  # amount of PFOA in storage compartment of urine
    # AEx_ur = 0.0,  # amount of PFOA in urine
    Afaeces = 0.0, # amount of PFOA in faeces
    ASk = 0.0,  # amount of PFOA in skin'
    AR = 0.0  # amount of PFOA in the rest of the body
  )))

  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Model for PFOA
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  PBPKmodPFOA <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      if(t<tchng){DoseOn=1} else{DoseOn=0}
      
      Input1 <- Oraldose/Tinput*(t %% tinterval<Tinput) 
      Input2 <- Dermdose/Tinput*(t %% tinterval<Tinput)

      
      # concentrations in the organs #
      CAFree <- APlas/VPlas  # free concentration of PFOA in plasma in (ug/L)
      CA <- CAFree/Free  # total concentration in plasma
      #CA <- APlas/VPlas # total concentration in plasma #RECODED
      CG <- AG/VG  # concentration of PFOA in gut (ug/L)
      CVG <- CG/PG # concentration of PFOA leaving gut (ug/L)
      CL <- AL/VL  # concentration of PFOA in liver (ug/L)
      CVL <- CL/PL  # concentration of PFOA leaving liver (ug/L)
      CF <- AF/VF  # concentration of PFOA in fat (ug/L)
      CVF <- CF/PF # concentration of PFOA leaving fat (ug/L)
      CK <- AK/VK  # concentration of PFOA in kidney (ug/L)
      CVK <- CK/PK  # concentration of PFOA leaving the kidney (ug/L)
      CSk <- ASk/VSk  # concentration of PFOA in skin (ug/L)
      CVSk <- CSk/PSk  # Concentration of PFOA leaving the skin (ug/L)
      CR <- AR/VR  # concentration of PFOA in rest of the body (ug/L)
      CVR <- CR/PR  # concentration of PFOA leaving the rest of the body (ug/L)
      Cfil <- Afil/Vfil # concentration of PFOA in filtrate compartment
      # CFil <- Afil/Vfil # UPDATED concentration of PFOA in filtrate compartment
      
      
      # # UPDATED SKIN & SKIN BARRIER ####
      # CSk <- ASk/VSk  # concentration of PFOA in skin (ug/L)
      # CSkb <- ASkb/VSkb # concentration of PFOA in skin barrier (ug/L)
      # 
      # 
      
      ## Plasma compartment
      RPlas <- QF*CF*FreeF+(QL+QG)*CL*FreeL+QR*CR*FreeR+QSk*CSk*FreeSk+
        QK*CK*FreeK-QCP*CA*Free #-Qfil*CA*Free 
      # Rate of PFOA amount change in plasma (ug/h)
      
      
      ## Gut compartment
      RG <- QG*CA*Free-QG*CG*FreeG+CL*FreeL*CLbiliary+Input1*DoseOn-CG*FreeG*CLfaeces
      # Rate of PFOA amount change in gut (ug/h)
      
      ## Liver compartment
      RL <- QL*CA*Free+QG*CG*FreeG-(QL+QG)*CL*FreeL-CL*FreeL*CLbiliary
      # Rate of PFOA amount change in the liver (ug/h)
      
      ## Fat compartment
      RF <- QF*(CA*Free-CF*FreeF) # Rate of PFOA amount change in fat (ug/h)
      
      ## Kidney compartment
      RK <- QK*(CA*Free-CK*FreeK)+Tm*Cfil/(Kt+Cfil)-kfil*CK*FreeK # Qfil*CK*Free was introduced to reflect clearance to filtrate compartment
      # Rate of PFOA amount change in kidney (ug/h)

      ## Filtrate compartment
      Rfil <- kfil*(CK*FreeK-Cfil)-Tm*Cfil/(Kt+Cfil) # changed from Qfil*CA*Free
      # Rate of PFOA amount change in filtrate compartment (ug/h)

      ## Storage compartment for urine
      Rdelay <- kfil*Cfil-CLurine*Adelay # ug/h

      ## Urine
      Rurine <- CLurine*Adelay # ug/h
      
      
      # ## UPDATED KIDNEY ####
      # 
      # ## Kidney compartment
      # ### Kidney: this is basically the kidney blood compartment
      # ### I think we're missing active secretion
      # RK <- QK*(CA*Free-CK*FreeK) + Vmax*CFil/Km+CFil - Qfil*CK*FreeK
      # #                   re-absorption       ultrAfiltration
      # 
      # ### Filtrate compartment: rate of formation of the filtrate(=urine) in the lumen; this is the urinary tract from Aude's lifestage equations
      # RFil <- Qfil*CK*FreeK - Vmax*CFil/Km+CFil - QUr*CFil 
      # #        ultrAfiltration  REABSORPTION      urine flow rate to the bladder
      # 
      # ### Urine: urine accumulation in the bladder; ### ??? why do we need this compartment? I would try having a simple compartment first; i.e not having an excretion compartment
      # RUr <- QUr*CFil - kUr*AUr #(ug/h)
      # #      urineflow  urine excretion
      # 
      # ### Excretion urinary: cumulative PFOA concentration in the urine
      # REx_ur <- kUr*AUr #(ug/h)
      
      
      
      ## Feaces
      
      Rfaeces <- CG*FreeG*CLfaeces #µg/h
      
      ## Skin compartment
      RSk <- QSk*(CA*Free-CSk*FreeSk)+Input2*DoseOn # Rate of PFOA amount change in skin
      
      # ## UPDATED ####
      # ## Skin barrier
      # RSkb <- Input2*DoseOn - kperm*CSkb #(ug/h)
      # ## Skin tissue
      # RSk <- kperm*CSkb + QSk*(CA*Free-CSk*FreeSk)
      
      
      ## Rest of the body
      RR <- QR*(CA*Free-CR*FreeR)
      
      # PFOAPlasTot <- APlas/(VPlasC*BW)/Free #microg/litre = ng/ml
      # PFOAUrine <- Aurine
      # PFOALiver <- AL/VL
      # PFOAKidney <- AK/VK
      
      
      return(list(c(Input1, Input2, RPlas, RG, RL, RF, RK, Rfil, Rdelay, Rurine, Rfaeces, RSk, RR))) # RFil, RUr, REx_ur,  RSkb, This will be changed to measure the uncertainty/sensitivity in other organs
    })
  }
  
  #outputs = c("PFOAPlasTot") # This will be changed to measure the uncertainty/sensitivity in other organs
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Solve the system of differential equations
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  v <- ode(y=yini, times=times, func=PBPKmodPFOA, parms=para, method="lsoda")
  
  head(v)
  
  
  PFOAamount <- as.data.frame(v)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # From amounts to concentrations
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  v[,"APlas"] <- v[,"APlas"]/(VPlasC*BW)/Free
  v[,"AG"] <- v[,"AG"]/VG
  v[,"AL"] <- v[,"AL"]/VL
  v[,"AF"] <- v[,"AF"]/VF
  v[,"AK"] <- v[,"AK"]/VK
  v[,"Afil"] <- v[,"Afil"]/Vfil
  v[,"AR"] <- v[,"AR"]/VR
  #v[, "ASkb"] <- v[, "ASkb"]/VSkb
  
  
  PFOAconc <- as.data.frame(v)
  
  PFOASerum[,(i)] <- PFOAconc$APlas
  PFOAFat[,(i)] <- PFOAconc$AF
  PFOAUrine[,(i)] <- PFOAconc$Aurine
  PFOAKidney[,(i)] <- PFOAconc$AK
  PFOALiver[,(i)] <- PFOAconc$AL
  
}

### END PBPK ###

PFOASerum <- as.data.frame(PFOASerum) # ug/L the same as ng/ml
PFOAFat <- as.data.frame(PFOAFat)
PFOAUrine <- as.data.frame(PFOAUrine)
PFOAKidney <- as.data.frame(PFOAKidney)
PFOALiver <- as.data.frame(PFOALiver)



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#### Check mass balance####
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Qbal = QCP-(QR+QL+QF+QK+QG+QSk) # Mass balance check for the cardiac output 

Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+VR)  # Mass balance check for the volumes


PFOA_bal <- sum(PFOAamount[,"Input1"]+ PFOAamount[,"Input2"]- PFOAamount[,"APlas"]- # Mass balance for PFOA
                  PFOAamount[,"AG"]-PFOAamount[,"AL"]-PFOAamount[,"AF"]-PFOAamount[,"AK"]-PFOAamount[,"AR"]-
                  PFOAamount[,"ASk"] - PFOAamount[,"Afil"]-PFOAamount[,"Aurine"]-PFOAamount[,"Adelay"]-PFOAamount[,"Afaeces"]
                )# - PFOAamount[,"ASkb"])


## MASS BALANCE & ERROR ADDED BY CHRYSANTHI ####

MB_dataf <- PFOAamount %>% 
  mutate(TOTAL = Input1 + Input2, 
         CALCULATED = rowSums(select(., APlas:AR))) %>% 
  mutate(MB = TOTAL - CALCULATED + 1) %>% 
  mutate(ER = (TOTAL - CALCULATED) / (TOTAL + 1*10^(-30)) * 100) %>% 
  mutate(time = time/(24*365)) #%>% 
# filter(time >= 25)


MB_ERROR <- MB_dataf %>% ggplot() +
  geom_line(aes(time, MB), color = "grey") +
  geom_line(aes(time, ER), color = "grey6") +
  labs(x = "Time (years)", y = "Value",
       caption = "Massbalance in grey, Error in black") +
  theme_CP()

MB_ERROR

ggsave(filename=file.path(Saveoutput, "MB_ERROR.png"),
       dpi = 300)



# add the time to the data frames

PFOASerum$time <- PFOAamount$time
PFOASerum$dose <- PFOAamount$Input1
PFOASerum$year <- PFOASerum[,"time"]/(24*365)

PFOAFat$time <- PFOAamount$time
PFOAFat$year <- PFOAFat[,"time"]/(24*365)

PFOAUrine$time <- PFOAamount$time
PFOAUrine$year <- PFOAUrine[,"time"]/(24*365)

PFOAKidney$time <- PFOAamount$time
PFOAKidney$year <- PFOAKidney[,"time"]/(24*365)

PFOALiver$time <- PFOAamount$time
PFOALiver$year <- PFOALiver[,"time"]/(24*365)






### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Plotting results####
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++####

# PFOA in plasma of one individual
Plot_PFOA_Plasm <- ggplot()+
  geom_path(data = PFOASerum, aes(x = year, y = V1))+
  theme_CP()+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  #scale_colour_hue()+
  ylab("Plasma (ng/ml)")

Plot_PFOA_Plasm


# PFOA in urine of one individual
Plot_PFOA_urine <- ggplot()+
  geom_path(data = PFOAUrine, aes(x = year, y = V1))+
  theme_CP()+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  #scale_colour_hue()+
  ylab("Urine (ng/ml)")

Plot_PFOA_urine


# PFOA in liver of one individual
Plot_PFOA_liver<- ggplot()+
  geom_path(data = PFOALiver, aes(x = year, y = V1))+
  theme_CP()+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  #scale_colour_hue()+
  ylab("Liver (ng/ml)")+
  xlab("Year")

Plot_PFOA_liver


# PFOA in kidney of one individual
Plot_PFOA_kidney <- ggplot()+
  geom_path(data = PFOAKidney, aes(x = year, y = V1))+
  theme_CP()+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  #scale_colour_hue()+
  ylab("Kidney (ng/ml)")+
  xlab("Year")

Plot_PFOA_kidney

# Combine plot from serum, urine, liver and kidney

Plot_combined_PFOA <- ggarrange(Plot_PFOA_Plasm, Plot_PFOA_urine, Plot_PFOA_liver, Plot_PFOA_kidney,
                                ncol = 2, nrow = 2)

Plot_combined_PFOA

# Save the results

ggsave(filename=file.path(Saveoutput, "Plot_combined_PFOA.png"),
       dpi = 300)
