# Westerhout et al., 2024 combined PFOA and PFOS PBK model
# Parameter values are obtained from [EFSA 2020] unless specified otherwise
# Remark: the doses are given in ug, the volumes in L, so concentrations are in ug/L
# From ug/L to uM is by dividing by MW

rm(list=ls()) # to clear out the global environment

# Set working directory
# setwd("Enter path here")

library(deSolve)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#### Settings ----
TSTART <- 0
TSTOP <- 365*50 # simulation duration (days); 50 years
exposure_duration = 365*50 # Duration of exposure (d); 50 years
DT <- 1 # time step (days)
TIME <- seq(TSTART,TSTOP,by=DT)
AGE <- 0 # to use as a starting point, e.g. 20 years of age instead of birth
age_start_derm_exp <- 14 # dermal and inhalation exposure starting age (e.g. occupational, or personal care products; exposure may start at a later age than 0)
age_stop_derm_exp <- 50 # dermal and inhalation exposure stopping age (e.g. occupational, or personal care products; exposure may stop at an earlier age than 80)
age_start_inh_exp <- 0 # dermal and inhalation exposure starting age (e.g. occupational, or personal care products; exposure may start at a later age than 0)
age_stop_inh_exp <- 50 # dermal and inhalation exposure stopping age (e.g. occupational, or personal care products; exposure may stop at an earlier age than 80)

Variables_df <- as.data.frame(list(TIME = TIME))
Variables_df <- Variables_df %>%
  mutate(day = TIME+1) %>%
  mutate(dayoftheweek = rep(seq(1:7),length.out = n())) %>%
  mutate(month = ceiling(day/(365/12))-(((ceiling(day/365))-1)*12)) %>%
  mutate(year = ceiling(day/(365))) %>%
  mutate(age = AGE + (TIME/(365)))

# In the Westerhout model, in comparison to the EFSA 2020 model, the plasma compartment is split into arterial and venous blood.
# In addition, the lung compartment is new and the skin compartment is now included.
# We did not know the appropriate starting values for these compartments before we did a run,
# so we first set them to 0 (except for the arterial and venous blood compartments, which we based on the plasma value). 
# The starting value for the rest compartment was also not given and was set to 0.

# The output of the the run with EFSA TWI values "Westerhout 2024 - Results.csv" was manually modified
# so that the initial amounts in the lungs, skin and remaining compartments at time 0 are representative 
# and the mass balance was in check.
# After doing so, the file was saved as 'Westerhout 2024 - EFSA TWI background amounts.csv'

# In the future runs, load the file here and use this for background amounts:
# Background_PFAS_amounts <- read.csv2(file = 'Westerhout 2024 - EFSA TWI background amounts.csv',na.strings=c(""," ","NA"),stringsAsFactors = FALSE)
# Background_PFAS_amounts[] <- lapply(Background_PFAS_amounts, as.numeric)
# Note that this file is also used to determine the starting values in case you start the simulation at a later age
# instead of birth (relevant for the occupational setting).

#### Physiological parameters ----
# Physiological parameters (from Brown, et al. 1997)
#fractional blood flows
QCC = 300 # 12.5*24   # Cardiac blood output (L/d/kg^0.75)
QFC = 0.052  # Fraction cardiac output going to fat
QLC = 0.069   # Fraction cardiac output going to liver, through hepatic artery
QKC = 0.175  # Fraction cardiac output going to kidney
QSkC = 0.058  # Fraction cardiac output going to skin
QGC = 0.181  # Fraction of cardiac output going to gut and in the liver via portal artery
# Not used #QfilC = 0.035  # Fraction cardiac output to the filtrate compartment (20% of kidney blood flow)

#fractional tissue volumes
VLC = 0.026  # Fraction liver volume
VFC = 0.214  # Fraction fat volume
VKC = 0.004  # Fraction kidney volume
VfilC = 0.0004  # Fraction filtrate compartment volume (10% of kidney volume) 
VGC = 0.0171  # Fraction gut volume 
VlunC = 0.007 # Fraction lung volume (ICRP 89)
VPlasC = 0.0428 # Fraction plasma volume (58% of blood) 
VvenC = 0.61*VPlasC # Fraction of venous plasma volume [Brown 1997]
VartC = 0.39*VPlasC # Fraction of arterial plasma volume [Brown 1997]
Htc = 0.44 # hematocrit 

# Mother-child parameters
BWbirth=3.68 # Body weight at birth in kg new add opinion 2020 

# Skin parameters 
# Dermal exposure 
# Dermexpo = 0 # 0 = NO, 1 = YES
# Dermconc = 0.0 # Dermal concentration (mg/mL) 
# Dermvol = 0.001 # Dermal exposure volume (mL); cannot be 0
# Dermdose = Dermconc*Dermvol*1000 # (ug) 
# Skinarea = 972 # Exposed area on skin (cm^2); surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
Skinthickness = 0.1 # Skin thickness (cm) 

fss <- 0.005 # fraction palm of hands (Sheridan et al., 1995, Rhodes et al., 2013)
hsurf <- 0.01 # applied layer thickness in cm, value is assumed to be 0.1 mm
hsc <- 0.0015 # Stratum corneum thickness (cm) # [Krüse 2007]
hve <- 0.0100 # Viable epidermis thickness (cm) # [Krüse 2007]
hcell <- 0.00005 # Blood vessel wall thickness (cm) # [Burton 1954] = 0.5 um capillary wall thickness

ffatsc <- 0.05 # Fraction of fat in stratum corneum [Polak 2012]
ffatve <- 0.02 # Fraction of fat in viable epidermis [Polak 2012]
ffatepi <- 0.02 # Fraction of fat in epidermis [Polak 2012]
ffatbl <- 0.007 # Fraction of fat in blood [Polak 2012]

Kpcell <- 0.93 # cm/h; Krüse model is normalized to 1 cm^2 # Permeability coefficient from arterial wall into the blood compartment [Krüse 2007]

#### PFOA ----
# Physicochemical properties
MW_PFOA = 414.07
logP_PFOA = 4.81
VP_PFOA = 2.34 # 10^0.37 # Pa [Zhang 2021]; 1 mmHg = 133.322368 Pa; 0.15 mmHg (https://haz-map.com/Agents/6596) = 20 Pa

#Kp_PFOA = 10^(0.71*logP_PFOA-0.0061*MW_PFOA-6.3)
Kscve_PFOA = ((1-ffatve)+(ffatve*(10^logP_PFOA)))/((1-ffatsc)+(ffatsc*(10^logP_PFOA)))
Kver_PFOA = ((1-ffatbl)+(ffatbl*(10^logP_PFOA)))/((1-ffatepi)+(ffatepi*(10^logP_PFOA)))
# Kpsc_PFOA = 10^(0.74*logP_PFOA - 0.006*MW_PFOA - 2.8) # Stratum corneum permeability coefficient (cm/h) # [Polak 2012; KrÃ¼se 2007; Bunge and Cleek 1995; Cleek and Bunge 1993; Potts and Guy 1992]
Kpsc_PFOA = 0.000088 # [Franko 2011]
# Kpsc_PFOA = 0.000000949 # [Fasano 2005]
# Kpve_PFOA <- 2.6/(MW_PFOA^0.5) # Viable epidermis permeability coefficient (cm/h) # [Polak 2012 -> EPA 2004]
Kpve_PFOA <- 0.000088 # [Franko 2011]
# Kpve_PFOA <- 0.000000949 # [Fasano 2005]

Pab_PFOA <- 1/(10^(6.96-(1.04*log10(VP_PFOA)) - 0.533*logP_PFOA - 0.00495*MW_PFOA)) # [Buist 2012]

# Oral exposure via food
Oralexpo_PFOA = 0.000187 # µg/kg/day [EFSA 2020; page 143]

# Drinking water exposure 
Drinkrate = 13 # Drinking water rate (mL/kg/day)
Drinkconc_PFOA = 0 # Drinking water concentration (ug/L or ppb) 

# Inhalation exposure
Cinh_PFOA = 0 # ug/L; 0.526 ug/m3 [Nilsson 2013]
Cinh_82FTOH = 0 # ug/L; 114 ug/m3 [Nilsson 2013]
Biotransformation = 0.003 # 8:2-FTOH to PFOA [Gomis 2016]
Cinh_PFOA_total_day = Cinh_PFOA + (Cinh_82FTOH*Biotransformation) # exposure spread across a 24h day
# In case of a temporary inhalation exposure during the day, multiply Cinh_PFOA with the fraction of the day that exposure took place (e.g. 6h = 25% -> Cinh_PFOA_total_day = Cinh_PFOA*0.25)

# Dermal exposure
Cdermal_PFOA = 0 # ug/kg = 0.680 ug/g [Freberg 2010]
Cdermal_PFOA_total_day = Cdermal_PFOA # exposure spread across a 24h day
# In case of a temporary dermal exposure during the day, multiply Cdermal_PFOA with the fraction of the day that exposure took place (e.g. 6h = 25% -> Cdermal_PFOA_total_day = Cdermal_PFOA*0.25)

# Chemical-specific parameters (PFOA) 
Tmc_PFOA = 144000 # 6000*24   # Maximum resorption rate 
Kt_PFOA = 55   # Resorption affinity# same as monkey
kurinec_PFOA = 0.0072 # 0.0003*24 # urinary elimination rate constant (/d/kg^-0.25)# estimated 
Free_PFOA = 0.02  # Free fraction of PFOA in plasma# same as monkey 
# Partition coefficients from Harada, et al 2005 
PL_PFOA = 2.2   # Liver/plasma partition coefficient of PFOA
PF_PFOA = 0.04   # Fat/plasma partition coefficient of PFOA
PK_PFOA = 1.05   # Kidney/plasma partition coefficient of PFOA
PSk_PFOA = 0.1   # Skin/plasma partition coefficient of PFOA
PR_PFOA = 0.12   # Rest of the body/plasma partition coefficient of PFOA
PG_PFOA = 0.05   # Gut/blood plasma coefficient of PFOA  
PLun_PFOA = 1.27   # Lung/plasma partition coefficient of PFOA; [Fabrega 2014]

# Free fraction of chemical in tissues 
FreeL_PFOA = Free_PFOA/PL_PFOA                      #liver 
FreeF_PFOA = Free_PFOA/PF_PFOA                     #fat 
FreeK_PFOA = Free_PFOA/PK_PFOA                     #kidney 
FreeSk_PFOA = Free_PFOA/PSk_PFOA                #skin 
FreeR_PFOA = Free_PFOA/PR_PFOA                   #rest of tissues 
FreeG_PFOA = Free_PFOA/PG_PFOA                    #gut 
FreeLun_PFOA = Free_PFOA/PLun_PFOA                    #lung 

# Mother-child parameters
# PFOA
PFOAmaternal = 2.0 # maternal concentration ng/mL at delivery [EFSA opinion 2020, p368]
PT_PFOA= 0.74 # placenta transfer for PFOA, new add opinion 2020 
Ratio_PFOA= 0.03 # milk concentration/maternal serum concentration for PFOA during breastfeeding, new add opinion 2020 
DECLINE_PFOA = 0.077 # decline of PFOA in milk was 7.7% per month, new add opinion 2020 

# Oral exposure via the mother
PFOaMilkconcentration = PFOAmaternal*Ratio_PFOA # initial Milk concentration at birth ng/mL 
Milkconsumption = 0.8 # milkconsumption L per day 
Intakemilka_PFOA = PFOaMilkconcentration*Milkconsumption# initial intake via breastfeeding first month
Intakemilkb_PFOA = Intakemilka_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeedingsecondmonth 
Intakemilkc_PFOA = Intakemilkb_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding3month
Intakemilkd_PFOA = Intakemilkc_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding4month 
Intakemilke_PFOA = Intakemilkd_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding5month
Intakemilkf_PFOA = Intakemilke_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding6month 
Intakemilkg_PFOA = Intakemilkf_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding7month 
Intakemilkh_PFOA = Intakemilkg_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding8month 
Intakemilki_PFOA = Intakemilkh_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding9month 
Intakemilkj_PFOA = Intakemilki_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding10month 
Intakemilkk_PFOA = Intakemilkj_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding11month 
Intakemilkl_PFOA = Intakemilkk_PFOA*(1-DECLINE_PFOA)#intakeviabreastfeeding12month 

CONCbirth_PFOA = PFOAmaternal*PT_PFOA # PFOA concentration at birth, new add in 2020 
# Due to the mass balance correction, the free fraction is removed from the following equation
# APlasbirth_PFOA = CONCbirth_PFOA*Free_PFOA*VPlasC*BWbirth # initial amount in plasma at birth
APlasbirth_PFOA = CONCbirth_PFOA*VPlasC*BWbirth # initial amount in plasma at birth
# The following starting values are slightly different from EFSA 2020 due to the addition of the skin and lung compartments
Avenbirth_PFOA = APlasbirth_PFOA*0.61 # New in Westerhout 2024; plasma divided into venous and arterial plasma [Brown 1997]
Aartbirth_PFOA = APlasbirth_PFOA*0.39 # New in Westerhout 2024; plasma divided into venous and arterial plasma [Brown 1997]
# Next, the free fraction is incorporated into the ratios at birth
Alunbirth_PFOA = APlasbirth_PFOA*Free_PFOA*8.965856746 # New in Westerhout 2024
Ascbirth_PFOA = APlasbirth_PFOA*0 # New in Westerhout 2024
Atransbirth_PFOA = APlasbirth_PFOA*0 # New in Westerhout 2024
Avebirth_PFOA = APlasbirth_PFOA*0 # New in Westerhout 2024
ASkbirth_PFOA = APlasbirth_PFOA*Free_PFOA*6.88526319 # New in Westerhout 2024
AGbirth_PFOA = APlasbirth_PFOA*Free_PFOA*0.99929853 # initial amount in gut at birth
ALbirth_PFOA = APlasbirth_PFOA*Free_PFOA*66.8360764 # initial amount in liver at birth
AFbirth_PFOA = APlasbirth_PFOA*Free_PFOA*10.000151 # initial amount in fat at birth
AKbirth_PFOA = APlasbirth_PFOA*Free_PFOA*5.88575746 # initial amount in kidney at birth
ARbirth_PFOA = APlasbirth_PFOA*Free_PFOA*74.9108792 # initial amount in rest of the body at birth
Afilbirth_PFOA = APlasbirth_PFOA*Free_PFOA*0.0000209606 # initial amount filtrate at birth
Adelaybirth_PFOA = APlasbirth_PFOA*Free_PFOA*42.27910714 # initial amount in storage at birth
Aurinebirth_PFOA = APlasbirth_PFOA*Free_PFOA*769.1347494 # initial amount in urine at birth

# PFOA age-dependent starting values
# APlas_PFOA_background = Background_PFAS_amounts$APlas_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aven_PFOA_background = Background_PFAS_amounts$Aven_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aart_PFOA_background = Background_PFAS_amounts$Aart_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AG_PFOA_background = Background_PFAS_amounts$AG_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AL_PFOA_background = Background_PFAS_amounts$AL_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AF_PFOA_background = Background_PFAS_amounts$AF_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AK_PFOA_background = Background_PFAS_amounts$AK_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Afil_PFOA_background = Background_PFAS_amounts$Afil_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Adelay_PFOA_background = Background_PFAS_amounts$Adelay_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aurine_PFOA_background = Background_PFAS_amounts$Aurine_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# ASk_PFOA_background = Background_PFAS_amounts$ASk_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Alun_PFOA_background = Background_PFAS_amounts$Alun_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AR_PFOA_background = Background_PFAS_amounts$AR_PFOA[which(Background_PFAS_amounts$age==round(AGE, digits=1))]

# PFOA starting values at age 0
Aven_PFOA_background = Avenbirth_PFOA
Aart_PFOA_background = Aartbirth_PFOA
AG_PFOA_background = AGbirth_PFOA
AL_PFOA_background = ALbirth_PFOA
AF_PFOA_background = AFbirth_PFOA
AK_PFOA_background = AKbirth_PFOA
Afil_PFOA_background = Afilbirth_PFOA
Adelay_PFOA_background = Adelaybirth_PFOA
Aurine_PFOA_background = Aurinebirth_PFOA
Asc_PFOA_background = Ascbirth_PFOA
Atrans_PFOA_background = Atransbirth_PFOA
Ave_PFOA_background = Avebirth_PFOA
ASk_PFOA_background = ASkbirth_PFOA
Alun_PFOA_background = Alunbirth_PFOA
AR_PFOA_background = ARbirth_PFOA

#### PFOS ----
# Physicochemical properties
MW_PFOS = 500.13
logP_PFOS = 4.49
VP_PFOS = 0.27 # Pa; 1 mmHg = 133.322368 Pa; 0.002 mmHg (https://pubchem.ncbi.nlm.nih.gov/source/hsdb/7099#section=Environmental-Fate-%26-Exposure)

#Kp_PFOS = 10^(0.71*logP_PFOS-0.0061*MW_PFOS-6.3)
Kscve_PFOS = ((1-ffatve)+(ffatve*(10^logP_PFOS)))/((1-ffatsc)+(ffatsc*(10^logP_PFOS)))
Kver_PFOS = ((1-ffatbl)+(ffatbl*(10^logP_PFOS)))/((1-ffatepi)+(ffatepi*(10^logP_PFOS)))
Kpsc_PFOS = 10^(0.74*logP_PFOS - 0.006*MW_PFOS - 2.8) # Stratum corneum permeability coefficient (cm/h) # [Polak 2012; KrÃ¼se 2007; Bunge and Cleek 1995; Cleek and Bunge 1993; Potts and Guy 1992]
Kpve_PFOS <- 2.6/(MW_PFOS^0.5) # Viable epidermis permeability coefficient (cm/h) # [Polak 2012 -> EPA 2004]

Pab_PFOS <- 1/(10^(6.96-(1.04*log10(VP_PFOS)) - 0.533*logP_PFOS - 0.00495*MW_PFOS)) # [Buist 2012]

# Oral exposure via food
Oralexpo_PFOS = 0.000444 # µg/kg/day [EFSA 2020; page 143]

# Drinking water exposure 
Drinkconc_PFOS = 0 # Drinking water concentration (ug/L or ppb) 

# Inhalation exposure
Cinh_PFOS = 0 # mg/L; 0 ug/m3
Cinh_PFOS_total_day = Cinh_PFOS # exposure spread across a 24h day
# In case of a temporary inhalation exposure during the day, multiply Cinh_PFOS with the fraction of the day that exposure took place (e.g. 6h = 25% -> Cinh_PFOS_total_day = Cinh_PFOS*0.25)

# Dermal exposure
Cdermal_PFOS = 0 # ug/kg = 0 ug/g [Freberg 2010]
Cdermal_PFOS_total_day = Cdermal_PFOS # exposure spread across a 24h day
# In case of a temporary dermal exposure during the day, multiply Cdermal_PFOS with the fraction of the day that exposure took place (e.g. 6h = 25% -> Cdermal_PFOS_total_day = Cdermal_PFOS*0.25)

# Chemical-specific parameters (PFOS) 
Tmc_PFOS = 84000 # 3500*24   # Maximum resorption rate 
Kt_PFOS = 23   # Resorption affinity# same as monkey
kurinec_PFOS = 0.024 # 0.001*24 # urinary elimination rate constant (/d/kg^-0.25)# estimated 
Free_PFOS = 0.025  # Free fraction of PFOS in plasma# same as monkey 
# Partition coefficients from Harada, et al 2005 
PL_PFOS = 3.72   # Liver/plasma partition coefficient of PFOS
PF_PFOS = 0.14   # Fat/plasma partition coefficient of PFOS
PK_PFOS = 0.8   # Kidney/plasma partition coefficient of PFOS
PSk_PFOS = 0.29   # Skin/plasma partition coefficient of PFOS
PR_PFOS = 0.2   # Rest of the body/plasma partition coefficient of PFOS
PG_PFOS = 0.57   # Gut/blood plasma coefficient of PFOS  
PLun_PFOS = 0.15   # Lung/plasma partition coefficient of PFOS; [Fabrega 2014]

# Free fraction of chemical in tissues 
FreeL_PFOS = Free_PFOS/PL_PFOS                      #liver 
FreeF_PFOS = Free_PFOS/PF_PFOS                     #fat 
FreeK_PFOS = Free_PFOS/PK_PFOS                     #kidney 
FreeSk_PFOS = Free_PFOS/PSk_PFOS                #skin 
FreeR_PFOS = Free_PFOS/PR_PFOS                   #rest of tissues 
FreeG_PFOS = Free_PFOS/PG_PFOS                    #gut 
FreeLun_PFOS = Free_PFOS/PLun_PFOS                    #lung 

# Mother-child parameters
# PFOS
PFOSmaternal = 4.89 # maternal concentration ng/mL at delivery [EFSA opinion 2020, p.372 (different from p.152)]
PT_PFOS= 0.36 # placenta transfer for PFOS, new add opinion 2020 
Ratio_PFOS= 0.015 # milk concentration/maternal serum concentration for PFOS during breastfeeding, new add opinion 2020 
DECLINE_PFOS = 0.031 # decline of PFOS in milk was 3.1% per month, new add opinion 2020 

# Oral exposure via the mother
PFOSMilkconcentration = PFOSmaternal*Ratio_PFOS # initial Milk concentration at birth ng/mL 
Milkconsumption = 0.8 # milkconsumption L per day 
Intakemilka_PFOS = PFOSMilkconcentration*Milkconsumption# initial intake via breastfeeding first month
Intakemilkb_PFOS = Intakemilka_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeedingsecondmonth 
Intakemilkc_PFOS = Intakemilkb_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding3month
Intakemilkd_PFOS = Intakemilkc_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding4month 
Intakemilke_PFOS = Intakemilkd_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding5month
Intakemilkf_PFOS = Intakemilke_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding6month 
Intakemilkg_PFOS = Intakemilkf_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding7month 
Intakemilkh_PFOS = Intakemilkg_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding8month 
Intakemilki_PFOS = Intakemilkh_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding9month 
Intakemilkj_PFOS = Intakemilki_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding10month 
Intakemilkk_PFOS = Intakemilkj_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding11month 
Intakemilkl_PFOS = Intakemilkk_PFOS*(1-DECLINE_PFOS)#intakeviabreastfeeding12month 

CONCbirth_PFOS = PFOSmaternal*PT_PFOS # PFOS concentration at birth, new add in 2020 
# Due to the mass balance correction, the free fraction is removed from the following equation
# APlasbirth_PFOS = CONCbirth_PFOS*Free_PFOS*VPlasC*BWbirth # initial amount in plasma at birth
APlasbirth_PFOS = CONCbirth_PFOS*VPlasC*BWbirth # initial amount in plasma at birth
# The following starting values are slightly different from EFSA 2020 due to the addition of the skin and lung compartments
Avenbirth_PFOS = APlasbirth_PFOS*0.61 # New in Westerhout 2024; plasma divided into venous and arterial plasma [Brown 1997]
Aartbirth_PFOS = APlasbirth_PFOS*0.39 # New in Westerhout 2024; plasma divided into venous and arterial plasma [Brown 1997]
# Next, the free fraction is incorporated into the ratios at birth
Alunbirth_PFOS = APlasbirth_PFOS*Free_PFOS*7.09776214 # New in Westerhout 2024
Ascbirth_PFOS = APlasbirth_PFOS*0 # New in Westerhout 2024
Atransbirth_PFOS = APlasbirth_PFOS*0 # New in Westerhout 2024
Avebirth_PFOS = APlasbirth_PFOS*0 # New in Westerhout 2024
ASkbirth_PFOS = APlasbirth_PFOS*Free_PFOS*16.44603423 # New in Westerhout 2024
AGbirth_PFOS = APlasbirth_PFOS*Free_PFOS*9.112506818 # initial amount in gut at birth 
ALbirth_PFOS = APlasbirth_PFOS*Free_PFOS*90.41427413 # initial amount in liver at birth 
AFbirth_PFOS = APlasbirth_PFOS*Free_PFOS*27.99955778 # initial amount in fat at birth 
AKbirth_PFOS = APlasbirth_PFOS*Free_PFOS*3.587821811 # initial amount in kidney at birth 
ARbirth_PFOS = APlasbirth_PFOS*Free_PFOS*100.1289877 # initial amount in rest of the body at birth 
Afilbirth_PFOS = APlasbirth_PFOS*Free_PFOS*0.0000150221 # initial amount filtrate at birth 
Adelaybirth_PFOS = APlasbirth_PFOS*Free_PFOS*9.156014881 # initial amount in storage at birth 
Aurinebirth_PFOS = APlasbirth_PFOS*Free_PFOS*565.1438095 # initial amount in urine at birth 

# PFOS age-dependent starting values
# APlas_PFOS_background = Background_PFAS_amounts$APlas_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aven_PFOS_background = Background_PFAS_amounts$Aven_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aart_PFOS_background = Background_PFAS_amounts$Aart_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AG_PFOS_background = Background_PFAS_amounts$AG_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AL_PFOS_background = Background_PFAS_amounts$AL_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AF_PFOS_background = Background_PFAS_amounts$AF_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AK_PFOS_background = Background_PFAS_amounts$AK_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Afil_PFOS_background = Background_PFAS_amounts$Afil_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Adelay_PFOS_background = Background_PFAS_amounts$Adelay_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Aurine_PFOS_background = Background_PFAS_amounts$Aurine_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# ASk_PFOS_background = Background_PFAS_amounts$ASk_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# Alun_PFOS_background = Background_PFAS_amounts$Alun_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]
# AR_PFOS_background = Background_PFAS_amounts$AR_PFOS[which(Background_PFAS_amounts$age==round(AGE, digits=1))]

# PFOS starting values at age 0
Aven_PFOS_background = Avenbirth_PFOS
Aart_PFOS_background = Aartbirth_PFOS
AG_PFOS_background = AGbirth_PFOS
AL_PFOS_background = ALbirth_PFOS
AF_PFOS_background = AFbirth_PFOS
AK_PFOS_background = AKbirth_PFOS
Afil_PFOS_background = Afilbirth_PFOS
Adelay_PFOS_background = Adelaybirth_PFOS
Aurine_PFOS_background = Aurinebirth_PFOS
Asc_PFOS_background = Ascbirth_PFOS
Atrans_PFOS_background = Atransbirth_PFOS
Ave_PFOS_background = Avebirth_PFOS
ASk_PFOS_background = ASkbirth_PFOS
Alun_PFOS_background = Alunbirth_PFOS
AR_PFOS_background = ARbirth_PFOS

#### Time-dependent variables ----
Variables_df <- Variables_df %>%
  mutate(BW = BWbirth + 4.47*age - 0.093*age^2 + 0.00061*age^3) %>%
  mutate(SkinTarea = 9.1*((BW*1000)**0.666)) %>% # Total area of skin (cm^2) 
  mutate(kurine_PFOA = kurinec_PFOA*BW**(-0.25)) %>%
  mutate(kurine_PFOS = kurinec_PFOS*BW**(-0.25)) %>%
  # Oral exposure
  mutate(Oraldose_PFOA = if_else(age < 0.083, Intakemilka_PFOA/BW,
                                 if_else(age >= 0.083 & age < 0.167, Intakemilkb_PFOA/BW,
                                         if_else(age >= 0.167 & age < 0.250, Intakemilkc_PFOA/BW,
                                                 if_else(age >= 0.250 & age < 0.333, Intakemilkd_PFOA/BW,
                                                         if_else(age >= 0.333 & age < 0.417, Intakemilke_PFOA/BW,
                                                                 if_else(age >= 0.417 & age < 0.500, Intakemilkf_PFOA/BW,
                                                                         if_else(age >= 0.500 & age < 0.583, Intakemilkg_PFOA/BW,
                                                                                 if_else(age >= 0.583 & age < 0.667, Intakemilkh_PFOA/BW,
                                                                                         if_else(age >= 0.667 & age < 0.750, Intakemilki_PFOA/BW,
                                                                                                 if_else(age >= 0.750 & age < 0.833, Intakemilkj_PFOA/BW,
                                                                                                         if_else(age >= 0.833 & age < 0.917, Intakemilkk_PFOA/BW,
                                                                                                                 if_else(age >= 0.917 & age < 1, Intakemilkl_PFOA/BW,
                                                                                                                         Oralexpo_PFOA))))))))))))) %>% # (ug/day)
  mutate(Oraldose_PFOA = if_else(day > exposure_duration,0,Oraldose_PFOA*BW)) %>% # stop exposure after exposure_duration
  mutate(Oraldose_PFOS = if_else(age < 0.083, Intakemilka_PFOS/BW,
                                 if_else(age >= 0.083 & age < 0.167, Intakemilkb_PFOS/BW,
                                         if_else(age >= 0.167 & age < 0.250, Intakemilkc_PFOS/BW,
                                                 if_else(age >= 0.250 & age < 0.333, Intakemilkd_PFOS/BW,
                                                         if_else(age >= 0.333 & age < 0.417, Intakemilke_PFOS/BW,
                                                                 if_else(age >= 0.417 & age < 0.500, Intakemilkf_PFOS/BW,
                                                                         if_else(age >= 0.500 & age < 0.583, Intakemilkg_PFOS/BW,
                                                                                 if_else(age >= 0.583 & age < 0.667, Intakemilkh_PFOS/BW,
                                                                                         if_else(age >= 0.667 & age < 0.750, Intakemilki_PFOS/BW,
                                                                                                 if_else(age >= 0.750 & age < 0.833, Intakemilkj_PFOS/BW,
                                                                                                         if_else(age >= 0.833 & age < 0.917, Intakemilkk_PFOS/BW,
                                                                                                                 if_else(age >= 0.917 & age < 1, Intakemilkl_PFOS/BW,
                                                                                                                         Oralexpo_PFOS))))))))))))) %>% # (ug/day)
  mutate(Oraldose_PFOS = if_else(day > exposure_duration,0,Oraldose_PFOS*BW)) %>% # stop exposure after exposure_duration
  mutate(Drinkdose_PFOA = (Drinkconc_PFOA*Drinkrate/1000)*BW) %>%   # (ug/day) 
  mutate(Drinkdose_PFOA = if_else(day > exposure_duration,0,Drinkdose_PFOA)) %>% # stop exposure after exposure_duration
  mutate(Drinkdose_PFOS = (Drinkconc_PFOS*Drinkrate/1000)*BW) %>%   # (ug/day) 
  mutate(Drinkdose_PFOS = if_else(day > exposure_duration,0,Drinkdose_PFOS)) %>% # stop exposure after exposure_duration
  
  # Inhalation exposure
  # Daily inhalation exposure
  mutate(Inhalation_PFOA = Cinh_PFOA_total_day) %>%
  mutate(Inhalation_PFOS = Cinh_PFOS_total_day) %>%
  
  # In case of age-dependent inhalation exposure, set the value to 0 for the ages that are not included
  mutate(Inhalation_PFOA = if_else(age >= age_start_inh_exp & age <= age_stop_inh_exp,Inhalation_PFOA,0)) %>%
  mutate(Inhalation_PFOS = if_else(age >= age_start_inh_exp & age <= age_stop_inh_exp,Inhalation_PFOS,0)) %>%
  
  # Blood flows
  mutate(QC = QCC*BW**0.75) %>% # Cardiac output (L/d) 
  mutate(QCP = QC*(1-Htc)) %>% # adjust for plasma flow 
  mutate(QL = QLC*QCP) %>% # Plasma flow to liver (L/d) 
  mutate(QF = QFC*QCP) %>% # Plasma flow to fat (L/d)
  mutate(QK = QKC*QCP) %>% # Plasma flow to kidney (L/d) 
  mutate(Qfil = 0.2*QK) %>% # Plasma flow to filtrate compartment (L/d)# 20% of QK 
  mutate(QG = QGC*QCP) %>% # Plasma flow to gut (L/d) 
  mutate(QSk = QSkC*QCP) %>% #*(Skinarea/SkinTarea)*Dermexpo) %>% #plasma flow to skin 
  mutate(QR = QCP - QL - QF - QK - QG - QSk) %>% # Plasma flow to rest of the body (L/d)
  mutate(Qbal = QCP - (QR+QL+QF+QK+QG+QSk)) %>% # balance check--better be 0 
  
  mutate(Qp = ((((1400-190)*(age)^2.5)/((age)^2.5 + 50)) + 190)*24) %>% # L/d [ICRP 89] manual curve fit of age dependency
  
  # Organ volumes
  mutate(VL = VLC*BW) %>% # Liver volume (L) 
  mutate(VF = VFC*BW) %>% # Fat volume (L) 
  mutate(VK = VKC*BW) %>% # Kidney volume (L) 
  mutate(Vfil = VfilC*BW) %>% # Filtrate compartment volume (L) 
  mutate(VG = VGC*BW) %>% # Gut volume (L) 
  mutate(VPlas = VPlasC*BW) %>% # Plasma volume (L) 
  mutate(Vart_plas = 0.39*VPlas) %>% # [Brown 1997]
  mutate(Vven_plas = 0.61*VPlas) %>% # [Brown 1997]
  mutate(VSk = (SkinTarea*Skinthickness)/1000) %>% # Skin volume (L) 
  mutate(Vlun = VlunC*BW) %>% # [ICRP 89]
  mutate(VFRC = 0.030*BW) %>% # Functional residual capacity # https://en.wikipedia.org/wiki/Functional_residual_capacity
  mutate(VT = 0.007*BW) %>% # Tidal volume # https://en.wikipedia.org/wiki/Tidal_volume
  mutate(Valv = VFRC + 0.5*VT) %>% # Alveolar volume
  
  mutate(VR = 0.84*BW - VL - VF - VK - Vfil - VG - VPlas - VSk - Vlun) %>% # Rest of the body volume (L) 
  mutate(Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+Vlun)) %>% # Balance check--better be 0 
  
  mutate(Vsurf = fss*SkinTarea*hsurf/1000) %>% # volume of dermal layer on exposed skin surface (L) 1 m^2 = 10,000 cm^2 * cm = cm^3 = ml / 1000 = L
  mutate(Dermaldose_PFOA_day = Cdermal_PFOA_total_day*Vsurf) %>%
  mutate(Dermaldose_PFOS_day = Cdermal_PFOS_total_day*Vsurf) %>%
  
  # Dermal exposure
  mutate(Csurf_PFOA = Cdermal_PFOA_total_day) %>%
  mutate(Csurf_PFOS = Cdermal_PFOS_total_day) %>%
  
  # In case of age-dependent dermal exposure, set the value to 0 for the ages that are not included
  mutate(Csurf_PFOA = if_else(age >= age_start_derm_exp & age <= age_stop_derm_exp,Csurf_PFOA,0)) %>%
  mutate(Csurf_PFOS = if_else(age >= age_start_derm_exp & age <= age_stop_derm_exp,Csurf_PFOS,0)) %>%
  
  mutate(Vsc = fss*SkinTarea*hsc/1000) %>% # volume of exposed stratum corneum (L) cm^2 * cm = cm^3 = ml / 1000 = L
  mutate(Vve = fss*SkinTarea*hve/1000) %>% # volume of exposed viable epidermis (L) cm^2 * cm = cm^3 = ml / 1000 = L
  
  mutate(CLsc_PFOA = 24*Kpsc_PFOA*fss*SkinTarea/1000) %>% # stratum corneum clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLsc_PFOS = 24*Kpsc_PFOS*fss*SkinTarea/1000) %>% # stratum corneum clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLve_PFOA = 24*Kpve_PFOA*fss*SkinTarea/1000) %>% # viable epidermis clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLve_PFOS = 24*Kpve_PFOS*fss*SkinTarea/1000) %>% # viable epidermis clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLcell = 24*Kpcell*fss*SkinTarea/1000) %>% # clearance from arterial wall into the blood compartment (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  
  mutate(Tm_PFOA = Tmc_PFOA*BW**0.75) %>% #transporter maximum 
  mutate(Tm_PFOS = Tmc_PFOS*BW**0.75) #transporter maximum 

varkurine_PFOA <- approxfun(Variables_df$TIME, Variables_df$kurine_PFOA, rule = 2)
varkurine_PFOS <- approxfun(Variables_df$TIME, Variables_df$kurine_PFOS, rule = 2)

varQCP <- approxfun(Variables_df$TIME, Variables_df$QCP, rule = 2)
varQG <- approxfun(Variables_df$TIME, Variables_df$QG, rule = 2)
varQL <- approxfun(Variables_df$TIME, Variables_df$QL, rule = 2)
varQF <- approxfun(Variables_df$TIME, Variables_df$QF, rule = 2)
varQK <- approxfun(Variables_df$TIME, Variables_df$QK, rule = 2)
varQfil <- approxfun(Variables_df$TIME, Variables_df$Qfil, rule = 2)
varQSk <- approxfun(Variables_df$TIME, Variables_df$QSk, rule = 2)
varQR <- approxfun(Variables_df$TIME, Variables_df$QR, rule = 2)

varQp <- approxfun(Variables_df$TIME, Variables_df$Qp, rule = 2)

varVPlas <- approxfun(Variables_df$TIME, Variables_df$VPlas, rule = 2)
varVart_plas <- approxfun(Variables_df$TIME, Variables_df$Vart_plas, rule = 2)
varVven_plas <- approxfun(Variables_df$TIME, Variables_df$Vven_plas, rule = 2)
varVG <- approxfun(Variables_df$TIME, Variables_df$VG, rule = 2)
varVL <- approxfun(Variables_df$TIME, Variables_df$VL, rule = 2)
varVF <- approxfun(Variables_df$TIME, Variables_df$VF, rule = 2)
varVK <- approxfun(Variables_df$TIME, Variables_df$VK, rule = 2)
varVfil <- approxfun(Variables_df$TIME, Variables_df$Vfil, rule = 2)
varVSk <- approxfun(Variables_df$TIME, Variables_df$VSk, rule = 2)
varVsurf <- approxfun(Variables_df$TIME, Variables_df$Vsurf, rule = 2)
varVsc <- approxfun(Variables_df$TIME, Variables_df$Vsc, rule = 2)
varVve <- approxfun(Variables_df$TIME, Variables_df$Vve, rule = 2)
varCLsc_PFOA <- approxfun(Variables_df$TIME, Variables_df$CLsc_PFOA, rule = 2)
varCLsc_PFOS <- approxfun(Variables_df$TIME, Variables_df$CLsc_PFOS, rule = 2)
varCLve_PFOA <- approxfun(Variables_df$TIME, Variables_df$CLve_PFOA, rule = 2)
varCLve_PFOS <- approxfun(Variables_df$TIME, Variables_df$CLve_PFOS, rule = 2)
varCLcell <- approxfun(Variables_df$TIME, Variables_df$CLcell, rule = 2)

varValv <- approxfun(Variables_df$TIME, Variables_df$Valv, rule = 2)
varVlun <- approxfun(Variables_df$TIME, Variables_df$Vlun, rule = 2)
varVR <- approxfun(Variables_df$TIME, Variables_df$VR, rule = 2)

varOraldose_PFOA <- approxfun(Variables_df$TIME, Variables_df$Oraldose_PFOA, rule = 2)
varOraldose_PFOS <- approxfun(Variables_df$TIME, Variables_df$Oraldose_PFOS, rule = 2)

varDrinkdose_PFOA <- approxfun(Variables_df$TIME, Variables_df$Drinkdose_PFOA, rule = 2)
varDrinkdose_PFOS <- approxfun(Variables_df$TIME, Variables_df$Drinkdose_PFOS, rule = 2)

varInhalation_PFOA <- approxfun(Variables_df$TIME, Variables_df$Inhalation_PFOA, rule = 2)
varInhalation_PFOS <- approxfun(Variables_df$TIME, Variables_df$Inhalation_PFOS, rule = 2)

varCsurf_PFOA <- approxfun(Variables_df$TIME, Variables_df$Csurf_PFOA, rule = 2)
varCsurf_PFOS <- approxfun(Variables_df$TIME, Variables_df$Csurf_PFOS, rule = 2)

varTm_PFOA <- approxfun(Variables_df$TIME, Variables_df$Tm_PFOA, rule = 2)
varTm_PFOS <- approxfun(Variables_df$TIME, Variables_df$Tm_PFOS, rule = 2)

#### MODELS ----
# manual adaptation of Westerhout 2024 model, removing cholesterol
#### PFOA and PFOS PBK model ----
PFAS_extended <- function(t, A, parms) {
  with(as.list(c(A, parms)), {
    ## PFOA kinetics --
    kurine_PFOA <- varkurine_PFOA(t)
    kurine_PFOS <- varkurine_PFOS(t)
    QCP <- varQCP(t)
    QG <- varQG(t)
    QL <- varQL(t)
    QF <- varQF(t)
    QK <- varQK(t)
    Qfil <- varQfil(t)
    QSk <- varQSk(t)
    QR <- varQR(t)
    Qp <- varQp(t)
    #VPlas <- varVPlas(t)
    Vart_plas <- varVart_plas(t)
    Vven_plas <- varVven_plas(t)
    VG <- varVG(t)
    VL <- varVL(t)
    VF <- varVF(t)
    VK <- varVK(t)
    Vfil <- varVfil(t)
    VSk <- varVSk(t)
    # Vsurf <- varVsurf(t)
    Vsc <- varVsc(t)
    Vve <- varVve(t)
    CLsc_PFOA <- varCLsc_PFOA(t)
    CLsc_PFOS <- varCLsc_PFOS(t)
    CLve_PFOA <- varCLve_PFOA(t)
    CLve_PFOS <- varCLve_PFOS(t)
    CLcell <- varCLcell(t)
    Valv <- varValv(t)
    Vlun <- varVlun(t)
    VR <- varVR(t)
    Oraldose_PFOA <- varOraldose_PFOA(t) # Dose expressed in ug/kg/day
    Oraldose_PFOS <- varOraldose_PFOS(t) # Dose expressed in ug/kg/day
    Drinkdose_PFOA <- varDrinkdose_PFOA(t) # Dose expressed in ug/kg/day
    Drinkdose_PFOS <- varDrinkdose_PFOS(t) # Dose expressed in ug/kg/day
    Inhalation_PFOA <- varInhalation_PFOA(t) # Dose expressed in ug/kg/day
    Inhalation_PFOS <- varInhalation_PFOS(t) # Dose expressed in ug/kg/day
    Csurf_PFOA <- varCsurf_PFOA(t) # Dose expressed in ug/kg/day
    Csurf_PFOS <- varCsurf_PFOS(t) # Dose expressed in ug/kg/day
    Tm_PFOA <- varTm_PFOA(t)
    Tm_PFOS <- varTm_PFOS(t)
    
    # Concentrations
    # PFOA
    #CAFree_PFOA <- APlas_PFOA/VPlas # free concentration of chemical in plasma µg/L (ng/mL) 
    #CA_PFOA <- CAFree_PFOA/Free_PFOA # total concentration of chemical in plasma 
    CG_PFOA <- AG_PFOA/VG # Concentration in gut (µg/L) 
    CVG_PFOA <- CG_PFOA/PG_PFOA # Concentration leaving gut (µg/L) 
    CL_PFOA = AL_PFOA/VL # Concentration in liver (µg/L) 
    CVL_PFOA = CL_PFOA/PL_PFOA # Concentration leaving liver (µg/L) 
    CF_PFOA = AF_PFOA/VF # Concentration in fat (µg/L)
    CVF_PFOA = CF_PFOA/PF_PFOA # Concentration leaving fat (µg/L) 
    CK_PFOA = AK_PFOA/VK # Concentration in kidneys (µg/L) 
    CVK_PFOA = CK_PFOA/PK_PFOA # Concentration leaving kidneys (µg/L) 
    Cfil_PFOA = Afil_PFOA/Vfil # Concentration in filtrate compartment (µg/L) 
    CSk_PFOA = ASk_PFOA/VSk # Concentration in skin compartment (µg/L) 
    CVSk_PFOA = CSk_PFOA/PSk_PFOA # Concentration leaving skin compartment (µg/L)  
    CR_PFOA = AR_PFOA/VR # Concentration in rest of the body (µg/L) 
    CVR_PFOA = CR_PFOA/PR_PFOA # Concentration leaving rest of the body (µg/L) 
    Clun_bl_PFOA <- Alun_PFOA/(Valv*Pab_PFOA + Vlun)
    Calv_PFOA <- (Alun_PFOA/(Valv*Pab_PFOA + Vlun))*Pab_PFOA
    # Csurf_PFOA <- Asurf_PFOA/Vsurf
    Csc_PFOA <- Asc_PFOA/Vsc
    Ctrans_PFOA <- Atrans_PFOA/0.0000001
    Cve_PFOA <- Ave_PFOA/Vve
    # Original EFSA 2020 code:
    # CvenFree_PFOA <- Aven_PFOA/Vven_plas
    # Cven_PFOA <- CvenFree_PFOA/Free_PFOA
    # CartFree_PFOA <- Aart_PFOA/Vart_plas
    # Cart_PFOA <- CartFree_PFOA/Free_PFOA
    # Correct mass balance:
    Cven_PFOA <- Aven_PFOA/Vven_plas
    CvenFree_PFOA <- Cven_PFOA*Free_PFOA
    Cart_PFOA <- Aart_PFOA/Vart_plas
    CartFree_PFOA <- Cart_PFOA*Free_PFOA
    
    # PFOS
    #CAFree_PFOS <- APlas_PFOS/VPlas # free concentration of chemical in plasma µg/L (ng/mL) 
    #CA_PFOS <- CAFree_PFOS/Free_PFOS # total concentration of chemical in plasma 
    CG_PFOS <- AG_PFOS/VG # Concentration in gut (µg/L) 
    CVG_PFOS <- CG_PFOS/PG_PFOS # Concentration leaving gut (µg/L) 
    CL_PFOS = AL_PFOS/VL # Concentration in liver (µg/L) 
    CVL_PFOS = CL_PFOS/PL_PFOS # Concentration leaving liver (µg/L) 
    CF_PFOS = AF_PFOS/VF # Concentration in fat (µg/L)
    CVF_PFOS = CF_PFOS/PF_PFOS # Concentration leaving fat (µg/L) 
    CK_PFOS = AK_PFOS/VK # Concentration in kidneys (µg/L) 
    CVK_PFOS = CK_PFOS/PK_PFOS # Concentration leaving kidneys (µg/L) 
    Cfil_PFOS = Afil_PFOS/Vfil # Concentration in filtrate compartment (µg/L) 
    CSk_PFOS = ASk_PFOS/VSk # Concentration in skin compartment (µg/L) 
    CVSk_PFOS = CSk_PFOS/PSk_PFOS # Concentration leaving skin compartment (µg/L)  
    CR_PFOS = AR_PFOS/VR # Concentration in rest of the body (µg/L) 
    CVR_PFOS = CR_PFOS/PR_PFOS # Concentration leaving rest of the body (µg/L) 
    Clun_bl_PFOS <- Alun_PFOS/(Valv*Pab_PFOS + Vlun)
    Calv_PFOS <- (Alun_PFOS/(Valv*Pab_PFOS + Vlun))*Pab_PFOS
    # Csurf_PFOS <- Asurf_PFOS/Vsurf
    Csc_PFOS <- Asc_PFOS/Vsc
    Ctrans_PFOS <- Atrans_PFOS/0.0000001
    Cve_PFOS <- Ave_PFOS/Vve
    # Original EFSA 2020 code:
    # CvenFree_PFOS <- Aven_PFOS/Vven_plas
    # Cven_PFOS <- CvenFree_PFOS/Free_PFOS
    # CartFree_PFOS <- Aart_PFOS/Vart_plas
    # Cart_PFOS <- CartFree_PFOS/Free_PFOS
    # Correct mass balance:
    Cven_PFOS <- Aven_PFOS/Vven_plas
    CvenFree_PFOS <- Cven_PFOS*Free_PFOS
    Cart_PFOS <- Aart_PFOS/Vart_plas
    CartFree_PFOS <- Cart_PFOS*Free_PFOS
    
    # Lung compartment
    dAlun_PFOA <- QCP*Cven_PFOA*Free_PFOA + Qp*Inhalation_PFOA - QCP*Clun_bl_PFOA*Free_PFOA #- Qp*Calv_PFOA
    
    dAlun_PFOS <- QCP*Cven_PFOS*Free_PFOS + Qp*Inhalation_PFOS - QCP*Clun_bl_PFOS*Free_PFOS #- Qp*Calv_PFOS
    
    # Skin compartment 
    # Skin surface
    # dAsurf_PFOA <- - 2*CLsc_PFOA*Csurf_PFOA + 2*CLsc_PFOA*Csc_PFOA - Remove*Asurf_PFOA
    # 
    # dAsurf_PFOS <- - 2*CLsc_PFOS*Csurf_PFOS + 2*CLsc_PFOS*Csc_PFOS - Remove*Asurf_PFOS
    
    # Stratum corneum
    dAsc_PFOA <- 2*CLsc_PFOA*Csurf_PFOA - 4*CLsc_PFOA*Csc_PFOA + 2*CLsc_PFOA*Ctrans_PFOA
    
    dAsc_PFOS <- 2*CLsc_PFOS*Csurf_PFOS - 4*CLsc_PFOS*Csc_PFOS + 2*CLsc_PFOS*Ctrans_PFOS
    
    # Transfer compartment
    dAtrans_PFOA <- 2*CLsc_PFOA*Csc_PFOA - 2*CLsc_PFOA*Ctrans_PFOA - 2*CLve_PFOA*(Ctrans_PFOA/Kscve_PFOA) + 2*CLve_PFOA*Cve_PFOA
    
    dAtrans_PFOS <- 2*CLsc_PFOS*Csc_PFOS - 2*CLsc_PFOS*Ctrans_PFOS - 2*CLve_PFOS*(Ctrans_PFOS/Kscve_PFOS) + 2*CLve_PFOS*Cve_PFOS
    
    # Viable epidermis
    dAve_PFOA <- 2*CLve_PFOA*(Ctrans_PFOA/Kscve_PFOA) - 2*CLve_PFOA*Cve_PFOA - CLcell*Cve_PFOA/Kver_PFOA
    
    dAve_PFOS <- 2*CLve_PFOS*(Ctrans_PFOS/Kscve_PFOS) - 2*CLve_PFOS*Cve_PFOS - CLcell*Cve_PFOS/Kver_PFOS
    
    # Skin compartment
    dASk_PFOA <- CLcell*Cve_PFOA/Kver_PFOA + QSk*(Cart_PFOA*Free_PFOA-CSk_PFOA*FreeSk_PFOA)    # Rate of change in skin(µg/h) 
    
    dASk_PFOS <- CLcell*Cve_PFOS/Kver_PFOS + QSk*(Cart_PFOS*Free_PFOS-CSk_PFOS*FreeSk_PFOS)    # Rate of change in skin(µg/h) 
    
    # Plasma compartment      
    # dAPlas_PFOA <- QF*CF_PFOA*FreeF_PFOA + (QL+QG)*CL_PFOA*FreeL_PFOA + QR*CR_PFOA*FreeR_PFOA + QSk*CSk_PFOA*FreeSk_PFOA + 
    #   QK*CK_PFOA*FreeK_PFOA - QCP*CA_PFOA*Free_PFOA - Qfil*CA_PFOA*Free_PFOA 
    # 
    # dAPlas_PFOS <- QF*CF_PFOS*FreeF_PFOS + (QL+QG)*CL_PFOS*FreeL_PFOS + QR*CR_PFOS*FreeR_PFOS + QSk*CSk_PFOS*FreeSk_PFOS + 
    #   QK*CK_PFOS*FreeK_PFOS - QCP*CA_PFOS*Free_PFOS - Qfil*CA_PFOS*Free_PFOS 
    # 
    # Venous blood (plasma) compartment      
    dAven_PFOA <- QF*CF_PFOA*FreeF_PFOA + (QL+QG)*CL_PFOA*FreeL_PFOA + QR*CR_PFOA*FreeR_PFOA + QSk*CSk_PFOA*FreeSk_PFOA + 
      QK*CK_PFOA*FreeK_PFOA - QCP*Cven_PFOA*Free_PFOA
    
    dAven_PFOS <- QF*CF_PFOS*FreeF_PFOS + (QL+QG)*CL_PFOS*FreeL_PFOS + QR*CR_PFOS*FreeR_PFOS + QSk*CSk_PFOS*FreeSk_PFOS + 
      QK*CK_PFOS*FreeK_PFOS - QCP*Cven_PFOS*Free_PFOS
    
    # Arterial blood (plasma) compartment      
    dAart_PFOA <- QCP*Clun_bl_PFOA*Free_PFOA - QCP*Cart_PFOA*Free_PFOA - Qfil*Cart_PFOA*Free_PFOA 
    
    dAart_PFOS <- QCP*Clun_bl_PFOS*Free_PFOS - QCP*Cart_PFOS*Free_PFOS - Qfil*Cart_PFOS*Free_PFOS 
    
    # Gut compartment 
    dAG_PFOA <- QG*(Cart_PFOA*Free_PFOA-CG_PFOA*FreeG_PFOA) + Oraldose_PFOA + Drinkdose_PFOA 
    
    dAG_PFOS <- QG*(Cart_PFOS*Free_PFOS-CG_PFOS*FreeG_PFOS) + Oraldose_PFOS + Drinkdose_PFOS 
    
    # Liver compartment 
    dAL_PFOA <- (QL*(Cart_PFOA*Free_PFOA)) + (QG*CG_PFOA*FreeG_PFOA) - ((QL+QG)*CL_PFOA*FreeL_PFOA) # Rate of change in liver (ug/h) 
    
    dAL_PFOS <- (QL*(Cart_PFOS*Free_PFOS)) + (QG*CG_PFOS*FreeG_PFOS) - ((QL+QG)*CL_PFOS*FreeL_PFOS) # Rate of change in liver (ug/h) 
    
    # Fat compartment 
    dAF_PFOA <- QF*(Cart_PFOA*Free_PFOA-CF_PFOA*FreeF_PFOA)   # Rate of change in fat (µg/h) 
    
    dAF_PFOS <- QF*(Cart_PFOS*Free_PFOS-CF_PFOS*FreeF_PFOS)   # Rate of change in fat (µg/h) 
    
    # Kidney compartment 
    dAK_PFOA <- QK*(Cart_PFOA*Free_PFOA-CK_PFOA*FreeK_PFOA) + (Tm_PFOA*Cfil_PFOA)/(Kt_PFOA+Cfil_PFOA) # Rate of change in kidneys (µg/h) 
    
    dAK_PFOS <- QK*(Cart_PFOS*Free_PFOS-CK_PFOS*FreeK_PFOS) + (Tm_PFOS*Cfil_PFOS)/(Kt_PFOS+Cfil_PFOS) # Rate of change in kidneys (µg/h) 
    
    # Filtrate compartment 
    dAfil_PFOA <- Qfil*(Cart_PFOA*Free_PFOA-Cfil_PFOA) - (Tm_PFOA*Cfil_PFOA)/(Kt_PFOA+Cfil_PFOA) # Rate of change in filtrate compartment (ug/h) 
    
    dAfil_PFOS <- Qfil*(Cart_PFOS*Free_PFOS-Cfil_PFOS) - (Tm_PFOS*Cfil_PFOS)/(Kt_PFOS+Cfil_PFOS) # Rate of change in filtrate compartment (ug/h) 
    
    # Storage compartment for urine 
    dAdelay_PFOA <- Qfil*Cfil_PFOA - kurine_PFOA*Adelay_PFOA   
    
    dAdelay_PFOS <- Qfil*Cfil_PFOS - kurine_PFOS*Adelay_PFOS   
    
    # Urine 
    dAurine_PFOA <- kurine_PFOA*Adelay_PFOA 
    
    dAurine_PFOS <- kurine_PFOS*Adelay_PFOS 
    
    # Rest of the body 
    dAR_PFOA <- QR*(Cart_PFOA*Free_PFOA-CR_PFOA*FreeR_PFOA)   # Rate of change in rest of the body (µg/h) 
    
    dAR_PFOS <- QR*(Cart_PFOS*Free_PFOS-CR_PFOS*FreeR_PFOS)   # Rate of change in rest of the body (µg/h) 
    
    list(c(dAlun_PFOA, dAsc_PFOA, dAtrans_PFOA, dAve_PFOA, dASk_PFOA, dAven_PFOA, dAart_PFOA,
           dAG_PFOA, dAL_PFOA, dAF_PFOA, dAK_PFOA, dAfil_PFOA, dAdelay_PFOA, dAurine_PFOA, dAR_PFOA,
           dAlun_PFOS, dAsc_PFOS, dAtrans_PFOS, dAve_PFOS, dASk_PFOS, dAven_PFOS, dAart_PFOS,
           dAG_PFOS, dAL_PFOS, dAF_PFOS, dAK_PFOS, dAfil_PFOS, dAdelay_PFOS, dAurine_PFOS, dAR_PFOS),
         c(Cart_PFOA=Cart_PFOA, Cven_PFOA=Cven_PFOA, Clun_bl_PFOA=Clun_bl_PFOA,
           Csurf_PFOA=Csurf_PFOA,Csc_PFOA=Csc_PFOA,Ctrans_PFOA=Ctrans_PFOA,Cve_PFOA=Cve_PFOA, CSk_PFOA=CSk_PFOA,
           CG_PFOA=CG_PFOA, CVG_PFOA=CVG_PFOA, CL_PFOA=CL_PFOA, CVL_PFOA=CVL_PFOA, CF_PFOA=CF_PFOA,CVF_PFOA=CVF_PFOA,
           CK_PFOA=CK_PFOA, CVK_PFOA=CVK_PFOA, Cfil_PFOA=Cfil_PFOA, CR_PFOA=CR_PFOA, CVR_PFOA=CVR_PFOA,
           Cart_PFOS=Cart_PFOS, Cven_PFOS=Cven_PFOS, Clun_bl_PFOS=Clun_bl_PFOS,
           Csurf_PFOS=Csurf_PFOS,Csc_PFOS=Csc_PFOS,Ctrans_PFOS=Ctrans_PFOS,Cve_PFOS=Cve_PFOS, CSk_PFOS=CSk_PFOS,
           CG_PFOS=CG_PFOS, CVG_PFOS=CVG_PFOS, CL_PFOS=CL_PFOS, CVL_PFOS=CVL_PFOS, CF_PFOS=CF_PFOS,CVF_PFOS=CVF_PFOS,
           CK_PFOS=CK_PFOS, CVK_PFOS=CVK_PFOS, Cfil_PFOS=Cfil_PFOS, CR_PFOS=CR_PFOS, CVR_PFOS=CVR_PFOS))
  })
}

#### parms ----
# Parameters used in the model
parms_PFAS_extended <- c(Kt_PFOA, Free_PFOA, FreeL_PFOA, FreeF_PFOA, FreeK_PFOA, FreeSk_PFOA, FreeR_PFOA, FreeG_PFOA,
                                     PL_PFOA, PF_PFOA, PK_PFOA, PSk_PFOA, PR_PFOA, PG_PFOA,
                                     Kt_PFOS, Free_PFOS, FreeL_PFOS, FreeF_PFOS, FreeK_PFOS, FreeSk_PFOS, FreeR_PFOS, FreeG_PFOS,
                                     PL_PFOS, PF_PFOS, PK_PFOS, PSk_PFOS, PR_PFOS, PG_PFOS,
                                     Pab_PFOA,
                                     Pab_PFOS,
                                     Kscve_PFOA,
                                     Kscve_PFOS,
                                     Kver_PFOA,
                                     Kver_PFOS,
                                     Remove = 0)

#### A_init ----
# Initial values
A_init_PFAS_extended <- c(Alun_PFOA = Alun_PFOA_background,
                                      Asc_PFOA = Asc_PFOA_background,
                                      Atrans_PFOA = Atrans_PFOA_background,
                                      Ave_PFOA = Ave_PFOA_background,
                                      ASk_PFOA = ASk_PFOA_background,
                                      Aven_PFOA = Aven_PFOA_background,
                                      Aart_PFOA = Aart_PFOA_background,
                                      AG_PFOA = AG_PFOA_background,
                                      AL_PFOA = AL_PFOA_background, 
                                      AF_PFOA = AF_PFOA_background,
                                      AK_PFOA = AK_PFOA_background,
                                      Afil_PFOA = Afil_PFOA_background,
                                      Adelay_PFOA = Adelay_PFOA_background,
                                      Aurine_PFOA = Aurine_PFOA_background,
                                      AR_PFOA = AR_PFOA_background,
                                      Alun_PFOS = Alun_PFOS_background,
                                      Asc_PFOS = Asc_PFOS_background,
                                      Atrans_PFOS = Atrans_PFOS_background,
                                      Ave_PFOS = Ave_PFOS_background,
                                      ASk_PFOS = ASk_PFOS_background,
                                      Aven_PFOS = Aven_PFOS_background,
                                      Aart_PFOS = Aart_PFOS_background,
                                      AG_PFOS = AG_PFOS_background,
                                      AL_PFOS = AL_PFOS_background, 
                                      AF_PFOS = AF_PFOS_background,
                                      AK_PFOS = AK_PFOS_background,
                                      Afil_PFOS = Afil_PFOS_background,
                                      Adelay_PFOS = Adelay_PFOS_background,
                                      Aurine_PFOS = Aurine_PFOS_background,
                                      AR_PFOS = AR_PFOS_background)

#### OUTPUT ----
# PFOA & PFOS
output_PFAS <- lsoda(A_init_PFAS_extended, TIME, PFAS_extended, parms_PFAS_extended)
output.df <- as.data.frame(output_PFAS)
output.df <- Variables_df %>%
  rename(time = TIME) %>%
  left_join(output.df)

# Depending on your settings in Excel, you could either save the comma separated file with a ',' or '.' as decimal separator.
# write.table(output.df, file = "Westerhout 2024 - correct mass balance results 20240618.csv", dec = ",", sep = ";", row.names = FALSE, na = "")
write.table(output.df, file = "Westerhout 2024 - correct mass balance results 20240618.csv", dec = ".", sep = ";", row.names = FALSE, na = "")

# Results_Amounts <- output.df %>%
#   select("year",starts_with("A"))
# 
# write.table(Results_Amounts, file = "Background amounts 20240306.csv", dec = ",", sep = ";", row.names = FALSE, na = "")

#### PLOTS ----
# # PFOA concentration in plasma and liver (ug/L)
# ggplot() +
#   geom_line(data=output.df, aes(x=age, y=Cart_PFOA, color="Cart_PFOA", lty="Cart_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=CL_PFOA, color="CL_PFOA", lty="CL_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=Cart_PFOS, color="Cart_PFOS", lty="Cart_PFOS")) +
#   geom_line(data=output.df, aes(x=age, y=CL_PFOS, color="CL_PFOS", lty="CL_PFOS")) +
#   
#   scale_colour_manual(name='Output',
#                       values=c('Cart_PFOA'='#FF0000',
#                                'CL_PFOA'='#0000FF',
#                                'Cart_PFOS'='#FF0000',
#                                'CL_PFOS'='#0000FF'),
#                       labels=c('Cart_PFOA'='PFOA concentration in plasma',
#                                'CL_PFOA'='PFOA concentration in liver',
#                                'Cart_PFOS'='PFOS concentration in plasma',
#                                'CL_PFOS'='PFOS concentration in liver')) +
#   
#   scale_linetype_manual(name='Output',
#                       values=c('Cart_PFOA'='solid',
#                                'CL_PFOA'='solid',
#                                'Cart_PFOS'='dashed',
#                                'CL_PFOS'='dashed'),
#                       labels=c('Cart_PFOA'='PFOA concentration in plasma',
#                                'CL_PFOA'='PFOA concentration in liver',
#                                'Cart_PFOS'='PFOS concentration in plasma',
#                                'CL_PFOS'='PFOS concentration in liver')) +
# 
#   theme_bw() +
#   labs(title="PFAS concentration-time profiles",x="\nAge (y)", y="Concentration (\u03BCg/L)\n") +
#   theme(plot.title = element_text(hjust = 0.5))

# PFOA concentration in plasma (ug/L)
Figure_PFAS <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Cart_PFOA, color="Cart_PFOA", lty="Cart_PFOA")) +
  geom_line(data=output.df, aes(x=age, y=Cart_PFOS, color="Cart_PFOS", lty="Cart_PFOS")) +

  scale_colour_manual(name='Output',
                      values=c('Cart_PFOA'='#000000',
                               'Cart_PFOS'='#000000'),
                      labels=c('Cart_PFOA'='PFOA concentration in plasma',
                               'Cart_PFOS'='PFOS concentration in plasma')) +

  scale_linetype_manual(name='Output',
                        values=c('Cart_PFOA'='solid',
                                 'Cart_PFOS'='dashed'),
                        labels=c('Cart_PFOA'='PFOA concentration in plasma',
                                 'Cart_PFOS'='PFOS concentration in plasma')) +

  theme_bw() +
  labs(title="PFAS concentration-time profiles",x="\nAge (y)", y="Concentration (\u03BCg/L)\n") +
  theme(plot.title = element_text(hjust = 0.5))

# tiff("Westerhout 2024 - Figure PFAS kinetics.tif",
#      res=600, compression = "lzw",
#      height=6,
#      width=8,
#      units="in")
# Figure_PFAS
# 
# dev.off()

Figure_PFAS

# PFOA concentration in liver, intestine and periphery (ug/L)
#### Figure PFAS ----
# Figure_PFAS <- ggplot() +
#   geom_line(data=output.df, aes(x=age, y=Cart_PFOA, color="Cart_PFOA", lty="Cart_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=CL_PFOA, color="CL_PFOA", lty="CL_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=CG_PFOA, color="CG_PFOA", lty="CG_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=CR_PFOA, color="CR_PFOA", lty="CR_PFOA")) +
#   geom_line(data=output.df, aes(x=age, y=Cart_PFOS, color="Cart_PFOS", lty="Cart_PFOS")) +
#   geom_line(data=output.df, aes(x=age, y=CL_PFOS, color="CL_PFOS", lty="CL_PFOS")) +
#   geom_line(data=output.df, aes(x=age, y=CG_PFOS, color="CG_PFOS", lty="CG_PFOS")) +
#   geom_line(data=output.df, aes(x=age, y=CR_PFOS, color="CR_PFOS", lty="CR_PFOS")) +
#   scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80),labels=c(0,10,20,30,40,50,60,70,80)) + 
#   scale_y_continuous(breaks=seq(0, 25, by = 5), limits=c(0,25)) +
#   coord_cartesian(xlim=c(0,80)) +
#   
#   scale_colour_manual(name='Output',
#                       values=c('Cart_PFOA'='#000000',
#                                'CL_PFOA'='#FF0000',
#                                'CG_PFOA'='#FFB700',
#                                'CR_PFOA'='#0022FF',
#                                'Cart_PFOS'='#000000',
#                                'CL_PFOS'='#FF0000',
#                                'CG_PFOS'='#FFB700',
#                                'CR_PFOS'='#0022FF'),
#                       labels=c('Cart_PFOA'='PFOA concentration in plasma',
#                                'CL_PFOA'='PFOA concentration in liver',
#                                'CG_PFOA'='PFOA concentration in intestine',
#                                'CR_PFOA'='PFOA concentration in periphery',
#                                'Cart_PFOS'='PFOS concentration in plasma',
#                                'CL_PFOS'='PFOS concentration in liver',
#                                'CG_PFOS'='PFOS concentration in intestine',
#                                'CR_PFOS'='PFOS concentration in periphery')) +
#   
#   scale_linetype_manual(name='Output',
#                         values=c('Cart_PFOA'='solid',
#                                  'CL_PFOA'='solid',
#                                  'CG_PFOA'='solid',
#                                  'CR_PFOA'='solid',
#                                  'Cart_PFOS'='dashed',
#                                  'CL_PFOS'='dashed',
#                                  'CG_PFOS'='dashed',
#                                  'CR_PFOS'='dashed'),
#                         labels=c('Cart_PFOA'='PFOA concentration in plasma',
#                                  'CL_PFOA'='PFOA concentration in liver',
#                                  'CG_PFOA'='PFOA concentration in intestine',
#                                  'CR_PFOA'='PFOA concentration in periphery',
#                                  'Cart_PFOS'='PFOS concentration in plasma',
#                                  'CL_PFOS'='PFOS concentration in liver',
#                                  'CG_PFOS'='PFOS concentration in intestine',
#                                  'CR_PFOS'='PFOS concentration in periphery')) +
#   
#   theme_bw() +
#   labs(title="PFAS concentration-time profiles",x="\nAge (y)", y="Concentration (\u03BCg/L)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_PFAS

# tiff("Figure PFAS kinetics - EFSA TWI.tif",
#      res=600, compression = "lzw",
#      height=6,
#      width=8,
#      units="in")
# Figure_PFAS
# 
# dev.off()
