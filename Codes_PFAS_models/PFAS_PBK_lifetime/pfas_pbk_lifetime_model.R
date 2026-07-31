#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# This function simulates PFAS concentrations over time using age-dependent physiology, oral exposure route, tissue distribution, renal reabsorption, and elimination. 
# The ODE system is solved with deSolve::lsoda().
#
# MAIN INPUT GROUPS
# -----------------
# Individual: age, BW, height, Gender and BWbirth.
# Simulation: TSTART, TSTOP and DT [days].
# Physiology: organ-size variability factors, heart rate, haematocrit, plasma volumes, blood-flow parameters and skin parameters.
# Exposure: Oralexpo, Drinkconc.
# Kinetics: Tmc, Kt, kurinec and CLfaeces.
# Distribution: Free and tissue partition coefficients (PL, PF, PK, PSk, PR, PG and PLun).
# Maternal/early life: maternal concentration, placental and milk-transfer factors, milk decline and initial compartment fractions.
#
# MODEL COMPARTMENTS
# ------------------
# Lung; stratum corneum; dermal transfer layer; viable epidermis; perfused skin; venous plasma; arterial plasma; gut; liver; fat; kidney; kidney filtrate;
# urinary delay compartment; urine; remaining body; menstrual loss; faeces.
#
# MAIN TIME-VARYING VARIABLES
# ---------------------------
# BW, height, BMI and BSA describe body growth.
# V* variables are compartment volumes; Q* variables are plasma or ventilation flows; 
# CL* variables are clearance terms.
# Cart and Cven are arterial and venous plasma concentrations.
# Atot is the total amount in the model, including cumulative elimination.
#
# EXPOSURE AND ELIMINATION
# ------------------------
# Exposure occurs through breast milk, oral intake, drinking water.
# Elimination occurs through urine, faeces and menstrual blood loss.
#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pbk.model<-function(age                   =NA, # age (in years) at day of plasma sampling
                    BW                    =NA, # body weight (kg)
                    height                =NA, # height (cm)     
                    Gender                =NA, # "F" or "M"
                    TSTART                =NA, # start of modelled period (days)
                    TSTOP                 =NA, # end of modelled period (days)
                    DT                    =NA, # modelling time step (days)
                    exposure_duration     =NA, # duration of exposure (d)
                    age_start_occ_exp     =NA, # occupational exposure starting age (dermal and inhalation)
                    age_stop_occ_exp      =NA, # occupational exposure stopping age (dermal and inhalation)
                    age_start_menstruation=NA, # years
                    age_stop_menstruation =NA, # years
                    BWbirth               =NA, # body weight at birth (kg)
                    MFV                   =NA, # menstrual fluid volume (L/d) recalculated from 28-day value
                    HR                    =NA, # heart beat rate (1/min)
                    SDLun                 =NA, # added standard deviation of the volume of lungs from Bosgra et al. (2012)
                    SDK                   =NA, # added standard deviation of the volume of kidney from Bosgra et al. (2012)
                    SDInt                 =NA, # added standard deviation of the volume of intestine from Bosgra et al. (2012)
                    SDL                   =NA, # added standard deviation of the volume of liver from Bosgra et al. (2012)
                    SDSk                  =NA, # added standard deviation of the volume of skin from Bosgra et al. (2012)
                    SDAdi                 =NA, # added standard deviation of the volume of lungs from Bosgra et al. (2012)
                    SDB                   =NA, # added standard deviation of the volume of blood from Bosgra et al. (2012)
                    QfilC                 =NA, # fraction of kidney plasma flow to filtrate
                    #fractional tissue volumes
                    Htc                   =NA, # hematocrit
                    VplasC                =NA, # plasma volume
                    VartC                 =NA, # arterial blood volume
                    VvenC                 =NA, # venous blood volume
                    #skin parameters
                    SkinThickness         =NA, # skin thickness (cm)
                    fss                   =NA, # fractional skin storage
                    hsc                   =NA, # thickness of stratum corneum (cm)
                    hve                   =NA, # thickness of viable epidermis (cm)
                    # fat‐fractions for Kscve/Kver:
                    ffatsc                =NA, # fraction of fat in stratum corneum
                    ffatve                =NA, # fraction of fat in viable epidermis
                    ffatepi               =NA, # fraction of fat in epidermis
                    ffatbl                =NA, # fraction of fat in blood
                    Kpcell                =NA, # cell membrane partition kinetic constant (cm/h)
                    Drinkrate             =NA, # drink rate (mL/kg/d)
                    MW                    =NA, # molecular weight (g/mol)
                    logP                  =NA, # logP
                    VP                    =NA, # vapour pressure
                    Kpve                  =NA, # partition koefficient for viable epidermis
                    Kpsc                  =NA, # partition coefficient for stratum corneum
                    # exposures
                    Oralexpo              =NA, # oral exposure (ng/kg/d)
                    Drinkconc             =NA, # drink concentration
                    Cinh                  =NA, # air concentration
                    Cdermal               =NA, # for this model, dermal & inhalation exposures are excluded
                    # kinetics
                    Tmc                   =NA, # maximum resorption rate
                    Kt                    =NA, # tubular transport coefficient
                    kurinec               =NA, # urine clearance coefficient
                    CLfaeces              =NA, # fecal clearance
                    # free & partitions
                    Free                  =NA, # free fraction of plasma
                    PL                    =NA, # partition coefficient of liver
                    PF                    =NA, # partition coefficient of adipose tissue
                    PK                    =NA, # partition coefficient of kidney
                    PSk                   =NA, # partition coefficient of skin
                    PR                    =NA, # partition coefficient of rest
                    PG                    =NA, # partition coefficient of intestine
                    PLun                  =NA, # partition coefficient of lungs
                    # maternal/breast milk
                    maternal              =NA, # maternal milk concentration
                    PT                    =NA, # placental transfer factor
                    Ratio                 =NA, # milk concentration/maternal serum concentration during breastfeeding
                    DECLINE               =NA, # milk decline
                    Milkconsumption       =NA, # maternal milk consumption (L/d)
                    # birth fractions
                    Aartbirth             =NA, # amount of PFAS in arterial blood at birth
                    Avenbirth             =NA, # amount of PFAS in venous blood at birth
                    AGbirth               =NA, # amount of PFAS in intestine at birth
                    ALbirth               =NA, # amount of PFAS in liver at birth
                    AFbirth               =NA, # amount of PFAS in adipose tissue at birth
                    AKbirth               =NA, # amount of PFAS in kidney at birth
                    ALunbirth             =NA, # amount of PFAS in lungs at birth
                    ASkbirth              =NA, # amount of PFAS in skin at birth
                    ARbirth               =NA, # amount of PFAS in rest of the body at birth
                    Afilbirth             =NA, # amount of PFAS in filtrate at birth
                    Adelaybirth           =NA, # amount of PFAS in  at birth
                    Aurinebirth           =NA  # amount of PFAS in urine at birth
) {
################################################################################################## 
#                                                                                                #
#  Group model inputs into settings, physiological parameters, and chemical-specific parameters  #
#                                                                                                #
################################################################################################## 
  
  settings<-list(
    "TSTART"=TSTART,
    "TSTOP"=TSTOP, # days
    "DT"=DT,
    "exposure_duration"=exposure_duration, # Duration of exposure (d)
    "age_start_occ_exp"=age_start_occ_exp, # occupational exposure starting age (dermal and inhalation)
    "age_stop_occ_exp"=age_stop_occ_exp, # occupational exposure stopping age (dermal and inhalation)
    "age_start_menstruation"=age_start_menstruation,
    "age_stop_menstruation"=age_stop_menstruation)
  
  physio <- list(
    #fractional blood flows
    "BWbirth"=BWbirth,     # kg
    "MFV"=MFV,  # L/28 days
    #"QCC"=QCC,      # L/d/kg^0.75
    #"QFC"=QFC,
    #"QLC"=QLC,
    #"QKC"=QKC,
    #"QSkC"=QSkC,     
    #"QGC"=QGC,       
    #fractional tissue volumes
    "Htc"=Htc,
    "VplasC"=VplasC,
    "VartC"=VartC,
    "VvenC"=VvenC,
    #skin parameters
    "SkinThickness"=SkinThickness,
    "fss"=fss,
    "hsc"=hsc,
    "hve"=hve,
    # fat‐fractions for Kscve/Kver:
    "ffatsc"=ffatsc,
    "ffatve"=ffatve,
    "ffatepi"=ffatepi,
    "ffatbl"=ffatbl,
    "Kpcell"=Kpcell,     # cm/h
    "Drinkrate"=Drinkrate)        # mL/kg/d
  
  
  p<-list(
    # physicochemical property
    "MW"=MW,
    "logP"=logP,
    "VP"=VP,
    "Kpsc"=Kpsc,
    "Kpve"=Kpve,
    # exposures
    "exposure"=list(
      "Oralexpo"=Oralexpo,
      "Drinkconc"=Drinkconc,
      "Cinh"=Cinh,
      "Cdermal"=Cdermal), #for this model, dermal & inhalation exposures are excluded
    
    # kinetics
    "kinetics" = list(
      "Tmc"=Tmc,
      "Kt"=Kt,
      "kurinec"=kurinec,
      "CLfaeces"=CLfaeces),
    
    # free & partitions
    "Free"=Free,
    "partition"=list(
      "PL"=PL,
      "PF"=PF,
      "PK"=PK,
      "PSk"=PSk,
      "PR"=PR,
      "PG"=PG,
      "PLun"=PLun),
    
    # maternal/breast milk
    "maternal_c"=list(
      "maternal"=maternal,
      "PT"=PT,
      "Ratio"=Ratio,
      "DECLINE"=DECLINE,
      "Milkconsumption"=Milkconsumption),
    
    # birth fractions
    "birth"=list(
      "Aartbirth"=Aartbirth,
      "Avenbirth"=Avenbirth,
      "AGbirth"=AGbirth,
      "ALbirth"=ALbirth,
      "AFbirth"=AFbirth,
      "AKbirth"=AKbirth,
      "ALunbirth"=ALunbirth,
      "ASkbirth"=ASkbirth,
      "ARbirth"=ARbirth,
      "Afilbirth"=Afilbirth,
      "Adelaybirth"=Adelaybirth,
      "Aurinebirth"=Aurinebirth))
  
#################################################################################  
#                                                                               #
#  Helper functions used to build time-varying inputs and solve the PBPK model  #
#                                                                               #
#################################################################################  
  
  # 1) build the time‐varying covariates for one individual of any selected PFAS
  #generates a data frame of every time point from TSTART to TSTOP, then computes the content.
  Variables_df <- function(age, BW, height, Gender) {
    TIME <- seq(settings$TSTART, settings$TSTOP, by = settings$DT)
    n    <- length(TIME)
    day  <- TIME+1
    
    BWtime<-(BWbirth+4.47*(TIME/365)-0.093*(TIME/365)^2+0.00061*(TIME/365)^3)*BW/(BWbirth+4.47*(TSTOP/365)-0.093*(TSTOP/365)^2+0.00061*(TSTOP/365)^3)         # change to fit BW
    
    base <- data.frame("TIME"          = TIME,
                       "day"           = day,
                       "dayoftheweek"  = rep(1:7, length.out = n),
                       "dayofthemonth" = rep(1:30, length.out = n),
                       "month"         = ceiling(day/(365/12)) - ((ceiling(day/365)-1)*12),
                       "year"          = ceiling(day/365),
                       "age"           = TIME/365,
                       "Gender"        = Gender,
                       "SkinTarea"     = 9.1 * ((BWtime * 1000)^0.666))                    # skinTarea later replaced by BSA (body-surface area) 
    
    
    base$Vsc<-physio$fss*base$SkinTarea*physio$hsc/1000
    base$Vve<-physio$fss*base$SkinTarea*physio$hve/1000
    
#########################################################################################################
#                                                                                                       #
#  Calculate age-dependent body size, organ volumes, and blood flows according to Bosgra et al. (2012)  #
#                                                                                                       #
#########################################################################################################
    
# Body-weight growth curves; adult points are added to smooth interpolation.
    originalweight<-BW
    age_curve<-c(0,0.004,0.047,((1:11)+0.5)/12,1.125,1.375,1.625,1.875,2.25,2.75,3.25,3.75,(4:18)+0.5,22.5,200)
    w_m<-c(3.4,3.5,3.7,4.8,5.8,6.4,7.1,7.7,8.0,8.7,8.9,9.5,9.6,10.0,10.4,11.3,11.8,12.8,13.5,14.9,15.7,16.7,18.8,20.8,24.2,27.0,30.4,33.6,37.5,41.3,47.0,52.4,58.8,64.2,67.5,70.0,72.2,75.0,75.0)
    w_f<-c(3.2,3.3,3.5,4.5,5.3,5.9,6.5,7.2,7.5,8.1,8.3,8.8,9.0, 9.4, 9.7,10.7,11.2,12.1,13.0,14.5,15.1,16.4,18.3,20.1,23.6,26.3,29.5,32.7,37.3,41.8,47.1,51.3,54.6,56.8,58.1,58.9,59.5,60.0,62.5)
    if(Gender=="M") {
      base$BW<-(approx(age_curve,w_m,base$age)$y-approx(age_curve,w_m,0)$y)/(approx(age_curve,w_m,max(base$day/365.24))$y-approx(age_curve,w_m,0)$y)*(originalweight-BWbirth)+BWbirth
    } else {
      base$BW<-(approx(age_curve,w_f,base$age)$y-approx(age_curve,w_f,0)$y)/(approx(age_curve,w_f,max(base$day/365.24))$y-approx(age_curve,w_f,0)$y)*(originalweight-BWbirth)+BWbirth
    }
    
    # Height growth curves.
    originalheight<-height
    age_curve<-c(0,0.25,0.5*c(1:36),200)
    h_m<-c(50.93,61.72,68.36,76.58,82.91,88.24,92.96,97.20,101.09,104.72,108.13,111.41,114.82,118.07,121.42,124.73,127.93,130.96,133.73,136.40,139.08,141.72,144.16,146.76,149.60,152.51,155.71,159.28,163.20,167.01,170.27,173.12,175.47,176.95,177.97,178.77,179.41,179.81,179.81)/100
    h_f<-c(50.21,60.26,66.57,75.06,81.64,86.84,91.47,95.90, 99.93,103.60,107.22,110.89,114.33,117.64,120.83,123.85,126.84,129.79,132.64,135.48,138.37,141.31,144.36,147.68,150.93,154.26,157.18,159.55,161.47,163.04,164.21,165.07,165.65,166.05,166.21,166.31,166.49,166.65,166.65)/100
    if(Gender=="M") {
      base$height<-originalheight/approx(age_curve,h_m,age)$y*approx(age_curve,h_m,base$age)$y
    } else {
      base$height<-originalheight/approx(age_curve,h_f,age)$y*approx(age_curve,h_f,base$age)$y
    }
    
    # Body mass index and body surface area.
    BMI<-base$BW/base$height^2
    
    BSA<-0.007184*base$BW^0.425*(100*base$height)^0.725
    
    # Organ masses estimated from body size.
    mLun<-exp(2.1 *log(base$height)-2.092 )*SDLun   # lungs
    mK  <-exp(1.93*log(base$height)-2.306 )*SDK     # kidneys
    mInt<-exp(2.47*log(base$height)-1.351 )*SDInt   # intestines
    mL  <-exp(1.98*log(base$height)-0.6786)*SDL     # liver
    mSk <-exp(1.64*BSA             -1.93  )*SDSk    # skin

    # Sex- and age-dependent brain, blood, gonad, and adipose masses.
    a00.01<-which(                       base$TIME<01*365.2425)
    a04.15<-which(base$TIME>=04*365.2425&base$TIME<15*365.2425)
    a18.99<-which(base$TIME>=18*365.2425                      )
    
    if(Gender=="M") {
      mB  <-3.33*BSA-0.81                                             *SDB                                                 # blood
      mAdi00.01<-(0.9084+0.706*base$BW[a00.01]+5.3*base$height[a00.01]-3.057*0.365*base$age[a00.01])/100*base$BW[a00.01]   # adipose between 0 and 1 year of age
      mAdi04.15<-(1.51*BMI[a04.15]-0.70*base$age[a04.15]- 3.6+1.4)/100*base$BW[a04.15]                                     # adipose between 4 and 15 years of age
      mAdi18.99<-(1.20*BMI[a18.99]+0.23*base$age[a18.99]-10.8-5.4)/100*base$BW[a18.99]                                     # adipose between 18 and 200 years of age
    } else {
      mB  <-2.66*BSA-0.46                                             *SDB                                                         # blood
      mAdi00.01<-(0.9084+0.706*base$BW[a00.01]+5.3*base$height[a00.01]+0.3585-3.057*0.365*base$age[a00.01])/100*base$BW[a00.01]    # adipose between 0 and 1 year of age
      mAdi04.15<-(1.51*BMI[a04.15]-0.70*base$age[a04.15]+1.4)/100*base$BW[a04.15]                                                  # adipose between 4 and 15 years of age
      mAdi18.99<-(1.20*BMI[a18.99]+0.23*base$age[a18.99]-5.4)/100*base$BW[a18.99]                                                  # adipose between 18 and 200 years of age
    }
    mAdi<-approx(c(0,a00.01,a04.15,a18.99,200*365.2425),c(mAdi00.01[1],mAdi00.01,mAdi04.15,mAdi18.99,mAdi18.99[length(mAdi18.99)]),base$age*365.2425)$y*SDAdi
    
    # Convert organ masses to volumes using tissue densities.
    VLun<-mLun/mean(c(1.04,1.056))
    VK  <-mK  /1.05
    VInt<-mInt/1.042
    VL  <-mL  /mean(c(1.045,1.056))
    VSk <-mSk /1.5
    VB  <-mB  /1.06
    VAdi<-mAdi/0.916
    
    # Calculate cardiac output and tissue plasma flows.
    sexF<-as.numeric(Gender=="F")
    EDV<-43.13+12.96*mB    # end diastolic volume ml (per beat)
    EF <-0.65-0.0018*BMI   # ejection fraction
    
    base$QC  =HR*EDV*EF/1000*60*24                                # cardiac output in L/d 
    base$QCP =base$QC*(1-physio$Htc)                              
    base$QL  =(25.53+1.30*sexF)/100*base$QCP
    base$QF  =exp(2.5-0.043*BMI+0.033*mAdi)/100*base$QCP          # according to Berton et al. (2022)
    base$QK  =(20.57-1.76*sexF)/100*base$QCP                      # according to Berton et al. (2022) 
    base$Qfil=QfilC*base$QK                                       # 0.2 according to EFSA 2020 and Lccisano, 2011
    base$QG  =(18.52+3.04*sexF-(0.20+0.06*sexF)*BMI+(0.0009+0.0004*sexF)*BMI^2)/100*base$QCP  # according to Berton et al. (2022)
    base$QSk =(5.68-0.034*BMI)/100*base$QCP                       # according to Berton et al. (2022)
    base$QR  =base$QCP-base$QL-base$QF-base$QK-base$QG-base$QSk   

    # Menstrual blood-loss clearance during the specified age interval.
    base$CLmenstruation=if_else(Gender=="F"&base$age>=settings$age_start_menstruation&base$age<=settings$age_stop_menstruation,
                                (physio$MFV*0.5*(1-Htc)+physio$MFV*0.5),0) # according to
    
    # Plasma, tissue, lung, and remaining-body volumes.
    base$Qp        =((((1400-190)*base$age^2.5)/(age^2.5+50))+190)*24
    base$Vplas     =physio$VplasC  *base$BW
    base$Vart_plas =physio$VartC   *base$Vplas
    base$Vven_plas =physio$VvenC   *base$Vplas
    #base$VL        =physio$VLC     *base$BW
    base$VL        =VL
    #base$VF        =physio$VFC     *base$BW
    base$VF        =VAdi                                          # fat = adipose
    #base$VK        =physio$VKC     *base$BW  
    base$VK        =VK
    #base$Vfil      =physio$VfilC   *base$BW                       
    base$Vfil      =VK*0.1                                        # according to EFSA, 2020; 10% of kidney volume
    #base$VG        =physio$VGC     *base$BW
    base$VG        =VInt                                          # gut = intestine
    #base$VSk       =(base$SkinTarea * physio$SkinThickness)/1000
    base$VSk       =VSk
    #base$Vlun      =physio$VlunC   *base$BW
    base$Vlun      =VLun
    base$VFRC      =0.03           *base$BW
    base$VT        =0.007          *base$BW
    base$Valv      =base$VFRC+0.5*base$VT
    base$VR        = 0.84*base$BW - base$VL-base$VF-base$VK-base$Vfil-base$VG-base$Vplas-base$VSk-base$Vlun
    #base$VR        = 
    
    # Skin-cell clearance.
    base$CLcell    =24*physio$Kpcell*physio$fss*base$SkinTarea/1000
    
    # Chemical-specific exposure and elimination rates.
    kin <-p$kinetics
    expo<-p$exposure
    mat <-p$maternal_c
    
    # Urinary and fecal elimination.
    base$kurine   <-kin$kurinec*base$BW^(-0.25)
    base$CLfaeces <-kin$CLfaeces*base$BW
    
    # Dermal transfer clearances.
    base$CLsc     <-24*p$Kpsc*physio$fss*base$SkinTarea/1000
    base$CLve     <-24*p$Kpve*physio$fss*base$SkinTarea/1000
    
    # Breast-milk exposure during infancy; oral exposure thereafter.
    base$month_idx<-floor(base$age*12)
    base$milk_conc<-mat$maternal*mat$Ratio*(1-mat$DECLINE)^base$month_idx
    base$milk_dose<-base$milk_conc*mat$Milkconsumption
    base$Oraldose <-if_else(base$age<1,base$milk_dose,expo$Oralexpo*base$BW)
    
    # Drinking-water exposure.
    base$Drinkdose<-(expo$Drinkconc*physio$Drinkrate/1000)*base$BW
    
    # Occupational inhalation exposure on working days and selected months.
    base$Inhalation<-if_else(base$dayoftheweek %in% 1:5 & base$month %in% c(1:4,11:12) & base$age >= settings$age_start_occ_exp & base$age <= settings$age_stop_occ_exp,
                             expo$Cinh/4,
                             0)
    # Occupational dermal exposure.
    base$Csurf<-if_else(base$dayoftheweek %in% 1:5 & base$month %in% c(1:4,11:12) & base$age >= settings$age_start_occ_exp & base$age <= settings$age_stop_occ_exp,
                        expo$Cdermal/4,
                        0)
    
    # Renal transporter capacity scaled to BW^0.75
    base$Tm<-kin$Tmc*base$BW^0.75
    
    return(base)
  }
  
  # 2) Convert time-varying variables to interpolation functions used by the ODE solver.
  
  make_interps <- function(Variables_df) {
    T  <- Variables_df$TIME
    
    # Time-varying model inputs.
    func_list <- list(
      #PFAS_params       = PFAS_params,
      physio            = physio,
      settings          = settings,
      varkurine         = approxfun(T, Variables_df$kurine,    rule=2),
      varCLfaeces       = approxfun(T, Variables_df$CLfaeces,  rule=2),
      varCLsc           = approxfun(T, Variables_df$CLsc,      rule=2),
      varCLve           = approxfun(T, Variables_df$CLve,      rule=2),
      varOraldose       = approxfun(T, Variables_df$Oraldose,  rule=2),
      varDrinkdose      = approxfun(T, Variables_df$Drinkdose, rule=2),
      varInhalation     = approxfun(T, Variables_df$Inhalation,rule=2),
      varCsurf          = approxfun(T, Variables_df$Csurf,     rule=2),
      varTm             = approxfun(T, Variables_df$Tm,        rule=2),
      
      #physiology functions
      varQCP           =approxfun(T, Variables_df$QCP,           rule=2),
      varQG            =approxfun(T, Variables_df$QG,            rule=2),
      varQL            =approxfun(T, Variables_df$QL,            rule=2),
      varQF            =approxfun(T, Variables_df$QF,            rule=2),
      varQK            =approxfun(T, Variables_df$QK,            rule=2),
      varQfil          =approxfun(T, Variables_df$Qfil,          rule=2),
      varQSk           =approxfun(T, Variables_df$QSk,           rule=2),
      varQR            =approxfun(T, Variables_df$QR,            rule=2),
      varCLmenstruation=approxfun(T, Variables_df$CLmenstruation,rule=2),
      varQp            =approxfun(T, Variables_df$Qp,            rule=2),
      varVart_plas     =approxfun(T, Variables_df$Vart_plas,     rule=2),
      varVven_plas     =approxfun(T, Variables_df$Vven_plas,     rule=2),
      varVG            =approxfun(T, Variables_df$VG,            rule=2),
      varVL            =approxfun(T, Variables_df$VL,            rule=2),
      varVF            =approxfun(T, Variables_df$VF,            rule=2),
      varVK            =approxfun(T, Variables_df$VK,            rule=2),
      varVfil          =approxfun(T, Variables_df$Vfil,          rule=2),
      varVSk           =approxfun(T, Variables_df$VSk,           rule=2),
      varVsc           =approxfun(T, Variables_df$Vsc,           rule=2),
      varVve           =approxfun(T, Variables_df$Vve,           rule=2),
      varCLcell        =approxfun(T, Variables_df$CLcell,        rule=2),
      varValv          =approxfun(T, Variables_df$Valv,          rule=2),
      varVlun          =approxfun(T, Variables_df$Vlun,          rule=2),
      varVR            =approxfun(T, Variables_df$VR,            rule=2))
    
    # Add static chemical parameters.
    static_parms <- parms_PFAS_extended()
    static_list  <- as.list(static_parms)
    return(c(func_list, static_list))
  }
  
  
  # 3) PBPK itself.
  PFAS_extended <- function(t, A, parms) {
    with(as.list(c(A, parms)), {
      ## kinetics --
      kurine <- varkurine(t)
      CLfaeces <- varCLfaeces(t)
      
      QCP <- varQCP(t)
      QG  <- varQG(t)
      QL  <- varQL(t)
      QF  <- varQF(t)
      QK  <- varQK(t)
      Qfil<- varQfil(t)
      QSk <- varQSk(t)
      QR  <- varQR(t)
      
      CLmenstruation <- varCLmenstruation(t)
      
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
      CLsc <- varCLsc(t)
      CLve <- varCLve(t)
      CLcell <- varCLcell(t)
      Valv <- varValv(t)
      Vlun <- varVlun(t)
      VR <- varVR(t)
      Oraldose <- varOraldose(t) # Dose expressed in ug/kg/day
      Drinkdose <- varDrinkdose(t) # Dose expressed in ug/kg/day
      Inhalation <- varInhalation(t) # Dose expressed in ug/kg/day
      Csurf <- varCsurf(t) # Dose expressed in ug/kg/day
      Tm <- varTm(t)
      
      # Concentrations
      #CAFree <- APlas/VPlas # free concentration of chemical in plasma µg/L (ng/mL) 
      #CA <- CAFree/Free # total concentration of chemical in plasma 
      CG <- AG/VG # Concentration in gut (µg/L) 
      CVG <- CG/PG # Concentration leaving gut (µg/L) 
      CL = AL/VL # Concentration in liver (µg/L) 
      CVL = CL/PL # Concentration leaving liver (µg/L) 
      CF = AF/VF # Concentration in fat (µg/L)
      CVF = CF/PF # Concentration leaving fat (µg/L) 
      CK = AK/VK # Concentration in kidneys (µg/L) 
      CVK = CK/PK # Concentration leaving kidneys (µg/L) 
      Cfil = Afil/Vfil # Concentration in filtrate compartment (µg/L) 
      CSk = ASk/VSk # Concentration in skin compartment (µg/L) 
      CVSk = CSk/PSk # Concentration leaving skin compartment (µg/L)  
      CR = AR/VR # Concentration in rest of the body (µg/L) 
      CVR = CR/PR # Concentration leaving rest of the body (µg/L) 
      Clun_bl <- Alun/(Valv*Pab + Vlun)
      Calv <- (Alun/(Valv*Pab + Vlun))*Pab
      # Csurf <- Asurf/Vsurf
      CvenFree <- Aven/Vven_plas
      Cven <- CvenFree/Free
      CartFree <- Aart/Vart_plas
      Cart <- CartFree/Free
      
      # Lung compartment
      dAlun <- QCP*Cven*Free + Qp*Inhalation - QCP*Clun_bl*Free #- Qp*Calv_PFOA
      
      # Skin compartment
      dASk <- QSk*(Cart*Free-CSk*FreeSk)    # Rate of change in skin(µg/h) 
      
      # Venous blood (plasma) compartment      
      dAven <- QF*CF*FreeF + (QL+QG)*CL*FreeL + QR*CR*FreeR + QSk*CSk*FreeSk + 
        QK*CK*FreeK - QCP*Cven*Free
      
      # Arterial blood (plasma) compartment      
      dAart <- QCP*Clun_bl*Free - QCP*Cart*Free - Qfil*Cart*Free - CLmenstruation*Cart*Free - CLfaeces*Cart*Free
      
      # Gut compartment 
      dAG <- QG*(Cart*Free-CG*FreeG) + Oraldose + Drinkdose #- CLfaeces*CG*FreeG
      
      # Liver compartment 
      dAL <- (QL*(Cart*Free)) + (QG*CG*FreeG) - ((QL+QG)*CL*FreeL) # Rate of change in liver (ug/h) 
      
      # Fat compartment 
      dAF <- QF*(Cart*Free-CF*FreeF)   # Rate of change in fat (µg/h) 
      
      # Kidney compartment 
      dAK <- QK*(Cart*Free-CK*FreeK) + (Tm*Cfil)/(Kt+Cfil) # Rate of change in kidneys (µg/h) 
      
      # Filtrate compartment 
      dAfil <- Qfil*(Cart*Free-Cfil) - (Tm*Cfil)/(Kt+Cfil) # Rate of change in filtrate compartment (ug/h) 
      
      # Storage compartment for urine 
      dAdelay <- Qfil*Cfil - kurine*Adelay   
      
      # Urine 
      dAurine <- kurine*Adelay 
      
      # Rest of the body 
      dAR <- QR*(Cart*Free-CR*FreeR)   # Rate of change in rest of the body (µg/h) 
      
      # Menstrual blood loss
      dAMenstruation <- CLmenstruation*Cart*Free
      
      # Fecal elimination
      dAfaeces <- CLfaeces*Cart*Free # CLfaeces*CG*FreeG # 
      
      # Mass balance
      Atot <- Alun + ASk + Aven + Aart +
        AG + AL + AF + AK + Afil + Adelay + Aurine + AR + AMenstruation + Afaeces
      
      list(c(dAlun, dASk, dAven, dAart,
             dAG, dAL, dAF, dAK, dAfil, dAdelay, dAurine, dAR, dAMenstruation, dAfaeces),
           c(Cart=Cart, Cven=Cven, Clun_bl=Clun_bl,
             Csurf=Csurf, CSk=CSk,
             CG=CG, CVG=CVG, CL=CL, CVL=CVL, CF=CF,CVF=CVF,
             CK=CK, CVK=CVK, Cfil=Cfil, CR=CR, CVR=CVR,
             Atot=Atot))
    })
  }
  
  
  parms_PFAS_extended <- function() {

    # Static transporter, free-fraction, and partition parameters.
    parms <- c(
      # transporter
      Kt = p$kinetics$Kt,
      # free + partitions
      Free = p$Free,
      FreeL = p$Free/p$partition$PL,
      FreeF = p$Free/p$partition$PF,
      FreeK = p$Free/p$partition$PK,
      FreeSk = p$Free/p$partition$PSk,
      FreeR = p$Free/p$partition$PR,
      FreeG = p$Free/p$partition$PG,
      FreeLun = p$Free/p$partition$PLun, #I added
      # partitions
      PL = p$partition$PL,
      PF = p$partition$PF,
      PK = p$partition$PK,
      PSk = p$partition$PSk,
      PR = p$partition$PR,
      PG = p$partition$PG,
      PLun = p$partition$PLun, #I added
      #Pab
      Pab = 1/(10^(6.96 - 1.04*log10(p$VP) - 0.533*p$logP - 0.00495*p$MW)),
      # Kscve & Kver
      Kscve = ((1-physio$ffatve) + physio$ffatve*10^p$logP) /
        ((1-physio$ffatsc) + physio$ffatsc*10^p$logP),
      Kver = ((1-physio$ffatbl) + physio$ffatbl*10^p$logP) /
        ((1-physio$ffatepi) + physio$ffatepi*10^p$logP),
      # and finally Tmc
      Tmc = p$kinetics$Tmc) #I added
  
    
    #names(parms) <- NULL
    return(parms)
  }
  
  # 5) Set initial chemical amounts in each compartment at birth.
  A_init_PFAS_extended <- function() {

    CONCbirth  <- p$maternal_c$maternal * p$maternal_c$PT             #CONCbirth
    APlasbirth <- CONCbirth * p$Free * physio$VplasC * physio$BWbirth #APlasbirth 
    art        <- p$birth$Aartbirth  * APlasbirth                     #I calculate here, Aartbirth
    ven        <- p$birth$Avenbirth  * APlasbirth                     #I calculate here, Avenbirth
    
    with(p$birth, {
      c(
        Alun   = ALunbirth   * APlasbirth, #ALunbirth
        ASk    = ASkbirth    * APlasbirth, #ASkbirth
        Aven   = ven,                      #Avenbirth
        Aart   = art,                      #Aartbirth
        AG     = AGbirth     * APlasbirth, #AGbirth
        AL     = ALbirth     * APlasbirth, #ALbirth
        AF     = AFbirth     * APlasbirth, #AFbirth
        AK     = AKbirth     * APlasbirth, #AKbirth
        Afil   = Afilbirth   * APlasbirth, #Afilbirth
        Adelay = Adelaybirth * APlasbirth, #Adelaybirth
        Aurine = Aurinebirth * APlasbirth, #Aurinebirth
        AR     = ARbirth     * APlasbirth, #ARbirth
        AMenstruation = 0,
        Afaeces       = 0)
    })
  }
  
  Variables_df<- Variables_df(age, BW, height, Gender)
  A_init    <- A_init_PFAS_extended()
  pars_list <- make_interps(Variables_df)
  
  out <- deSolve::lsoda(
    y      = A_init,
    times  = Variables_df$TIME,
    func   = PFAS_extended,
    parms  = pars_list)
  
  out %>%
    as.data.frame() #%>%
    #as_tibble()
  
  # Return the arterial plasma concentration at the final simulated time.
  return(out[nrow(out),"Cart"])
}