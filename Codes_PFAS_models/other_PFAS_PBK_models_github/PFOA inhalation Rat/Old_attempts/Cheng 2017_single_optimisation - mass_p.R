library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    # BW in kg
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
    #assuming the density of tissue is 1 g/mL.
    
    #units conversion from Cheng 2017R, time-> h, PFOA mass->ng, tissues mass-> g
    
    #======Table S1=======#    
    
    #Blood
    PVB <- 54e-3 #13.5 mL/244 g=0.055 mL/g~55e-3 mL/g (kg=L)
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 31.2e-3 
    Vplasma <- PVplasma * BW #plasma volume kg=L
    Vven <- BW*11.3/250 	#volume of venous plasma (L); from doi:10.1007/bf02353860
    Vart <- BW*5.6/250	#volume of arterial plasma (L); from doi:10.1007/bf02353860
    #Kidney
    PVK <- 7.3e-3 
    VK <- PVK * BW #kidney volume kg=L
    PVKB <- 0.16 
    VKB <- PVKB * PVK * BW #kidney blood volume kg=L
    PVKF <- 0.13 
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume kg=L
    VKT <- VK - VKF #kidney tissue volume kg=L
    VFil <- 0.25 #renal filtrate volume kg=L
    #Liver
    PVL <- 3.66e-2 
    VL <- PVL * BW #liver volume kg=L
    PVLB <- 0.21 
    VLB <- PVLB * PVL * BW #liver blood volume kg=L
    PVLF <- 0.049 
    VLF <- PVLF * PVL* BW #liver interstitial fluid volume kg=L
    VLT <- VL - VLF #liver tissue volume kg=L
    PVbile <- 0.004 
    Vbile <- PVbile * PVL * BW #bile volume kg=L
    #Gut
    PVG <- 2.69e-2 
    VG <- PVG * BW #gut volume kg=L
    PVGB <- 0.034
    VGB <- PVGB * PVG * BW #gut blood volume kg=L
    PVGF <- 0.28
    VGF <- PVGF * PVG * BW #gut interstitial fluid volume kg=L
    VGT <- VG - VGF #gut tissue volume kg=L
    PVGL <- 4.5e-2 
    VGL <- PVGL * BW #gut lumen volume kg=L
    #Muscle
    PVM <- 40.43e-2 
    VM <- PVM * BW #muscle volume kg=L
    PVMB <- 0.04 
    VMB <- PVMB * PVM * BW #muscle blood volume kg=L
    PVMF <- 0.054
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume kg=L
    VMT <- VM - VMF #muscle tissue volume kg=L
    #Adipose
    PVA <- 7e-2 
    VA <- PVA * BW #adipose volume kg=L
    PVAB <- 0.02 
    VAB <- PVAB * PVA * BW #% adipose blood volume kg=L
    PVAF <- 0.174
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume kg=L
    VAT <- VA - VAF #adipose tissue volume kg=L
    #RoB
    PVR <- 1 - PVB - PVK - PVL - PVG - PVM - PVA
    VR <- PVR * BW #volume of the rest of the body kg=LL
    PVRB <- 0.036 
    VRB <- PVRB * PVR * BW #volume of the blood of the rest of body kg=L
    PVRF <- 0.18  
    VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body kg=L
    VRT <- VR - VRF #tissue volume of the rest of body kg=L
    #Lymph nodes
    PVLN <- 1.15/280 #Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLN <- PVLN*BW#lymph fluid volume
    #Lung
    PVLu <- 0.48 #Brown et al. 1997, p 418-19 mean of values for male rat
    VLu <- PVLu * BW
    PVLuB <- 0.09*0.55 #Brown et al. 1997, p 459 --> capillary blood occupied 9% of the lung volume
    VLuB <-PVLuB*BW #volume of the blood of lung kg=L
    PVLuF <- 0.263/280 #Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLuF <- PVLuF*BW #lung interstitial fluid volume
    PVLuAF <- 0.4/1000/275 #0.4 mL Leslie et al, 1989 https://doi.org/10.1164/ajrccm/139.2.360 --> Watkins & Rannels 1979 https://doi.org/10.1152/jappl.1979.47.2.325  
    VLuAF <- PVLuAF*BW #lung alveolar lining fluid volume kg=LL
    PVLuT <- VLu - VLuF - VLuAF # einai ok ???
    VLuT <- PVLuT * BW #lung tissue volume kg=L
      
    #Capillary surface area for each tissue (Ai) as percentage of body weight
    #or weight of corresponding tissue (PAi, m^2/g) and surface area (m^2)
    
    PAK <- 350e-4 #cm2/g   cm2*L/g-> 1e4 m2*1000
    AK <- PAK * VK * 1e7 #kidney surface area (m^2)
    PAKG <- 68.90e-4 
    AKG <- PAKG * VK * 1e7 #the surface area of glomerular capillary (m^2)
    PAL <- 250e-4 
    AL <- PAL * VL * 1e7 #liver surface area (m^2)
    PAG <- 100e-4 
    AG <- PAG * VG * 1e7 #gut surface area (m^2)
    PAGL <- 4.14 #"The surface area of gut lumen would be 4.14 m2/kg" !!!!???
    AGL <- PAGL * BW #gut lumen surface area (m^2)
    PAM <- 70e-4 
    AM <- PAM * VM * 1e7 #muscle surface area (m^2)
    PAA <- 70e-4
    AA <- PAA * VA * 1e7 #adipose surface area (m^2)
    PAR <- 100e-4
    AR <- PAR * VR * 1e7 #surface area of rest of body (m^2)
    PLu <- 250e-4
    ALu <- PLu* VLu * 1e7 #surface area of rest of body (m^2)
    
    #Effective permeability (Peff, in m/h) for blood (B), liver(L), kidney(K),
    #gut(G),adipose(A), muscle(M), rest of body(R)

    PeffB <- 4.98e-8*3600
    PeffK <- 4.38e-8*3600
    PeffL <- 5.15e-8*3600
    PeffG <- 2.65e-8*3600
    PeffA <- 2.65e-8*3600
    PeffM <- 2.65e-8*3600
    PeffR <- 2.65e-8*3600
    PeffLu <- 2.65e-8*3600 #assumption
    
    
    
    #Blood flow rates (QBi, in L/h) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac L/h), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 0.235 * (BW^0.75) *60 #L/min->*60-> L/h
    PQBK <- 14.1/100 
    QBK <- PQBK * Qcardiac #L/h
    PQBG <- 15.1/100 
    QBG <- PQBG * Qcardiac #L/h
    PQBL <- 2.4/100 
    QBL <- (PQBL+PQBG) * Qcardiac #L/h
    PQBM <- 27.8/100 
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 7/100 
    QBA <- PQBA * Qcardiac #L/h
    PQBR = 1 - PQBK - PQBG - PQBL - PQBM - PQBA
    QBR = PQBR * Qcardiac #L/h
    PQBLu <- 1#2.1/100
    QBLu <- PQBLu * Qcardiac #L/h
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    
    Qfeces <- 5.63/1000/24 #mL water/d --> L/h
    PQbile <- 90/1000/24 #mL/d/kg BW --> L/d/kg BW
    Qbile <- PQbile * BW #L/h
    PQurine <- 200/1000/24 #mL/d/kg BW --> L/h/kg
    Qurine <- PQurine * BW #L/h
    PQGFR <- 10.74*60/1000 #mL/min/kg  -->L/h/kg    
    QGFR <- PQGFR * BW #L/h
    
    #======Table S2=======#
    #Albumin concentration in blood and interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)

    CalbB <- 281e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbKF <- 243e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbLF <- 243e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbGF <- 146e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbMF <- 146e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbAF <- 73e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbRF <- 73e-3*7.8 #n=7.8 binding sites (mol/m3)
    CalbLN <- CalbGF #assumption based on https://doi.org/10.1016/j.jconrel.2020.07.046
    
    CalbLuF <- CalbGF #assumption 
    CalbLuAF <- 10/100 * CalbB #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
      
    #Alpha2mu-globulin concentration in kidney tissue (mol/m3)

    Ca2uKT <- 110e-3

    #LFABP concentration in kidney and liver tissue (mol/m^3)

    CLfabpKT <- 2.65e-3 * 3 #n=3 binding sites (mol/m3)
    CLfabpLT <- 133e-3 * 3 #n=3 binding sites (mol/m3)
  

    #======Table S2=======#
    #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
    #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details

    Ka <-  24.18    # 3.1*7.8 m3/mol (M-1) multiplying by number of binding sites (Cheng et al. 2021)
    KLfabp <- 135  # 45.0*3 geo_mean of 3 binding affinity, may try normal mean (Cheng et al. 2021)
    Ka2u <- 0.5 #m3/mol

    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #passive diffusion rates
    
    ClLFT_unscaled= 67.8#uL/min/10^6 cells, Han et al. 2008
    ClKFT_unscaled= 17.5 #uL/min/mg protein, Yang et al. 2010
    ClGFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    ClMFT_unscaled= 18.1 #same as ClGFT
    ClAFT_unscaled= 18.1 #same as ClGFT
    ClRFT_unscaled= 18.1 #same as ClGFT
    ClLuFT_unscaled= 18.1 #same as ClGFT
    
    
    #For all CMTs
    MW <- 414.07 #ug/umol, PFOA molecular weight
    Acell = 4000 #um^2/cell
    
    #Kidney
    kidney_protein_per_rat <- 1000*(0.218+0.225+0.212)/3#mg of protein per rat  (Addis 1936)
    rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
    rat_kidney_weight_addis <- rat_weight_addis*0.0073 # kidney fraction to BW, Brown (1997)
    kidney_protein_per_gram <- kidney_protein_per_rat/rat_kidney_weight_addis #mg of protein/g kidney
    
    kidney_cells = 1.47e07 #cells/g https://doi.org/10.1038/s41598-024-53270-2
    kidney_cells_total <- kidney_cells* (1000*VKT)
    kidney_protein_total <- kidney_protein_per_gram* (1000*VKT) #mg
    ClKFT <- ClKFT_unscaled * kidney_protein_total #uL/min for the whole kidney
    # Uptake for the whole kidney tissue 
    kKFKT <- (60*ClKFT)/1e06 #L/h
    
    #Cmedium_K = 4*MW# Yang et al. 2010 4uM umol/L -->  ug/L
    # PeffKT = (60*ClKFT/1000) /Acell #m/hour/cell
    # kKFKT <- PeffKT*(Acell/10^6)*kidney_cells_total*1000 #L/h
    
    RAFOatp_k <- estimated_params[1]
    RAFOat1 <- estimated_params[2]
    RAFOat3 <-  estimated_params[3]
    
    #kBKF <- (((1/QBK) + 1/(1000*PeffB * AK))^(-1)) #multiplication is to convert m3 -> L
    n <- 5 #enlargement factor of apical membrane of proximal tubule
    kFKT <- PeffK * AK * n
    #kBF <- PeffB * AKG
    
    
    #Oatp kidney
    VmK_Oatp_in_vitro <- 9.3 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oatp_scaled <- 60*VmK_Oatp_in_vitro*MW*kidney_protein_per_gram/1000  #physiologically scaled to in vivo, ug/L/h
    VmK_Oatp <- VmK_Oatp_scaled*RAFOatp_k #in vivo value, in  ug/L/h
    KmK_Oatp=  126.4*MW#umol/L (Weaver et al. 2010) --> ug/L
    
    #VmK_Osta= 9.3*2.0e-6*75960 #nmol/min (Weaver et al. 2010 - same as Oatp1)
    #KmK_Osta= 52.3 #mg/mL (Weaver et al. 2010 - mean of Oat1,3,Oatp)
    
    #oat1 kidney
    VmK_Oat1_in_vitro= 2.6 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oat1_scaled = 60*VmK_Oat1_in_vitro*MW*kidney_protein_per_gram/1000  #physiologically scaled to in vivo, ug/L/h
    VmK_Oat1= VmK_Oat1_scaled*RAFOat1 #in vivo value, in  ug/L/h
    KmK_Oat1= 43.2 * MW #umol/L (Weaver et al. 2010) --> ug/L
    
    #oat3 kidney
    VmK_Oat3_in_vitro= 3.8 #nmol/mg protein/min  (Weaver et al. 2010)
    VmK_Oat3_scaled = 60*VmK_Oat3_in_vitro*MW*kidney_protein_per_gram/1000  #physiologically scaled to in vivo, ug/L/h
    VmK_Oat3 = VmK_Oat3_scaled*RAFOat3 #in vivo value, in  ug/L/h
    KmK_Oat3= 65.7 * MW #umol/L (Weaver et al. 2010) --> ug/L

    #Liver
    liver_protein_per_rat <- 1000*(1.52+1.53+1.52)/3#mg of protein per rat  (Addis 1936)
    rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
    rat_liver_weight_addis <- rat_weight_addis*0.0366 # liver fraction to BW, Brown (1997)
    liver_protein_per_gram <- liver_protein_per_rat/rat_liver_weight_addis #mg or protein/g liver
    liver_cells = 117*10^6 #hepatocytes per g of liver (Sohlenius-Sternbeck et al. 2006) (2e09 cells: https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=110895)
    
    ClLFT <- ClLFT_unscaled*(liver_cells/(10^6))*(VLT*1000) #uL/min for the whole liver
    kLFLT <-  (60*ClLFT)/1e06 #L/h
    #Cmedium_L = 1*MW# Han et al. 2008 1uM umol/L -->  ug/L
    #PeffLT = 1000*60*ClLFT*liver_protein_per_rat/(Acell*Cmedium_L) #m/h
    #kLFLT <- PeffLT*Acell*liver_cells/1000 #L/h
    #kBLF <- ((1/QBL) + 1/(1000*PeffB * AL))^(-1) #multiplication is to convert m3 -> L
    kbileLT <- PeffL *estimated_params[6]* AL
    
    RAFOatp_l <- estimated_params[4]
    RAFNtcp <- estimated_params[5]
    
    # oatp-liver
    VmL_Oatp_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
    VmL_Oatp_scaled = 60*VmL_Oatp_in_vitro*MW*liver_protein_per_gram/1000  #physiologically scaled to in vivo, ug/L/h
    VmL_Oatp = VmL_Oatp_scaled*RAFOatp_l #in vivo value, in  ug/L/h
    KmL_Oatp = KmK_Oatp #same as kidney
    
    #Ntcp liver
    VmL_Ntcp_in_vitro= 3#nmol/mg protein/min  (Weaver et al. 2010)
    VmL_Ntcp_scaled = 60*VmL_Ntcp_in_vitro*MW*liver_protein_per_gram/1000  #physiologically scaled to in vivo, ug/L/h
    VmL_Ntcp = VmL_Ntcp_scaled*RAFNtcp #in vivo value, in  ug/L/h
    KmL_Ntcp= 20 * MW #umol/L, Ruggiero et al. 2021 --> ug/L
    
    #Muscle
    muscle_cells= NA
    Cmedium_M = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    muscle_protein <- 158.45 #mg/g muscle (protein data from Cheek et al.,1971 and muscle mass from Caster et al.,1956)
    muscle_protein_tot <- muscle_protein * (VMT*1000)
    ClMFT <- ClMFT_unscaled * muscle_protein_tot#uL/min for the whole muscle comp
    kMFMT <-  (60*ClMFT)/1e06 #L/h
    
    #PeffMT = 1000*60*ClMFT*muscle_protein/(Acell*Cmedium_M) #m/h
    #kBMF <- ((1/QBM) + 1/(1000*PeffB * AM))^(-1) #multiplication is to convert m3 -> L
    #kMFMT <- ClMFT*muscle_protein*VMT/60
    
    
    #Gut
    gut_cells = NA
    Cmedium_G = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    gut_protein <- muscle_protein#NA#5034
    gut_protein_total <- gut_protein*(1000*VGT)
    ClGFT <- ClGFT_unscaled *gut_protein_total#uL/min for the whole gut compartment
    kGFGT <-  (60*ClGFT)/1e06 #L/h
    
    #PeffGT = 1000*60*ClGFT*gut_protein/(Acell*Cmedium_G) #m/h
    #kGFGT <- PeffGT*Acell*gut_cells/1000 #L/h
    #kBGF <- ((1/QBG) + 1/(1000*PeffB * AG))^(-1) #multiplication is to convert m3 -> L
    kGLGT <- PeffG * AGL
   
    #Adipose
    adipose_cells = NA
    Cmedium_A = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    adipose_protein <- muscle_protein#NA#5034
    adipose_protein_total <- adipose_protein * (1000*VAT)
    ClAFT <- ClAFT_unscaled * adipose_protein_total#uL/min
    kAFAT <-  (60*ClAFT)/1e06 #L/h
    
    #PeffAT = 1000*60*ClAFT*adipose_protein/(Acell*Cmedium_A) #m/h
    #kAFAT <- PeffAT*Acell*adipose_cells/1000 #L/h
    #kBAF <- ((1/QBA) + 1/(1000*PeffB * AA))^(-1) #multiplication is to convert m3 -> L
    
    #Rest of body
    RoB_cells = NA
    Cmedium_R = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    RoB_protein <- muscle_protein#18456
    RoB_protein_total <- RoB_protein * (1000* VRT)
    ClRFT <- ClRFT_unscaled * RoB_protein_total#uL/min
    kRFRT <-  (60*ClRFT)/1e06 #L/h
    
    #PeffRT = 1000*60*ClRFT*RoB_protein/(Acell*Cmedium_R) #m/h
    #kRFRT <- PeffRT*Acell*RoB_cells/1000 #L/h
    #kBRF <- ((1/QBR) + 1/(1000*PeffB *AR))^(-1) #multiplication is to convert m3 -> L
    
    
    
    #Rest of body
    Lung_protein <- muscle_protein#18456
    Lung_protein_total <- Lung_protein * (1000* VRT)
    ClLuFT <- ClLuFT_unscaled * Lung_protein_total#uL/min
    kLuFLuT <-  (60*ClLuFT)/1e06 #L/h
    
    return(list('PVB'=PVB, 'VB'=VB, 'PVplsma'=PVplasma, 
                'Vplasma'=Vplasma, 'PVK'=PVK, 'VK'=VK, 'PVKB'=PVKB, 'VKB'=VKB, 
                'PVKF'=PVKF, 'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'PVL'=PVL, 'VL'=VL, 'PVLB'=PVLB, 'VLB'=VLB, 'PVLF'=PVLF, 
                'VLF'=VLF, 'VLT'=VLT, 'PVbile'=PVbile, 'Vbile'=Vbile, 'PVG'=PVG, 
                'VG'=VG, 'PVGB'=PVGB, 'VGB'=VGB,
                'PVGF'=PVGF,'VGF'=VGF, 'VGT'=VGT, 'PVGL'=PVGL, 'VGL'=VGL,
                'PVM'=PVM, 'VM'=VM, 'PVMB'=PVMB, 'VMB'=VMB, 'PVMF'=PVMF, 'VMF'=VMF, 
                'VMT'=VMT,'PVA'=PVA, 'VA'=VA, 'PVAB'=PVAB, 'VAB'=VAB, 'PVAF'=PVAF, 
                'VAF'=VAF, 'VAT'=VAT, 'PVR'=PVR, 'VR'=VR, 'PVRB'=PVRB, 'VRB'=VRB, 
                'PVRF'=PVRF, 'VRF'=VRF, 'VRT'=VRT, 'VLN' = VLN, 'Vven' = Vven,
                'Vart' = Vart, 'VLu'=VLu, 'PVLuB'=PVLuB, 'VLuB'=VLuB, 'PVLuF'=PVLuF, 'VLuF'=VLuF,
                'PVLuAF'=PVLuAF, 'VLuAF'=VLuAF, 'PVLuT'=PVLuT, 'VLuT'=VLuT,
                
                'AK'=AK,'AKG'=AKG, 'AL'=AL, 'AG'=AG, 'AGL'=AGL,
                'AM'=AM,  'AA'=AA,  'AR'=AR, 'ALu' = ALu, 
                
                'PeffB'=PeffB, 'PeffK'=PeffK, 'PeffL'=PeffL, 'PeffG'=PeffG, 
                'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
                
                'Qcardiac'=Qcardiac, 'PQBK'=PQBK, 'QBK'=QBK, 'PQBG'=PQBG, 'QBG'=QBG, 
                'PQBL'=PQBL, 'QBL'=QBL, 'PQBM'=PQBM, 'QBM'=QBM, 'PQBA'=PQBA, 'QBA'=QBA,
                'PQBR'=PQBR, 'QBR'=QBR, 'PQBLu'=PQBLu, 'QBLu'=QBLu,
                'Qfeces'=Qfeces, 'PQbile'=PQbile, 'Qbile'=Qbile, 'PQurine'=PQurine,
                'PQGFR'=PQGFR, 'QGFR'=QGFR,'Qurine'=Qurine, 'QGFR'=QGFR,
                
                'CalbB'=CalbB, 'CalbKF'=CalbKF, 'CalbLF'=CalbLF, 'CalbGF'=CalbGF, 
                'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF,'CalbLN' = CalbLN,
                'CalbLuF' = CalbLuF,
                
                'Ca2uKT'=Ca2uKT,'CLfabpKT'=CLfabpKT,'CLfabpLT'=CLfabpLT, 
                
                'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,
                
                'ClLFT'=ClLFT, 'ClKFT'=ClKFT, 'ClGFT'=ClGFT, 'ClMFT'=ClMFT,
                'ClAFT'=ClAFT, 'ClRFT'=ClRFT, 
                
                'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp,
                'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 
                'KmK_Oatp'=KmK_Oatp,
                #'VmK_Osta'=VmK_Osta,'KmK_Osta'=KmK_Osta, 
                'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
                'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
                  
                
                'kKFKT'=kKFKT, 'kFKT'=kFKT,  
                'kLFLT'=kLFLT, 'kbileLT'=kbileLT,  'kAFAT'=kAFAT, 
                'kRFRT'=kRFRT,  'kGFGT'=kGFGT,'kGLGT'=kGLGT, 
               'kMFMT'=kMFMT, 'kLuFLuT' =kLuFLuT,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type
                ))
  
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{

    #====================PFOA mass balance at each tissue or fluid compartment==============================     
      
    # Concentrations in ug/L
    
    # Blood concentration
    CBven <- MBven/Vven
    CBart <- MBart/Vart
    
    #Lymph node
    CLN <- MLN /VLN
    
    # Kidney 
    CKF = MKF/VKF #interstitial fluid concentration
    CKB = MBK/VKB # blood concentration
    CKT = MKT/VKT # tissue concentration
    CFil = MFil/VFil# filtrate concentration
    
    #Liver
    CLF = MLF/VLF #interstitial fluid concentration
    CLB = MBL/VLB # blood concentration
    CLT = MLT/VLT # tissue concentration
    Cbile = Mbile/Vbile 
    
    #Gut
    CGB = MBG/VGB # blood concentration
    CGF = MGF/VGF  #interstitial fluid concentration
    CGT = MGT/VGT # tissue concentration
    CGL = MGL/VGL # Lumen concentration
    
    # Muscle
    CMB = MBM/VMB # blood concentration
    CMF = MMF/VMF  #interstitial fluid concentration
    CMT = MMT/VMT # tissue concentration
    
    #Adipose
    CAB = MBA/VAB # blood concentration
    CAF = MAF/VAF  #interstitial fluid concentration
    CAT = MAT/VAT # tissue concentration
    
    #Rest-of-the-body
    CRB = MBR/VRB # blood concentration
    CRF = MRF/VRF  #interstitial fluid concentration
    CRT = MRT/VRT # tissue concentration
    
    #Lung
    CLuB = MBLu/VLuB # blood concentration
    CLuF = MLuF/VLuF  #interstitial fluid concentration
    CLuT = MLuT/VLuT # tissue concentration
    CLuAF = MLuAF/VLuAF #alveolar lining fluid concentration
    
    #Cfree calculation using the expression of free fraction ff
    CBfart = CBart * 1.0 / (1.0 + CalbB * Ka)
    CBfven= CBven * 1.0 / (1.0 + CalbB * Ka)
    CLNf = CLN * 1.0 / (1.0 + CalbB * Ka)  # find albumin concentration in Lymph Nodes
    
    #Calculation of free concentrations in organ blood
    CKBf = CKB * 1.0 / (1.0 + CalbB * Ka)
    CLBf = CLB * 1.0 / (1.0 + CalbB * Ka)
    CGBf = CGB * 1.0 / (1.0 + CalbB * Ka)
    CMBf = CMB * 1.0 / (1.0 + CalbB * Ka)
    CABf = CAB * 1.0 / (1.0 + CalbB * Ka)
    CRBf = CRB * 1.0 / (1.0 + CalbB * Ka)
    CLuBf = CLuB * 1.0 / (1.0 + CalbB * Ka)
    
    #Calculation of free concentrations in organ interstitial fluid
    CKFf = CKF * 1.0 / (1.0 + CalbKF * Ka)
    CLFf = CLF * 1.0 / (1.0 + CalbLF * Ka)
    CGFf = CGF * 1.0 / (1.0 + CalbGF * Ka)
    CMFf = CMF * 1.0 / (1.0 + CalbMF * Ka)
    CAFf = CAF * 1.0 / (1.0 + CalbAF * Ka)
    CRFf = CRF * 1.0 / (1.0 + CalbRF * Ka)
    CLuFf = CLuF * 1.0 / (1.0 + CalbLuF * Ka)
    
    #Calculation of free concentrations in organ where we have tissue binding
    CKTf = CKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp)
    CLTf = CLT * 1.0 / (1.0 + CLfabpLT * KLfabp)
    
    #Arterial Blood
    dMBart = (QBLu-QBLu/500)*CLuBf - CBfart*(QBK+QBL+ QBG+ QBM+ QBA+ QBR)-QGFR*CBfart
    
    #Venous Blood
    dMBven = - CBfven *QBLu +  CLNf*(QBK/500+QBL/500+QBG/500+QBM/500+QBA/500+QBR/500+QBLu/500)+
              (QBK-QBK/500)*CKBf+(QBL-QBL/500)*CLBf+(QBG-QBG/500)*CGBf +
               (QBM-QBM/500)*CMBf+ (QBA/500)*CABf+(QBR-QBR/500)*CRBf
      
    #Lymph nodes
    dMLN = CKFf*QBK/500 + CLFf*QBL/500 +  CGFf*QBG/500 + CMFf*QBM/500 + CAFf*QBA/500+
      CRFf*QBR/500 +  CLuFf*QBLu/500 - CLNf*(QBK/500+QBL/500+QBG/500+QBM/500+QBA/500+QBR/500+QBLu/500)
    
    
    #Kidney
    
    #blood subcompartment
    dMBK = QBK*CBfart - (QBK-QBK/500)*CKBf - PeffK*AK*(CKBf-CKFf) - CKBf*QBK/500 
    #interstitial fluid subcompartment
    dMKF = CKBf*QBK/500 - CKFf*QBK/500 + PeffK*AK*(CKBf-CKFf) - kKFKT*(CKFf-CKTf) -
            (VmK_Oat1*CKFf/KmK_Oat1+CKFf)*VKF - (VmK_Oat3*CKFf/KmK_Oat3+CKFf)*VKF #+ (VmK_Osta*CKTf/KmK_Osta+CKTf)
    #Kidney proximal tubule cells subcompartment
    dMKT = kKFKT*(CKFf-CKTf) - kFKT*(CKTf - CFil) + (VmK_Oatp*CFil/KmK_Oatp+CFil)*VFil +
            (VmK_Oat1*CKFf/KmK_Oat1+CKFf)*VKF + (VmK_Oat3*CKFf/KmK_Oat3+CKFf)*VKF # + (VmK_Osta*CKTf/KmK_Osta+CKTf)
    dMFil =  QGFR*CBfart+ kFKT*(CKTf - CFil) - (VmK_Oatp*CFil/KmK_Oatp+CFil)*VFil - (Qurine/VFil)*CFil
    dMurine = (Qurine/VFil)*CFil
    
    #Liver
    
    #blood subcompartment
    dMBL = QBL*CBfart - (QBL-QBL/500)*CLBf - PeffL*AL*(CLBf-CLFf) - CLBf*QBL/500
    #interstitial fluid subcompartment 
    dMLF = CLBf*QBL/500 - CLFf*QBL/500 + PeffL*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
          (VmL_Oatp*CLFf/KmL_Oatp+CLFf)*VLF - (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf)*VLF
    #Liver tissue subcompartment
    dMLT = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/KmL_Oatp+CLFf)*VLF + 
                     (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf)*VLF - kbileLT*(CLTf-Cbile)
    dMbile = kbileLT*(CLTf-Cbile) - (Qbile/Vbile)*Cbile
    
    
    #Gut
    #blood subcompartment
    dMBG = QBG*CBfart - (QBG-QBG/500)*CGBf - PeffG*AG*(CGBf-CGFf) - CGBf*QBG/500
    #interstitial fluid subcompartment 
    dMGF = CGBf*QBG/500 - CGFf*QBG/500 + PeffG*AG*(CGBf-CGFf) - kGFGT*(CGFf-CGT) 
    #Gut tissue subcompartment
    dMGT = kGFGT*(CGFf-CGT)  - kGLGT*(CGT-CGL)
    # Gut lumen
    dMGL = kGLGT*(CGT-CGL) + (Qbile/Vbile)*Cbile - (Qfeces/VGL)*CGL
    # Feces
    dMfeces = (Qfeces/VGL)*CGL
    
    #Muscle
    
    #blood subcompartment
    dMBM = QBM*CBfart - (QBM-QBM/500)*CMBf - PeffM*AM*(CMBf-CMFf) - CMBf*QBM/500
    #interstitial fluid subcompartment 
    dMMF = CMBf*QBM/500 - CMFf*QBM/500 +PeffM*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMT)
    #Muscle tissue subcompartment 
    dMMT = kMFMT*(CMFf- CMT)
    
    #Adipose
    
    #blood subcompartment
    dMBA = QBA*CBfart - (QBA/500)*CABf - PeffA*AA*(CABf-CAFf) - CABf*QBA/500
    #interstitial fluid subcompartment 
    dMAF = CABf*QBA/500 - CAFf*QBA/500 +  PeffA*AA*(CABf-CAFf) - kAFAT*(CAFf-CAT) 
    #Adipose tissue subcompartment 
    dMAT =  kAFAT*(CAFf-CAT) 
    
    #Rest of body
    
    #blood subcompartment
    dMBR = QBR*CBfart - (QBR-QBR/500)*CRBf - PeffR*AR*(CRBf-CRFf) - CRBf*QBR/500
    #interstitial fluid subcompartment 
    dMRF = CRBf*QBR/500 - CRFf*QBR/500 + PeffR*AR*(CRBf-CRFf) - kRFRT*(CRFf -CRT) 
    #Rest of body tissue subcompartment 
    dMRT = kRFRT*(CRFf -CRT) 
    
    
    #Lung  

    #blood subcompartment
    dMBLu = CBfven *QBLu - (QBLu-QBLu/500)*CLuBf - PeffLu*ALu*(CLuBf-CLuFf) - CLuBf*QBLu/500
    #interstitial fluid subcompartment
    dMLuF = CLuBf*QBLu/500 - CLuFf*QBLu/500 + PeffLu*ALu*(CLuBf-CLuFf) - kLuFLuT*(CLuFf-CLuT)
    #Lung tissue
    dMLuT =  kLuFLuT*(CLuFf-CLuT) #- kLuTLuAF * (CLuAF-CLuT)
    #Alveolar lining fluid
    dMLuAF = 0#kLuTLuAF * (CLuAF-CLuT) - Qp * CLuAF/PLungA + IVR*Cair*Cpfoa*dfalveolar
    #depot
    #dMDep = Qp * CLuAF/PLungA
  
    
    Cven <- CBven
    Cart <- CBart
    Cblood <- (MBven  + MBart)/ VB
    Ckidney <- (MBK + MKF+ MKT)/(VKB+VKF+VKT)
    Cliver <- (MBL + MLF+ MLT)/(VLB+VLF+VKT)
    Cgut <-  (MGL + MGF+ MGT)/(VGB+VGF+VGT)
    Cmuscle <-  (MBM + MMF+ MMT)/(VMB+VMF+VMT)
    Cadipose <-  (MBA + MAF+ MAT)/(VAB+VAF+VAT)
    Clungs <-  (MBLu + MLuF+ MLuT)/VLu
    
    Crest <-  (MBR + MRF+ MRT)/(VRB+VRF+VRT)
    Ccarcass <- (MBM + MMF+ MMT+MBA + MAF+ MAT +MBR + MRF+ MRT)/(VM+VA+VR) 
    Cfeces <- CGL
    Cbile <- Cbile
    Curine <- CFil
    
    list(c('dMBart'=dMBart, 'dMBven'=dMBven, 'dMLN'=dMLN, 'dMBK'=dMBK, 
            'dMKF'=dMKF, 'dMKT'=dMKT,
                  'dMFil'=dMFil, 'dMurine'=dMurine, 'dMBL'=dMBL, 
           'dMLF'=dMLF, 'dMLT'=dMLT, 'dMbile'=dMbile,
                  'dMBG'=dMBG, 'dMGF'=dMGF, 'dMGT'=dMGT,
           'dMGL'=dMGL, 'dMfeces'=dMfeces, 'dMBM'=dMBM, 
                  'dMMF'=dMMF, 'dMMT'=dMMT, 'dMBA'=dMBA, 'dMAF'=dMAF, 
           'dMAT'=dMAT, 'dMBR'=dMBR, 'dMRF'=dMRF,'dMRT'=dMRT,
                  'dMBLu'=dMBLu, 'dMLuF'=dMLuF,'dMLuT'=dMLuT,'dMLuAF' = dMLuAF),
                  
                  'CKFf'=CKFf, 'CLNf'=CLNf, 'CLFf'=CLFf, 'CGFf'=CGFf, 'CMFf'=CMFf,
         'CAFf'=CAFf, 'CRFf'=CRFf, 'CBfart'=CBfart, 
                  'CKBf'=CKBf, 'CLBf'=CLBf, 'CGBf'=CGBf, 'CMBf'=CMBf, 'CABf'=CABf,
         'CRBf'=CRBf, 'CFil'=CFil, 'Cbile'=Cbile, 
                  'Cfeces'=Cfeces, 'CKTf'=CKTf, 'CLTf'=CLTf, 'CGT'=CGT, 'CGL'=CGL,
         'CMT'=CMT, 'CAT'=CAT, 'CRT'=CRT, 
                  
                  'Cven'=Cven, 'Cart' = Cart,'Cblood' = Cblood,
                'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cgut'=Cgut,
                  'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose, 
         'Clungs' = Clungs, 'Crest'=Crest,'Ccarcass' = Ccarcass,
                  'Cfeces'=Cfeces, 'Cbile'=Cbile, 'Curine'=Curine)
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBart <- 0; MBven <-0;  MLN <-0; MBK <-0; MKF <-0; MKT <-0; MFil <-0; Murine <-0; MBL <-0
    MLF <-0; MLT <-0; Mbile <-0; MBG <-0; MGF <-0; MGT <-0
    MGL <-0; Mfeces <-0; MBM <-0; MMF <-0; MMT <-0; MBA <-0; MAF <-0
    MAT <-0; MBR <-0; MRF <-0; MRT <-0; MLuAF<- 0;MBLu <- 0; MLuF <- 0;MLuT <- 0 
    
    return(c('MBart'=MBart, 'MBven'=MBven, 'MLN'=MLN, 'MBK'=MBK, 'MKF'=MKF, 'MKT'=MKT,
             'MFil'=MFil, 'Murine'=Murine, 'MBL'=MBL, 'MLF'=MLF, 'MLT'=MLT, 'Mbile'=Mbile,
             'MBG'=MBG, 'MGF'=MGF, 'MGT'=MGT, 'MGL'=MGL, 'Mfeces'=Mfeces, 'MBM'=MBM, 
             'MMF'=MMF, 'MMT'=MMT, 'MBA'=MBA, 'MAF'=MAF, 'MAT'=MAT, 'MBR'=MBR, 'MRF'=MRF,
             'MRT'=MRT,'MBLu'=MBLu, 'MLuF'=MLuF,'MLuT'=MLuT, 'MLuAF' = MLuAF  
    ))
    
    
  })
}


create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "iv"){
        events <- list(data = rbind(data.frame(var = c("MBven"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("MGL"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }
    }
    return(events)
  })
}


#  absolute average fold error
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

#  average fold error
AFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- (log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

#############################
#--------------------------
# Objective function
#-------------------------
obj.func <- function(x, dataset){
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  
  # Set up simulations for the first case, i.e. kudo (2007) high dose, tissues
  BW <- 0.290  # body weight (kg)
  admin.dose_per_g <- 16.56 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0.01#0.5/60 # time when doses are administered, in hours
  admin.type <- "iv"

  estimated_params <- exp(x)
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,2,0.01)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  exp_data <- dataset$df1 # retrieve data of Kudo et al. 2007 high dose
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kudo_high <- solution[solution$time %in% unique(exp_data$time), c("Cblood",
                                                  "Cliver","Ckidney","Cgut",
                                                  "Ccarcass", 'Clungs' )]
  preds_kudo_high <- preds_kudo_high /1000 #convert ug/kg to ug/g
  # Estimate gut concentration from stomach and intestines
  V_stomach <- 0.0046*BW*1000
  V_int <- 0.0124*BW*1000
  C_stomach <- exp_data[exp_data$Tissue == "Stomach", "concentration"] 
  C_int <- exp_data[exp_data$Tissue == "Intestine", "concentration"]   
  exp_gut <- (C_stomach*V_stomach+C_int*V_int)/(V_stomach+V_int)
  
  obs_kudo_high <- c(exp_data[exp_data$Tissue == "Blood", "concentration"],
           exp_data[exp_data$Tissue == "Liver", "concentration"],
           exp_data[exp_data$Tissue == "Kidney", "concentration"],
           exp_gut,
           exp_data[exp_data$Tissue == "Carcass", "concentration"],
           exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  # Set up simulations for the second case, i.e. kudo (2007) low dose, tissues
  #preds_kudo_low
  #obs_kudo_low
  #obs <- c(obs_kudo_high, obs_kudo_low)
  #preds <- c(preds_kudo_high, preds_kudo_low)
  
  # Initiate a variable to save the score of goodness of fit for each variable
  score <- AAFE(predictions = preds, observations = obs)
 
  return(score)
}

################################################################################
setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")
# Read data
kudo_high_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_high_kudo_2007.xlsx")
dataset <- list("df1" = kudo_high_dose)

#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX",#"NLOPT_LN_NEWUOA","NLOPT_LN_SBPLX"
              "xtol_rel" = 1e-07,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 500, 
              "print_level" = 1)
# Create initial conditions (zero initialisation)
N_pars <- 6# Number of parameters to be fitted
fit <- log(rep(1,N_pars))

# Run the optimization algorithmm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                               eval_f = obj.func,
                               lb	= rep(log(0.01), N_pars),
                               ub = rep(log(100), N_pars),
                               opts = opts,
                               dataset = dataset)

estimated_params <- exp(optimizer$solution)

# Set up simulations for the first case, i.e. kudo (2007) high dose, tissues
BW <- 0.290  # body weight (kg)
admin.dose_per_g <- 16.56 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0.01#0.5/60 # time when doses are administered, in hours
admin.type <- "iv"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,2,0.01)
 solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                     y = inits, parms = params,events = events,
                                     method="lsodes",rtol = 1e-05, atol = 1e-05))

 preds <- solution[solution$time %in% unique(dataset$df1$Time_hours), c("Cblood",
                                                               "Cliver","Ckidney","Cgut",
                                                               "Ccarcass", 'Clungs' )]
 preds <- preds /1000 #convert ug/kg to ug/g
 print(preds)
 # ######################################################################################
# 
# 
# # Plot with ggplot2
# library (ggplot2)
# compartment <- c('Cblood')
# color_codes <- scales::hue_pal()(length(compartment))
# 
# plot <- ggplot(data = solution)+
#   geom_line( aes(x = time, y = Cblood, color='Cblood'), size=1.3)+
#   scale_y_log10()+       
#   scale_x_continuous()+
#  
#   labs(title = 'Predicted values',
#        y = 'Concentration (ng/g)', x = "Time (hours)")+
#   theme(plot.title = element_text(hjust = 0.5,size=30),
#         axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
#         axis.text.y=element_text(size=22),
#         axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
#         axis.text.x=element_text(size=22),
#         legend.title=element_text(hjust = 0.5,size=25),
#         legend.text=element_text(size=22),
#         panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
#   
#   scale_color_manual("Compartment", values=color_codes)+
#   theme(legend.key.size = unit(1.5, 'cm'),
#         legend.title = element_text(size=14),
#         legend.text = element_text(size=14),
#         axis.text = element_text(size = 14))
# print(plot)
