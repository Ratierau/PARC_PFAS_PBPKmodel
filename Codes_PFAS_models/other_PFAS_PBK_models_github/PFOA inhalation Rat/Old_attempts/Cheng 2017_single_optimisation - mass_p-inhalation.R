
#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d

library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    # BW in kg
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
    #assuming the density of tissue is 1 g/mL.
    # Estimated parameters
    if (sex == "M"){
      RAFOatp_k <- 1.090359e-19
      RAFOat1 <- 4.163065e+01
      RAFOat3 <- 7.840798e+00
      #RAFOat3 <-  estimated_params[3]
      RAFOatp_l <- 3.556025e-18
      RAFNtcp <- 1.504512e+14 
      #RAFNtcp <- estimated_params[5]
    }else if(sex == "F"){
      RAFOatp_k <- 5.283413e-18
      RAFOat1 <- 2.967585e+07 
      RAFOat3 <- 3.059358e+04
     # RAFOat3 <-  estimated_params[8]
      RAFOatp_l <- 4.563770e+11 
      RAFNtcp <- RAFOatp_l
     # RAFNtcp <- estimated_params[10]
    }
    
    #permeabilities correction factor
    CF_Peff <- 1 
    # Absorption rate per area
    kabs <- 1
    # Bile correction factor
    #bile_correction_factor <- 1
    
    kabsUA <- estimated_params[1]
    kCLEal <- estimated_params[2]
    kCLEua <- estimated_params[3]
    kLuAF <- estimated_params[4]
    
    # urine correction factor
    CFurine <- 1#estimated_params[9] 
    # Correction factor for female GFR
    CF_GFR_female <-1# estimated_params[11]

    
    #units conversion from Cheng 2017R, time-> h, PFOA mass->ng, tissues mass-> g
    
    Hct <- 0.41 #hematocrit for rats, https://doi.org/10.1080/13685538.2017.1350156 mean value for both males and females
    
    #======Table S1=======#    
    
    #Blood
    VB <- (0.06*(BW*1000)+0.77)/1000 # Total blood volume (l), from Lee & Blaufox (1985), https://pubmed.ncbi.nlm.nih.gov/3965655/
    PVB <- 54e-3 #13.5 mL/244 g=0.055 mL/g~55e-3 mL/g (kg=L), Brown et al. 1997
    PVplasma <- 31.2e-3 
    Vplasma <- PVplasma * BW #plasma volume kg=L
    Vven <- BW*11.3/250 	#volume of venous plasma (L); from doi:10.1007/bf02353860
    Vart <- BW*5.6/250	#volume of arterial plasma (L); from doi:10.1007/bf02353860
    
    #Kidney
    PVK <- 7.3e-3 #Brown et al. 1997
    VK <- PVK * BW #kidney volume kg=L
    PVKB <- 0.16 #Brown et al. 1997
    VKB <- PVKB * PVK * BW #kidney blood volume kg=L
    PVKF <- 0.13 #Larson et al., 1984
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume kg=L
    VKT <- VK - VKF - VKB #kidney tissue volume kg=L
    VFil <- 0.25/1000 #0.25 mL renal filtrate volume Cheng et al., 2017 (from Arthur, 1986; Bonvalet, 1981)
    
    #Liver
    PVL <- 3.66e-2 #Brown et al. 1997
    VL <- PVL * BW #liver volume kg=L
    PVLB <- 0.21  #Brown et al. 1997
    VLB <- PVLB * PVL * BW #liver blood volume kg=L
    PVLF <- 0.049  #Blouin et al. 1977
    VLF <- PVLF * PVL* BW #liver interstitial fluid volume kg=L
    VLT <- VL - VLF #liver tissue volume kg=L
    # PVbile <- 0.004 #Blouin et al. 1977
    # Vbile <- 0.001 #PVbile * VL #bile volume kg=L
   
    #Intestine (small and large)
    PVIN <- 2.24e-2 #Brown et al. 1997, p 416, Table 5: 1.4+0.84
    VIN <- PVIN * BW #intestine volume kg=L
    # In the following we add the plasma and blood cell volumes of the small and large intestine from Shah and betts, (2012)
    PVINB <- ((0.0795+0.0651)+(0.0458+0.0375))/280 # From https://doi.org/10.2170/jjphysiol.33.1019
    VINB <- PVINB * VIN #intestine  blood volume kg=L
    PVINF <- (0.867+0.5)/280 # From https://doi.org/10.2170/jjphysiol.33.1019
    VINF <- PVINF * PVIN * BW #intestine interstitial fluid volume kg=L
    VINT <- VIN - VINF #intestine tissue volume kg=L
    
    #Stomach
    PVST <- 0.46e-2 #Brown et al. 1997, p 416, Table 5
    VST <- PVST * BW #stomach volume kg=L
    #In the next calculations, Nose et al. (1983) provide plasma volume per organ and we convert to blood by taking into account the hematocrit 
    #VSTB <- (0.089/(1-Hct))*((0.98*(1000*BW/100))/1000) - VINB # GI (stomach + small + large) blood volume in ml/g of dry tissue, dry tissue mass = 0.98 g/100g BW, https://doi.org/10.2170/jjphysiol.33.1019
    #VSTF <- 0.62*((0.98*(1000*BW/100))/1000) - VINF # GI (stomach + small + large) IF volume in ml/g of dry tissue, dry tissue mass = 0.98 g/100g BW, https://doi.org/10.2170/jjphysiol.33.1019
    PVSTB <- 0.032 #from pkSim
    VSTB <- PVSTB * VST
    PVSTF <-  0.10 # from pkSim
    VSTF <- PVSTF * VST
    VSTT <- VST - VSTF #stomach tissue volume kg=L
    
    #Stomach and intestine lumen
    PVSTL <- 3.4/175 # Connell et al., 2008, https://doi.org/10.1211/jpp.60.1.0008
    VSTL <- PVSTL * BW #stomach lumen volume kg=L
    PVINL <- (0.894+0.792+0.678+0.598+0.442)/230 #Funai et al., 2023 https://doi.org/10.1038/s41598-023-44742-y --> Figure 3C
    VINL <- PVINL * BW #intestine lumen volume kg=L
     
    #Muscle
    PVM <- 40.43e-2 #Brown et al. 1997
    VM <- PVM * BW #muscle volume kg=L
    PVMB <- 0.04 #Brown et al. 1997
    VMB <- PVMB * PVM * BW #muscle blood volume kg=L
    PVMF <- 0.054 #Ng, 2013
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume kg=L
    VMT <- VM - VMF #muscle tissue volume kg=L
    
    #Adipose
    PVA <- 7e-2 #Brown et al. 1997
    VA <- PVA * BW #adipose volume kg=L
    PVAB <- 0.02 #Brown et al. 1997
    VAB <- PVAB * PVA * BW #% adipose blood volume kg=L
    PVAF <- 0.174 #Ng, 2013
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume kg=L
    VAT <- VA - VAF #adipose tissue volume kg=L
    
    #Lymph nodes
    PVLN <- 1.15/280 #Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLN <- PVLN*BW#lymph fluid volume
    
    #Upper airways
    PVUA <- 257e-3/288 #mm^3/g, for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    VUA <- PVUA * BW  #total volume of nasal cavity 
    
    #respiratory area of nasal cavity
    RA_area = 623 #mm^2
    RA_capillary_density = 362 #capillaries/mm^2
    capillary_area = RA_area/(RA_area*RA_capillary_density) #mm^2, area of one capillary
    
    #Lung
    PVLu <- 0.48e-2 #Brown et al. 1997, p 418-19 mean of values for male rat
    VLu <- PVLu * BW
    PVLuB <- 0.09*VLu #Brown et al. 1997, p 459 --> capillary blood occupied 9% of the lung volume
    VLuB <- PVLuB * PVLu * BW #volume of the blood of lung kg=L
    PVLuF <- 0.263/280 #0.263 ml, Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLuF <- PVLuF * PVLu * BW #lung interstitial fluid volume
    PVLuAF <- 0.4/275 #0.4 mL Leslie et al, 1989 https://doi.org/10.1164/ajrccm/139.2.360 --> Watkins & Rannels 1979 https://doi.org/10.1152/jappl.1979.47.2.325  
    VLuAF <- PVLuAF * PVLu * BW #lung alveolar lining fluid volume kg=LL
    PVLuT <- VLu - VLuF - VLuAF 
    VLuT <- PVLuT * PVLu* BW #lung tissue volume kg=L
    
    IVR=0
    EC=0
    Dua=0
    
    
    #Spleen
    PVSP <- 0.2e-2  #Brown et al. 1997, p 416, Table 5
    VSP <- PVSP * BW
    PVSPB <- 0.22 #Brown et al. 1997, p 458, Table 30
    VSPB <- PVSPB * PVSP * BW #volume of the blood of spleen kg=L
    PVSPF <- 0.554/280 #ml/g -> Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VSPF <- PVSPF * PVSP * BW #spleen interstitial fluid volume kg=L
    VSPT <- VSP - VSPF #spleen tissue volume kg=L
    
    #Heart
    PVH <- 0.33e-2  #Brown et al. 1997, p 416, Table 5
    VH <- PVH * BW
    PVHB <- 0.26 #Brown et al. 1997, p 458, Table 30
    VHB <- PVHB * PVH * BW #volume of the blood of heart kg=L
    PVHF <- PVH * 0.22 #https://doi.org/10.1016/S0022-2828(87)80559-X --> IS 22% of dry tissue
    VHF <- PVHF * PVH * BW #heart interstitial fluid volume kg=L
    VHT <- VH - VHF #heart tissue volume kg=L
    
    #Brain
    PVBr <- 0.57e-2  #Brown et al. 1997, p 416, Table 5
    VBr <- PVBr * BW
    PVBrB <- 0.03 #Brown et al. 1997, p 458, Table 30
    VBrB <- PVBrB * PVBr * BW #volume of the blood of brain kg=L
    PVBrF <- 17.5 * VBrB #https://doi.org/10.1016/j.pneurobio.2015.12.007 --> The ISS occupies 15% to 20% of the total brain volume 
    VBrF <- PVBrF * PVBr * BW #brain interstitial fluid volume kg=L
    VBrT <- VBr - VBrF #brain tissue volume kg=L
    
    #Testis
    PVT <- 1.24/239.03 #Table 01 https://doi.org/10.1590/S1516-89132012000100013
    VT <- PVT * BW
    PVTB <- 4.75e-2/1.49/243 #Figure 2 https://doi.org/10.1007/BF00410278 --> blood weight/g of testis/BW
    VTB <-PVTB * PVT * BW #volume of the blood of testis kg=L
    PVTF <- 0.16/239.03 #Table 02 https://doi.org/10.1590/S1516-89132012000100013
    VTF <- PVTF * PVT * BW #testis interstitial fluid volume kg=L
    VTT <- VT - VTF #testis tissue volume kg=L

    #Skin
    PVSK <- 19.03/100 #Brown et al. 1997, p 416, Table 5
    VSK <- PVSK * BW
    PVSKB <- 0.02 #Brown et al. 1997, p 458, Table 3
    VSKB <-PVSKB * PVSK * BW #volume of the blood of testis kg=L
    PVSKF <- 40*VSK/100/0.225 # https://doi.org/10.1111/j.1748-1716.1981.tb06901.x 40 mL/100 g tissue, BW = 200-250 g
    VSKF <- PVSKF * PVSK * BW #testis interstitial fluid volume kg=L
    VSKT <- VSK - VSKF #testis tissue volume kg=L
    
    #RoB
    PVR <- 1 - PVB - PVK - PVL - PVM - PVA - PVSP - PVH - PVBr - PVT - PVST - PVIN - PVSK
    #PVR <- 1 - PVB - PVK - PVL - PVM - PVA - PVSP - PVH - PVBr - PVT - PVST - PVIN
    VR <- PVR * BW #volume of the rest of the body kg=LL
    PVRB <- 0.036 
    VRB <- PVRB * PVR * BW #volume of the blood of the rest of body kg=L
    PVRF <- 0.18  
    VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body kg=L
    VRT <- VR - VRF #tissue volume of the rest of body kg=L
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight
    #or weight of corresponding tissue (PAi, m^2/g) and surface area (m^2)
    
    PAK <- 350e-4 #cm2/g  --> *e-4 --> m^2/g  https://doi.org/10.1111/j.1748-1716.1963.tb02652.x
    AK <- PAK * VK * 1e3 #kidney capillary surface area (m^2)
    PAKG <- 68.90e-4 
    AKG <- PAKG * VK * 1e3 #the surface area of glomerular capillary (m^2)
    PAL <- 250e-4 
    AL <- PAL * VL * 1e3 #liver capillary surface area (m^2)
    
    PAST <- 100e-4 #Cheng et al., 2017 value for gut
    AST <- PAST * VST * 1e3 #stomach capillary surface area (m^2)
    PASTL<- 100e-4 #Cheng et al., 2017 value for gut
    ASTL<- PASTL * VSTL #stomach lumen capillary surface area (m^2)
    PAIN <- 100e-4 #Cheng et al., 2017 value for gut
    AIN <- PAIN * VIN * 1e3 #intestine capillary surface area (m^2)
    # PAINL<- 100e-4 #Cheng et al., 2017 value for gut
    # AINL<- PAINL * VINL * 1e3 #intestine lumen capillary surface area (m^2)
    
    
    #Calculations of rat intestinal lumen surface area based on Kothari et al. (2020),https://doi.org/10.1002/btm2.10146
    # Lengths
    L_duodenum <- 9.6# cm
    L_jejunum <- 26 #cm
    L_ileum <- 34.4# cm
    # Inner diameters
    d_duodenum <- 2.21 #cm
    d_jejunum <- 2.56 #cm
    d_ileum <- 3.36 #cm
    
    SA_ref <- 2*pi*(d_duodenum/2)*L_duodenum + 2*pi*(d_jejunum/2)*L_jejunum + 2*pi*(d_ileum/2)*L_ileum #cm^2
    #Allometric scaling
    SA <- SA_ref* (BW/0.19)^(2/3) # O.19 kg is the BW of rats used in Kothari et al. (2020)
    n <- 5 #enlargement factor of apical membrane of proximal tubule or enterocytes
    AINL <- n * SA * 1e-4 #m^2
    
    PAM <- 70e-4 #m2/g tissue
    AM <- PAM * VM *1e03 #muscle surface area (m^2)
    PAA <- 70e-4 #m2/g tissue
    AA <- PAA * VA *1e03 #adipose surface area (m^2)
    PAR <- 100e-4#m2/g tissue
    AR <- PAR * VR *1e03 #surface area of rest of body (m^2)
    
    
    PLu <- 250e-4  #m2/g tissue      #https://doi.org/10.1111/j.1748-1716.1963.tb02652.x
    ALu <- PLu* VLu*1e03 #lung surface area (m^2)
    PSP <- 70e-4 #m2/g tissue
    ASP <- PSP* VSP *1e03 #spleen surface area (m^2), same as muscle #assumption
    PH <-  67.8# mm^2/mm^3 <--> m2/L of tissue, from group III in https://doi.org/10.1016/S0022-2828(87)80533-3
    AH <- PH* VH #heart surface area (m^2)
    PBr <- 240e-4#m2/g tissue     #https://doi.org/10.1111/j.1748-1716.1963.tb02652.x 
    ABr <- PBr* VBr *1e03 #brain surface area (m^2)
    PT <- 70e-4#m2/g tissue
    AT <- PT* VT *1e03 #testis surface area (m^2), same as muscle #assumption
    PSK <- 70e-4#m2/g tissue
    ASK <- PSK* VSK*1e03 #skin surface area (m^2), same as muscle #assumption
    
    #Organ surface area
    
    PAUA <- ((1343.5+76.8)*1e-06)/((156+257)*1e-03) #m2/g tissue, total surface area 1343.5 mm^2 for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    #plus nasopharynx area 76.8 mm^2 calculated by Ménache et al., 1997, https://doi.org/10.1080/00984109708984003
    AUA <- PAUA * VUA*1e03 #total surface area of upper airways (m^2) !think of separate nose in different areas of absorption
    PALF <- 1.1529e-8/173e-10 #11,529 μm2/173e3 μm3--> 1.1529e-8 m2/g tissue, Mercer et al.,  https://doi.org/10.1152/jappl.1987.62.4.1480
    ALF <- PALF* VLF*1e03 #alveolar lining fluid area (m^2)
      
    
    
    #Effective permeability (Peff, in m/h) for blood (B), liver(L), kidney(K),
    #stomach(ST),intestine (IN), adipose(A), muscle(M), spleen (SP), heart (H), 
    #brain (Br), testis (T), rest of body(R)

    #PeffB <- 4.98e-8*3600*CF_Peff1
    #PeffK <- 4.38e-8*3600*CF_Peff1
    #PeffL <- 5.15e-8*3600*CF_Peff1
    PeffK <- 5e-8*3600*CF_Peff
    PeffL <- 5e-8*3600*CF_Peff
    #PeffG <- 2.65e-8*3600
    PeffST <- 5e-8*3600*CF_Peff #assumption
    PeffIN <- 5e-8*3600*CF_Peff#assumption
    PeffA <- 5e-8*3600*CF_Peff
    PeffM <- 5e-8*3600*CF_Peff
    PeffR <- 5e-8*3600*CF_Peff
    PeffLu <- 5e-8*3600*CF_Peff #assumption
    PeffSP <- 5e-8*3600*CF_Peff #assumption
    PeffH <- 5e-8*3600*CF_Peff #assumption
    PeffBr <- 5e-8*3600*CF_Peff #assumption
    PeffT <- 5e-8*3600*CF_Peff #assumption
    PeffSK <- 5e-8*3600*CF_Peff #assumption
    
    
    #Blood flow rates (QBi, in L/h) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac L/h), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 0.235 * (BW^0.75) *60 #L/min->*60-> L/h
    PQBK <- 14.1/100 #Brown et al. 1997, p 438, Table 23
    QBK <- PQBK * Qcardiac #L/h
    PQBL <- 2.1/100 
    QBL <- PQBL * Qcardiac #L/h Brown et al. 1997, p 438, Table 23
    PQBST <- 0.16/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBST = (8/1000)/Qcard *100
    QBST <- PQBST * Qcardiac #L/h
    PQBIN <- 9/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBIN = (451/1000)/Qcard *100
    QBIN <- PQBIN * Qcardiac #L/h
   
    QGE<- 0.54/BW^0.25 #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013
    #QGE<- VSTL * 2.03 #L/h https://doi.org/10.1124/dmd.118.085902, Table 5   
    
    PQBM <- 27.8/100 #Brown et al. 1997, p 438, Table 23
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 7/100 #Brown et al. 1997, p 438, Table 23
    QBA <- PQBA * Qcardiac #L/h
    
    PQBLu <- 1 
    QBLu <- PQBLu * Qcardiac #L/h
    QBUA <- PQBUA * Qcardiac #L/h
    
    
    PQBSP <- 0.75/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBSP = (37.5/1000)/Qcard *100
    QBSP <- PQBSP * Qcardiac #L/h
    PQBH <- 5.1/100 #Brown et al. 1997, p 438, Table 23
    QBH <- PQBH * Qcardiac #L/h
    PQBBr <- 2.0/100 #Brown et al. 1997, p 438, Table 23
    QBBr <- PQBBr * Qcardiac #L/h 
    PQBT <- 0.28/100 #https://doi.org/10.1152/ajpregu.1987.253.2.R228 Qcard=0.235*(0.335^0.75)*60 (L/h) and PQBT = (0.295*60/1000)/Qcard *100
    QBT <- PQBT * Qcardiac #L/h
    PQBSK <- 1.35/100 #https://doi.org/10.1111/1523-1747.ep12277181 Qcard=0.235*(0.2^0.75)*60 (L/h) and PQBSK = (0.95*60/1000)/Qcard *100
    QBSK <- PQBSK * Qcardiac #L/h
    
    QBLtot <- QBL+QBSP+QBIN+QBST
    
    PQBR = 1 - PQBK - PQBST - PQBIN - PQBL - PQBM - PQBA - PQBH - PQBSK - PQBSP - PQBT - PQBBr
    QBR <- PQBR * Qcardiac #L/h
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    
    #Qfeces <- 5.63/1000/24 #mL water/d --> L/h
    PQbile <- 90/1000/24 # [mL/d/kg BW]/1000/24 --> L/h/kg BW
    Qbile <- PQbile * BW #L/h
    #PQurine <- 200/1000/24 #mL/d/kg BW --> L/h/kg
    #Qurine <- PQurine * BW #L/h
    Qfeces <- (8.18/0.21)/24*BW #g/kg BW, based on Cui et al.(2010)
    feces_density <- 1.29 #g/cm^3 --> g/mL from Lupton 1986, Fig 1. Fiber free control diet, https://doi.org/10.1093/jn/116.1.164

    if (sex == "M"){
      PQGFR <- 62.1  #L/h/kg   Corley et al., 2005 https://doi.org/10.1093/toxsci/kfi119, GFRC --> scaled fraction of kidney weight
      QGFR <- PQGFR * VK #L/h
      Qurine <- (60/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, 10.1002/nau.1006
      }else if(sex == "F"){
      PQGFR <- 41.04  #L/h/kg  Corley et al., 2005 https://doi.org/10.1093/toxsci/kfi119, GFRC --> scaled fraction of kidney weight
      QGFR <- PQGFR * VK #L/h
      Qurine <- (85/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, 10.1002/nau.1006  
    }
    
    # if (sex == "M"){
    #   PQGFR <- 1.08e-3*60*10/1000 #L/h/kg   1.08 mL/min/100 g BW --> L/h/ kg Bw Sadick et al., 2011, https://doi.org/10.1093/ndt/gfr148
    #   QGFR <- PQGFR * BW#L/h
    # }else if(sex == "F"){
    #   PQGFR <- 1.08e-3*60*10/1000 #L/h/kg   1.08 mL/min/100 g BW --> L/h/ kg Bw Sadick et al., 2011, https://doi.org/10.1093/ndt/gfr148
    #   QGFR <- PQGFR  * BW * CF_GFR_female #L/h
    #   Qurine <- (85/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, 10.1002/nau.1006  
    # }
    
    #======Table S2=======#
    #Albumin concentration in blood and interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)

    CalbB <- 486*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbKF <- 243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbLF <- 243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    
    CalbGF <- 146*1e-6 #[umol/L]*1e-6 -->(mol/m3), from Cheng et al. (2017)
    
    CalbSTF <- 146*1e-6 # [umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
    CalbINF <- 146*1e-6 #[umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
    
    CalbMF <- 146*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbAF <- 73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbRF <- 73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    
    CalbLN <- CalbGF #assumption based on https://doi.org/10.1016/j.jconrel.2020.07.046
    CalbLuF <- CalbGF #assumption 
    CalbLuAF <- 10/100 * CalbB #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
    
    MW_albumin <- 66500#g/mol
    
    CalbSPF <- 243e-3 #[umol/L]*1e-6 -->(mol/L), same as liver (assumption)
    CalbTF <- 41/MW_albumin #mg/mL-->  (mol/L) from https://doi.org/10.1210/endo-116-5-1983 --> 41 mg/mL, MW=65 kg/mol
    CalbHF <- 65/MW_albumin##mg/mL--> (mol/L) https://doi.org/10.1007/s12291-010-0042-x --> 6.5 g/100 g tissue, MW=65 kg/mol
    CalbBrF <- 8e-2/MW_albumin  ##mg/mL--> (mol/L) https://doi.org/10.1016/0014-4886(90)90158-O --> 0.08 g/L, MW=65 kg/mol 
    CalbSKF <- 21/MW_albumin ##mg/mL-->  (mol/L) https://doi.org/10.1111/j.1748-1716.1973.tb05464.x -->Table 2: 2.1 g/100 mL
    
    
    #Alpha2mu-globulin concentration in kidney tissue (mol/L)
    if (sex == "M"){
      Ca2uKT <- 110e-3 #mol/m3
    }else if(sex == "F"){
      Ca2uKT <- 0 #mol/m3
    }
    
    
    #FABP concentration in kidney and liver tissue (mol/m^3)
    CLfabpKT <- 2.65*1e-6  #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CLfabpLT <- 133*1e-6  #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    
    
    #======Table S2=======#
    #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
    #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details

    #Ka <-  24.18 #3.1*7.8 m3/mol multiplying by number of binding sites (Cheng et al. 2021)
    Ka <-  1e05/1e03 #3.1*7.8 m3/mol multiplying by number of binding sites (Cheng et al. 2021)
    
    KLfabp <- 135  #45.0*3 geo_mean of 3 binding affinity, may try normal mean (Cheng et al. 2021)
    Ka2u <- 0.5 #m3/mol

    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #passive diffusion rates
    
    ClLFT_unscaled= 67.8#uL/min/10^6 cells, Han et al. 2008
    ClKFT_unscaled= 17.5 #uL/min/mg protein, Yang et al. 2010
    
    ClGFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    
    ClSTFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    
    ClMFT_unscaled= 18.1 #same as ClGFT
    ClAFT_unscaled= 18.1 #same as ClGFT
    ClRFT_unscaled= 18.1 #same as ClGFT
    ClLuFT_unscaled= 18.1 #same as ClGFT
    
    ClSPFT_unscaled= 18.1 #same as ClGFT
    ClHFT_unscaled= 18.1 #same as ClGFT
    ClBrFT_unscaled= 18.1 #same as ClGFT
    ClTFT_unscaled= 18.1 #same as ClGFT
    ClSKFT_unscaled= 18.1 #same as ClGFT
    
    
    #For all CMTs
    MW = 414.07 #g/mol, PFOA molecular weight
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
    
    n <- 5 #enlargement factor of apical membrane of proximal tubule
    kFKT <- PeffK * AK * n
    
    #Oatp kidney
    VmK_Oatp_in_vitro <- 9.3 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oatp_scaled <- 60*VmK_Oatp_in_vitro*MW*kidney_protein_total/1000  #physiologically scaled to in vivo, ug/L/h
    VmK_Oatp <- VmK_Oatp_scaled*RAFOatp_k #in vivo value, in  ug/L/h
    KmK_Oatp=  126.4*MW#umol/L (Weaver et al. 2010) --> ug/L
    
    #oat1 kidney
    VmK_Oat1_in_vitro= 2.6 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oat1_scaled = 60*VmK_Oat1_in_vitro*MW*kidney_protein_total/1000  #physiologically scaled to in vivo, ug/L/h
    VmK_Oat1= VmK_Oat1_scaled*RAFOat1 #in vivo value, in  ug/L/h
    KmK_Oat1= 43.2 * MW #umol/L (Weaver et al. 2010) --> ug/L
    
    #oat3 kidney
    VmK_Oat3_in_vitro= 3.8 #nmol/mg protein/min  (Weaver et al. 2010)
    VmK_Oat3_scaled = 60*VmK_Oat3_in_vitro*MW*kidney_protein_total/1000  #physiologically scaled to in vivo, ug/L/h
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
    
    #kbileLT <- PeffL * AL * bile_correction_factor
    
    # oatp-liver
    VmL_Oatp_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
    VmL_Oatp_scaled = 60*VmL_Oatp_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000  #physiologically scaled to in vivo, ug/L/h
    VmL_Oatp = VmL_Oatp_scaled*RAFOatp_l #in vivo value, in  ug/L/h
    KmL_Oatp = KmK_Oatp #same as kidney
    
    #Ntcp liver
    VmL_Ntcp_in_vitro= 3#nmol/mg protein/min  (Weaver et al. 2010)
    VmL_Ntcp_scaled = 60*VmL_Ntcp_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000  #physiologically scaled to in vivo, ug/L/h
    VmL_Ntcp = VmL_Ntcp_scaled*RAFNtcp #in vivo value, in  ug/L/h
    KmL_Ntcp= 20 * MW #umol/L, Ruggiero et al. 2021 --> ug/L
    
    #Muscle
    muscle_cells= NA
    Cmedium_M = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    muscle_protein <- 158.45 #mg/g muscle (protein data from Cheek et al.,1971 and muscle mass from Caster et al.,1956)
    muscle_protein_tot <- muscle_protein * (VMT*1000)
    ClMFT <- ClMFT_unscaled * muscle_protein_tot#uL/min for the whole muscle comp
    kMFMT <-  (60*ClMFT)/1e06 #L/h
    
    #Intestine
    intestine_cells = 93.1 * 2.365e6 # https://doi.org/10.1152/ajpgi.00290.2013 cells per crypt Figure 3B, number of crypts: Figure 1, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1271340/ (mean value of rat 300 and 400 g)
    Cmedium_G = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    intestine_protein <- muscle_protein#NA#5034
    intestine_protein_total <- intestine_protein*(1000*VINT)
    ClINFT <- ClINFT_unscaled *intestine_protein_total#uL/min for the whole gut compartment
    kINFINT <-  (60*ClINFT)/1e06 #L/h
    kabIN <- kabs * AINL # AINL in m^2/h
    
    #Stomach
    stomach_cells = NA
    Cmedium_G = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    stomach_protein <- muscle_protein#NA#5034
    stomach_protein_total <- stomach_protein*(1000*VSTT)
    ClSTFT <- ClSTFT_unscaled *stomach_protein_total #uL/min for the whole gut compartment
    kSTFSTT <-  (60*ClSTFT)/1e06 #L/h
    # For identifiability reasons we assume that absorption is realised only through the intestines
    #kabST <- kabs * ASTL
    kabST <- (kabs* ASTL)*1000
    
    #Adipose
    adipose_cells = NA
    Cmedium_A = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    adipose_protein <- muscle_protein#NA#5034
    adipose_protein_total <- adipose_protein * (1000*VAT)
    ClAFT <- ClAFT_unscaled * adipose_protein_total#uL/min
    kAFAT <-  (60*ClAFT)/1e06 #L/h
    
    #Rest of body
    RoB_cells = NA
    Cmedium_R = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    RoB_protein <- muscle_protein#18456
    RoB_protein_total <- RoB_protein * (1000* VRT)
    ClRFT <- ClRFT_unscaled * RoB_protein_total#uL/min
    kRFRT <-  (60*ClRFT)/1e06 #L/h
    
    #Lung
    Lung_protein <- muscle_protein#18456
    Lung_protein_total <- Lung_protein * (1000* VRT)
    ClLuFT <- ClLuFT_unscaled * Lung_protein_total#uL/min
    kLuTLuF <-  (60*ClLuFT)/1e06 #L/h
    kUAB <- kabsUA * RA_area #transport rate from upper airways to blood
    CLEal <- kCLEal * ALF #clearance rate from alveolar lining fluid, ALf in m^2
    CLEua <- kCLEua * AUA #clearance rate rate from upper airways to stomach, AUA in m^2
    kLuAFLuT <- kLuAF * ALF #transport rate from alveolar lining fluid to lung tissue, ALF in m^2
    
    
    #Spleen
    spleen_cells = 1.76e8 #cells/g tissue https://doi.org/10.1371/journal.pone.0059602 --> Figure 5C
    Cmedium_SP = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    spleen_protein <- muscle_protein#NA#5034
    spleen_protein_total <- spleen_protein * (1000*VSPT)
    ClSPFT <- ClSPFT_unscaled * spleen_protein_total#uL/min
    kSPFSPT <-  (60*ClSPFT)/1e06 #L/h
    
    #Heart
    heart_cells = 3.3e8 #cells/g tissue  https://doi.org/10.1620/tjem.95.177 
    Cmedium_H = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    heart_protein <- muscle_protein#NA#5034
    heart_protein_total <- heart_protein * (1000*VHT)
    ClHFT <- ClHFT_unscaled * heart_protein_total#uL/min
    kHFHT <-  (60*ClHFT)/1e06 #L/h
    
    #Brain
    brain_cells = 3.3e8/1.8 #cells/g tissue https://doi.org/10.1523/JNEUROSCI.4526-04.2005 --> Table 1
    Cmedium_Br = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    brain_protein <- muscle_protein#NA#5034
    brain_protein_total <- brain_protein * (1000*VBrT)
    ClBrFT <- ClBrFT_unscaled * brain_protein_total#uL/min
    kBrFBrT <-  (60*ClBrFT)/1e06 #L/h
    
    #Testis
    testis_cells = 1.85e7+1.58e7 #Sertoli and Leydig cells/g tissue  https://doi.org/10.1007/BF00297504 --> Figure 2, LC and SC cells
    Cmedium_T = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    testis_protein <- muscle_protein#NA#5034
    testis_protein_total <- testis_protein * (1000*VTT)
    ClTFT <- ClTFT_unscaled * testis_protein_total#uL/min
    kTFTT <-  (60*ClTFT)/1e06 #L/h

    #Skin
    skin_cells = NA
    Cmedium_SK = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
    skin_protein <- muscle_protein#NA#5034
    skin_protein_total <- skin_protein * (1000*VSKT)
    ClSKFT <- ClSKFT_unscaled * skin_protein_total#uL/min
    kSKFSKT <-  (60*ClSKFT)/1e06 #L/h
    
    return(list('VB'=VB, 'Vplasma'=Vplasma, 'VK'=VK, 'VKB'=VKB, 
                'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'VL'=VL, 'VLB'=VLB, 'VLF'=VLF, 'VLT'=VLT, 
                'VM'=VM, 'VMB'=VMB, 'VMF'=VMF, 'VMT'=VMT, 'VA'=VA, 'VAB'=VAB, 
                'VAF'=VAF, 'VAT'=VAT, 'VR'=VR, 'VRB'=VRB, 
                'VRF'=VRF, 'VRT'=VRT, 'VLN' = VLN, 'VUA'=VUA, 'Vven' = Vven,
                'Vart' = Vart, 'VLu'=VLu, 'VLuB'=VLuB, 'VLuF'=VLuF,
                'VLuAF'=VLuAF, 'VLuT'=VLuT, 
                'IVR'=IVR, 'EC'=EC, 'Dua'=Dua,
                'VSP'=VSP, 'VSPB'=VSPB, 'VSPF'=VSPF, 'VSPT'=VSPT,
                'VH'=VH, 'VHB'=VHB, 'VHF'=VHF, 'VHT'=VHT,
                'VBr'=VBr, 'VBrB'=VBrB, 'VBrF'=VBrF, 'VBrT'=VBrT,
                'VT'=VT, 'VTB'=VTB, 'VTF'=VTF, 'VTT'=VTT,
                'VIN'=VIN, 'VINB'=VINB, 'VINF'=VINF, 'VINT'=VINT,
                'VST'=VST, 'VSTB'=VSTB, 'VSTF'=VSTF, 'VSTT'=VSTT,
                'VSTL'=VSTL, 'VINL'=VINL,
                'VSK'=VSK,'VSKB'=VSKB, 'VSKF'=VSKF, 'VSKT'=VSKT,
                
                'AK'=AK, 'AKG'=AKG, 'AL'=AL, 
                'AM'=AM, 'AA'=AA, 'AR'=AR, 'AUA'=AUA, 'ALu'= ALu, 
                'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
                'AIN'=AIN, 'AT'=AT,
                'ASK'= ASK, 'ALF'=ALF,
                
                #'PeffB'=PeffB, 
                'PeffK'=PeffK, 'PeffL'=PeffL, 
                'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
                'PeffSP'=PeffSP, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 'PeffST' = PeffST,
                'PeffIN'=PeffIN, 'PeffT'=PeffT,
                'PeffSK' = PeffSK,
                
                'Qcardiac'=Qcardiac, 'QBK'=QBK, 
                'QBL'=QBL, 'QBLtot'=QBLtot,
                'QBM'=QBM, 'QBA'=QBA,
                'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces, 'feces_density'=feces_density,
                'Qbile'=Qbile, 'QGFR'=QGFR,'Qurine'=Qurine,
                'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
                'QBIN'=QBIN, 'QGE'=QGE,
                'QBT'=QBT,
                'QBSK'=QBSK,
                
                'CalbB'= CalbB, 'CalbKF'=CalbKF, 'CalbLF'=CalbLF,
                'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF,'CalbLN' =CalbLN,
                'CalbLuF' =CalbLuF, 'CalbLuAF' =CalbLuAF,
                'CalbSPF' =CalbSPF, 'CalbTF' =CalbTF, 'CalbHF' =CalbHF,
                'CalbBrF' =CalbBrF, 'CalbSTF' =CalbSTF, 'CalbINF' =CalbINF,
                'CalbSKF' =CalbSKF, 
                
                'Ca2uKT'=Ca2uKT,'CLfabpKT'=CLfabpKT,'CLfabpLT'=CLfabpLT, 
                
                'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,
                
                'ClLFT'=ClLFT, 'ClKFT'=ClKFT,
                'ClMFT'=ClMFT,
                'ClAFT'=ClAFT, 'ClRFT'=ClRFT, 'ClLuFT'=ClLuFT, 'ClSTFT'=ClSTFT,
                'ClINFT'=ClINFT, 'ClSPFT'=ClSPFT, 'ClHFT'=ClHFT, 'ClBrFT'=ClBrFT, 
                'ClTFT'=ClTFT,
                'ClSKFT'=ClSKFT,
                
                
                'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp,
                'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 
                'KmK_Oatp'=KmK_Oatp,
                #'VmK_Osta'=VmK_Osta,'KmK_Osta'=KmK_Osta, 
                'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
                'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
                  
                
                'kKFKT'=kKFKT, 'kFKT'=kFKT,  
                'kLFLT'=kLFLT, 'kAFAT'=kAFAT, 
                'kRFRT'=kRFRT,
                'kabST'=kabST, 'kabIN'=kabIN,
                'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF,
                'kUAB'= kUAB, 'CLEal'=CLEal, 'CLEua'=CLEua, 'kLuAFLuT'=kLuAFLuT,
                'kSPFSPT' =kSPFSPT,'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
                'kBrFBrT' =kBrFBrT, 'kTFTT' =kTFTT,
                'kSKFSKT' =kSKFSKT,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, "capillary_area"=capillary_area
                ))
  
  })
}  

estimate_BF_TV <- function(sex){
  
  #Breathing frequency
  if (sex == "M"){
    PBF <- 168*60/0.219  #1/min -->  1/h/kg, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    BF <- PBF * BW  # 1/h
  }else if(sex == "F"){
    PBF <- 177*60/0.145  #1/min -->  1/h/kg, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    BF <- PBF * BW  # 1/h
  }
  
  #Tidal Volume
  if (sex == "M"){
    PTV <- 1.6/219   #mL/g --> mL=g, unitless,Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    TV <- PTV * BW  #L
  }else if(sex == "F"){
    PTV <- 1.2/145   #mL/g --> mL=g, unitless,Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    TV <- PTV * BW  #L
  }
  return(c("BF"=BF, "TV"=TV))
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
    CKB = MBK/VKB # blood concentration
    CKF = MKF/VKF #interstitial fluid concentration
    CKT = MKT/VKT # tissue concentration
    CFil = MFil/VFil# filtrate concentration
    
    #Liver
    CLB = MBL/VLB # blood concentration
    CLF = MLF/VLF #interstitial fluid concentration
    CLT = MLT/VLT # tissue concentration
    
    
    #Stomach
    CSTB = MBST/VSTB # blood concentration
    CSTF = MSTF/VSTF  #interstitial fluid concentration
    CSTT = MSTT/VSTT # tissue concentration
    
    #Intestine
    CINB = MBIN/VINB # blood concentration
    CINF = MINF/VINF  #interstitial fluid concentration
    CINT = MINT/VINT # tissue concentration
    
    #Stomach and Intestine lumens
    CSTL = MSTL/VSTL # Stomach Lumen concentration
    CINL = MINL/VINL # Intestine Lumen concentration
    
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
    
    #Upper airways
    CUA = MUA/VUA # upper airways concentration
    
    #Lung
    CLuB = MBLu/VLuB # blood concentration
    CLuF = MLuF/VLuF  #interstitial fluid concentration
    CLuT = MLuT/VLuT #lung tissue concentration
    CLuAF = MLuAF/VLuAF # alveolar lining fluid concentration
    
    #Spleen
    CSPB = MBSP/VSPB # blood concentration
    CSPF = MSPF/VSPF  #interstitial fluid concentration
    CSPT = MSPT/VSPT # tissue concentration
    
    #Heart
    CHB = MBH/VHB # blood concentration
    CHF = MHF/VHF  #interstitial fluid concentration
    CHT = MHT/VHT # tissue concentration
    
    #Brain
    CBrB = MBBr/VBrB # blood concentration
    CBrF = MBrF/VBrF  #interstitial fluid concentration
    CBrT = MBrT/VBrT # tissue concentration
    
    #Testis
    CTB = MBT/VTB # blood concentration
    CTF = MTF/VTF  #interstitial fluid concentration
    CTT = MTT/VTT # tissue concentration
    
    #Skin
    CSKB = MBSK/VSKB # blood concentration
    CSKF = MSKF/VSKF  #interstitial fluid concentration
    CSKT = MSKT/VSKT # tissue concentration
    
    #Cfree calculation using the expression of free fraction ff
    CBfart = CBart * 1.0 / (1.0 + CalbB * Ka)
    CBfven= CBven * 1.0 / (1.0 + CalbB * Ka)
    CLNf = CLN * 1.0 / (1.0 + CalbB * Ka)  # find albumin concentration in Lymph Nodes
    
    #Calculation of free concentrations in organ blood
    CKBf = CKB * 1.0 / (1.0 + CalbB * Ka)
    CLBf = CLB * 1.0 / (1.0 + CalbB * Ka)
    
    CSTBf = CSTB * 1.0 / (1.0 + CalbB * Ka)
    CINBf = CINB * 1.0 / (1.0 + CalbB * Ka)
    
    CMBf = CMB * 1.0 / (1.0 + CalbB * Ka)
    CABf = CAB * 1.0 / (1.0 + CalbB * Ka)
    CRBf = CRB * 1.0 / (1.0 + CalbB * Ka)
    CLuBf = CLuB * 1.0 / (1.0 + CalbB * Ka)
   
    CSPBf = CSPB * 1.0 / (1.0 + CalbB * Ka)
    CHBf = CHB * 1.0 / (1.0 + CalbB * Ka)
    CBrBf = CBrB * 1.0 / (1.0 + CalbB * Ka)
    CTBf = CTB * 1.0 / (1.0 + CalbB * Ka)
    CSKBf = CSKB * 1.0 / (1.0 + CalbB * Ka)
   
    
    #Calculation of free concentrations in organ interstitial fluid
    CKFf = CKF * 1.0 / (1.0 + CalbKF * Ka)
    CLFf = CLF * 1.0 / (1.0 + CalbLF * Ka)
    
    CSTFf = CSTF * 1.0 / (1.0 + CalbSTF * Ka)
    CINFf = CINF * 1.0 / (1.0 + CalbINF * Ka)
    
    CMFf = CMF * 1.0 / (1.0 + CalbMF * Ka)
    CAFf = CAF * 1.0 / (1.0 + CalbAF * Ka)
    CRFf = CRF * 1.0 / (1.0 + CalbRF * Ka)
    CLuFf = CLuF * 1.0 / (1.0 + CalbLuF * Ka)
    
    CSPFf = CSPF * 1.0 / (1.0 + CalbSPF * Ka)
    CHFf = CHF * 1.0 / (1.0 + CalbHF * Ka)
    CBrFf = CBrF * 1.0 / (1.0 + CalbBrF * Ka)
    CTFf = CTF * 1.0 / (1.0 + CalbTF * Ka)
    
    CSKFf = CSKF * 1.0 / (1.0 + CalbSKF * Ka)
    
    #Calculation of free concentrations in alveolar lining fluid 
    CLuAFf = CLuAF * 1.0 / (1.0 + CalbLuAF * Ka)
    
    #Calculation of free concentrations in organ where we have tissue binding
    CKTf = CKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp)
    CLTf = CLT * 1.0 / (1.0 + CLfabpLT * KLfabp)
    
    QBLtot = QBL + (QBSP-QBSP/500) + (QBIN-QBIN/500) + (QBST-QBST/500)
    
    #Arterial Blood
    dMBart = (QBLu-QBLu/500)*CLuBf - CBfart*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+QBST+QBIN+QBT+QBSK)-QGFR*CBfart
    
    #Venous Blood
    dMBven = kUAB * CUA - CBfven *QBLu +  CLNf*(QBK/500+QBLtot/500+QBM/500+QBA/500+QBR/500+QBLu/500+QBH/500+QBT/500+QBSK/500+
                                       QBST/500+QBSP/500+QBIN/500)+
                                    (QBK-QBK/500)*CKBf + (QBLtot-QBLtot/500)*CLBf + 
                                    (QBM-QBM/500)*CMBf + (QBA-QBA/500)*CABf + (QBR-QBR/500)*CRBf+
                                    (QBH-QBH/500)*CHBf + QBBr*CBrBf+ (QBT-QBT/500)*CTBf + (QBSK-QBSK/500)*CSKBf
    
    
    #Lymph nodes
    dMLN = CKFf*QBK/500 + CLFf*QBLtot/500 + CMFf*QBM/500 + CAFf*QBA/500+
           CRFf*QBR/500 + CLuFf*QBLu/500 + CHFf*QBH/500 + 
           CTFf*QBT/500 + CSKFf*QBSK/500 + CSTFf*QBST/500 + CINFf*QBIN/500 + CSPFf*QBSP/500-
           CLNf*(QBK/500+QBLtot/500+QBM/500+QBA/500+QBR/500+QBLu/500+QBH/500+QBT/500+QBSK/500+
                 QBST/500+QBSP/500+QBIN/500)
   
    
    #Kidney
    
    #blood subcompartment
    dMBK = QBK*CBfart - (QBK-QBK/500)*CKBf - PeffK*AK*(CKBf-CKFf) - CKBf*QBK/500 
    #interstitial fluid subcompartment
    dMKF = CKBf*QBK/500 - CKFf*QBK/500 + PeffK*AK*(CKBf-CKFf) - kKFKT*(CKFf-CKTf) -
            (VmK_Oat1*CKFf/KmK_Oat1+CKFf) - (VmK_Oat3*CKFf/KmK_Oat3+CKFf) #+ (VmK_Osta*CKTf/KmK_Osta+CKTf)
    #Kidney proximal tubule cells subcompartment
    dMKT = kKFKT*(CKFf-CKTf) - kFKT*(CKTf - CFil) + (VmK_Oatp*CFil/KmK_Oatp+CFil) +
            (VmK_Oat1*CKFf/KmK_Oat1+CKFf) + (VmK_Oat3*CKFf/KmK_Oat3+CKFf) # + (VmK_Osta*CKTf/KmK_Osta+CKTf)
    dMFil =  QGFR*CBfart+ kFKT*(CKTf - CFil) - (VmK_Oatp*CFil/KmK_Oatp+CFil) - (Qurine*CFil)
    dMurine = Qurine*CFil
    
    #Liver
    
    #blood subcompartment
    dMBL = QBL*CBfart + (QBSP-QBSP/500)*CSPBf + (QBIN-QBIN/500)*CINBf + (QBST-QBST/500)*CSTBf - PeffL*AL*(CLBf-CLFf) - CLBf*QBLtot/500 - (QBLtot-QBLtot/500)*CLBf
    #interstitial fluid subcompartment 
    dMLF = CLBf*QBLtot/500 - CLFf*QBLtot/500 + PeffL*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
          (VmL_Oatp*CLFf/KmL_Oatp+CLFf) - (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf)
    #Liver tissue subcompartment
    dMLT = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/KmL_Oatp+CLFf) + 
                     (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf) - Qbile*CLTf
    
    # Feces
    dMfeces = Qfeces*CINL
    
    #Stomach
    #blood subcompartment
    dMBST = QBST*CBfart - (QBST-QBST/500)*CSTBf - PeffST*AST*(CSTBf-CSTFf) - CSTBf*QBST/500
    #interstitial fluid subcompartment 
    dMSTF = CSTBf*QBST/500 - CSTFf*QBST/500 + PeffST*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTT)
    #Stomach tissue subcompartment
    dMSTT = kSTFSTT*(CSTFf-CSTT) + kabST*CSTL
    #Stomach lumen
    dMSTL = -kabST*CSTL - QGE*CSTL + CLEal*CLuAFf + CLEua*CUA
    
    #Intestine
    #blood subcompartment
    dMBIN = QBIN*CBfart - (QBIN-QBIN/500)*CINBf - PeffIN*AIN*(CINBf-CINFf) - CINBf*QBIN/500
    #interstitial fluid subcompartment 
    dMINF = CINBf*QBIN/500 - CINFf*QBIN/500 + PeffIN*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINT) 
    #Intestine tissue subcompartment
    dMINT = kINFINT*(CINFf-CINT) + kabIN*CINL
    #Intestine lumen
    dMINL = QGE*CSTL - (Qfeces*CINL) - kabIN*CINL + Qbile*CLTf

    
    #Muscle
    
    #blood subcompartment
    dMBM = QBM*CBfart - (QBM-QBM/500)*CMBf - PeffM*AM*(CMBf-CMFf) - CMBf*QBM/500
    #interstitial fluid subcompartment 
    dMMF = CMBf*QBM/500 - CMFf*QBM/500 +PeffM*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMT)
    #Muscle tissue subcompartment 
    dMMT = kMFMT*(CMFf- CMT)
    
    #Adipose
    
    #blood subcompartment
    dMBA = QBA*CBfart - (QBA-QBA/500)*CABf - PeffA*AA*(CABf-CAFf) - CABf*QBA/500
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
    
    #Upper airways
    
    dMUA = IVR*EC*Dua - CLEua*CUA  - kUAB * CUA 
    
    
    #Lung
    
    #blood subcompartment
    dMBLu = CBfven *QBLu - (QBLu-QBLu/500)*CLuBf - CLuBf*QBLu/500 - PeffLu*ALu*(CLuBf-CLuFf)
    #interstitial fluid subcompartment (it has albumin???)
    dMLuF = CLuBf*QBLu/500 - CLuFf*QBLu/500 + PeffLu*ALu*(CLuBf-CLuFf) + kLuTLuF*(CLuT-CLuFf)
    #Lung tissue subcompartment
    dMLuT = - kLuTLuF*(CLuT-CLuFf) + kLuAFLuT*CLuAFf  
    #Alveolar lining fluid
    dMLuAF = - CLEal*CLuAFf - kLuAFLuT*CLuAFf #+ IVR*EC*Dal 
    
    #Spleen
    
    #blood subcompartment
    dMBSP = QBSP*CBfart - (QBSP-QBSP/500)*CSPBf - PeffSP*ASP*(CSPBf-CSPFf) - CSPBf*QBSP/500
    #interstitial fluid subcompartment 
    dMSPF = CSPBf*QBSP/500 - CSPFf*QBSP/500 + PeffSP*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf -CSPT) 
    #Spleen tissue subcompartment 
    dMSPT = kSPFSPT*(CSPFf -CSPT) 
    
    #Heart
    
    #blood subcompartment
    dMBH = QBH*CBfart - (QBH-QBH/500)*CHBf - PeffH*AH*(CHBf-CHFf) - CHBf*QBH/500 
    #interstitial fluid subcompartment 
    dMHF = CHBf*QBH/500 - CHFf*QBH/500 + PeffH*AH*(CHBf-CHFf) - kHFHT*(CHFf -CHT) 
    #Heart tissue subcompartment 
    dMHT = kHFHT*(CHFf -CHT) 
    
    #Brain
    
    #blood subcompartment
    dMBBr = QBBr*CBfart - QBBr*CBrBf - PeffBr*ABr*(CBrBf-CBrFf)
    #interstitial fluid subcompartment 
    dMBrF = PeffBr*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf -CBrT) 
    #Brain tissue subcompartment 
    dMBrT = kBrFBrT*(CBrFf -CBrT) 
    
    #Testis
    
    #blood subcompartment
    dMBT = QBT*CBfart - (QBT-QBT/500)*CTBf - PeffT*AT*(CTBf-CTFf) - CTBf*QBT/500
    #interstitial fluid subcompartment 
    dMTF = CTBf*QBT/500 - CTFf*QBT/500 + PeffT*AT*(CTBf-CTFf) - kTFTT*(CTFf -CTT) 
    #Testis tissue subcompartment 
    dMTT = kTFTT*(CTFf -CTT) 
    
    #Skin
    
    #blood subcompartment
    dMBSK = QBSK*CBfart - (QBSK-QBSK/500)*CSKBf - PeffSK*ASK*(CSKBf-CSKFf) - CSKBf*QBSK/500
    #interstitial fluid subcompartment
    dMSKF = CSKBf*QBSK/500 - CSKFf*QBSK/500 + PeffSK*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf -CSKT)
    #Skin tissue subcompartment
    dMSKT = kSKFSKT*(CSKFf -CSKT)
    
    #Vurine, Vfeces
    dVurine = Qurine
    dVfeces = Qfeces
    
    #Concentration calculation in each compartment 
    Cven <- CBven
    Cart <- CBart
    Cblood <- (MBven + MBart)/(Vven+Vart)
    Mblood <- MBven + MBart
    Cplasma <- (MBven + MBart)/ Vplasma
    
    Ckidney <- (MBK + MKF+ MKT)/(VKB+VKF+VKT)
    Mkidney <- MBK + MKF+ MKT
    
    Cliver <- (MBL + MLF+ MLT)/(VLB+VLF+VKT)
    Mliver <- MBL + MLF+ MLT
    
    Cstomach <-  (MBST + MSTF+ MSTT)/(VSTB+VSTF+VSTT)
    Cintestine <-  (MBIN + MINF+ MINT+MINL)/(VINB+VINF+VINT+VINL)
    Cmuscle <-  (MBM + MMF+ MMT)/(VMB+VMF+VMT)
    Cadipose <-  (MBA + MAF+ MAT)/(VAB+VAF+VAT)
    
    CUpperair <- MUA/VUA
    CalveolarLF <- MAF/VLuAF
    Clungs <-  (MBLu + MLuF+ MLuAF+ MLuT)/(VLuB+VLuF+VLuAF+VLuT)
    
    Crest <-  (MBR + MRF+ MRT)/(VRB+VRF+VRT)
    Ccarcass <- (MBM + MMF+ MMT+MBA + MAF+ MAT +MBR + MRF+ MRT)/(VM+VA+VR)
    Cfeces <- Mfeces/(Vfeces*feces_density)
    Curine <- Murine/Vurine
    Cspleen <-  (MBSP + MSPF+ MSPT)/(VSPB+VSPF+VSPT)
    Cheart <-  (MBH + MHF+ MHT)/(VHB+VHF+VHT)
    
    Cbrain <-  (MBBr + MBrF+ MBrT)/(VBrB+VBrF+VBrT)
    Mbrain <- MBBr + MBrF+ MBrT
    
    Ctestis <-  (MBT + MTF+ MTT)/(VTB+VTF+VTT)
    Cskin <-  (MBSK + MSKF+ MSKT)/(VSKB+VSKF+VSKT)
   
    list(c( 'dMBart'=dMBart, 'dMBven'=dMBven, 'dMLN'=dMLN, 'dMBK'=dMBK, 
            'dMKF'=dMKF, 'dMKT'=dMKT,
            'dMFil'=dMFil, 'dMurine'=dMurine, 'dMBL'=dMBL, 
            'dMLF'=dMLF, 'dMLT'=dMLT, 
            'dMSTL'=dMSTL,'dMINL'=dMINL,'dMfeces'=dMfeces,
            
            'dMBST'=dMBST, 'dMSTF'=dMSTF, 'dMSTT'=dMSTT,
            'dMBIN'=dMBIN, 'dMINF'=dMINF, 'dMINT'=dMINT,
           
            'dMBM'=dMBM, 'dMMF'=dMMF, 'dMMT'=dMMT,
            'dMBA'=dMBA, 'dMAF'=dMAF, 
            'dMAT'=dMAT, 'dMBR'=dMBR, 'dMRF'=dMRF,'dMRT'=dMRT,
            'dMUA'=dMUA, 'dMBLu'=dMBLu, 'dMLuF'=dMLuF,'dMLuT'=dMLuT,'dMLuAF' = dMLuAF,
           
            'dMBSP'=dMBSP, 'dMSPF'=dMSPF, 'dMSPT'=dMSPT,
            'dMBH'=dMBH, 'dMHF'=dMHF, 'dMHT'=dMHT,
            'dMBBr'=dMBBr, 'dMBrF'=dMBrF, 'dMBrT'=dMBrT,
            'dMBT'=dMBT, 'dMTF'=dMTF, 'dMTT'=dMTT,
            'dMBSK'=dMBSK, 'dMSKF'=dMSKF, 'dMSKT'=dMSKT,
            'dVurine'=dVurine, 'dVfeces'=dVfeces
            ),
                  
           'CKFf'=CKFf, 'CLNf'=CLNf, 'CLFf'=CLFf,
           'CMFf'=CMFf,'CAFf'=CAFf, 'CRFf'=CRFf,
           'CLuFf'=CLuFf, 'CLuAFf'=CLuAFf,
           'CBfart'=CBfart, 'CKBf'=CKBf, 'CLBf'=CLBf, 'CMBf'=CMBf, 'CABf'=CABf,
           'CRBf'=CRBf, 'CFil'=CFil, 
           'CKTf'=CKTf, 'CLTf'=CLTf,
           'CSTL'=CSTL,'CINL'=CINL,
           'CMT'=CMT, 'CAT'=CAT, 'CRT'=CRT, 
                  
            'Cven'= Cven, 'Cart' = Cart,'Cblood' = Cblood, 'Mblood'= Mblood, 'Cplasma'= Cplasma,
            'Ckidney'= Ckidney, 'Mkidney'= Mkidney,
            'Cliver'= Cliver, 'Mliver'= Mliver,
            'Cstomach'= Cstomach, 'Cintestine'= Cintestine, 'CUpperair'= CUpperair, 'CalveolarLF'=CalveolarLF,
            'Cmuscle'= Cmuscle, 'Cadipose'= Cadipose, 
            'Clungs' = Clungs, 'Crest'= Crest,'Ccarcass' = Ccarcass,
            'Cfeces'= Cfeces, 'Curine'= Curine,
            'Cspleen'= Cspleen, 'Cheart'= Cheart,
            'Cbrain'= Cbrain, 'Mbrain'= Mbrain,
            'Ctestis'= Ctestis, 'Cskin'= Cskin
                
         )
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBart <- 0; MBven <-0;  MLN <-0; MBK <-0; MKF <-0; MKT <-0; MFil <-0; Murine <-0; MBL <-0;
    MLF <-0; MLT <-0; 
    MSTL <-0;  MINL <-0;
    Mfeces <-0; MBST <-0; MSTF <-0; MSTT <-0; MBIN <-0; MINF <-0; MINT <-0;
    MBM <-0; MMF <-0; MMT <-0; MBA <-0; MAF <-0
    MAT <-0; MBR <-0; MRF <-0; MRT <-0; MUA <-0; MBLu <- 0; MLuF <- 0; MLuT <- 0;MLuAF <- 0;
    MBSP <-0; MSPF <-0; MSPT <-0; MBH <-0; MHF <-0; MHT <-0; 
    MBBr <-0; MBrF <-0; MBrT <-0; MBT <-0; MTF <-0; MTT <-0;
    MBSK <-0; MSKF <-0; MSKT <-0; Vurine <-0; Vfeces <-0
    
    return(c('MBart'=MBart, 'MBven'=MBven, 'MLN'=MLN, 'MBK'=MBK, 'MKF'=MKF, 'MKT'=MKT,
             'MFil'=MFil, 'Murine'=Murine, 'MBL'=MBL, 'MLF'=MLF, 'MLT'=MLT, 
             'MSTL'=MSTL, 'MINL'=MINL, 'Mfeces'=Mfeces, 
             
             'MBST'=MBST, 'MSTF'=MSTF, 'MSTT'=MSTT,
             'MBIN'=MBIN, 'MINF'=MINF, 'MINT'=MINT,
             
             'MBM'=MBM, 'MMF'=MMF, 'MMT'=MMT,
             'MBA'=MBA, 'MAF'=MAF, 'MAT'=MAT,
             'MBR'=MBR, 'MRF'=MRF, 'MRT'=MRT, 'MUA'=MUA,
             'MBLu'=MBLu, 'MLuF'=MLuF,'MLuT' = MLuT,'MLuAF'=MLuAF, 
             
             'MBSP'=MBSP, 'MSPF'=MSPF, 'MSPT'=MSPT,
             'MBH'=MBH, 'MHF'=MHF, 'MHT'=MHT,
             'MBBr'=MBBr, 'MBrF'=MBrF, 'MBrT'=MBrT,
             'MBT'=MBT, 'MTF'=MTF, 'MTT'=MTT,
             'MBSK'=MBSK, 'MSKF'=MSKF, 'MSKT'=MSKT,
             'Vurine'=Vurine, 'Vfeces'=Vfeces
             
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
        events <- list(data = rbind(data.frame(var = c("MSTL"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
        
      }else if (admin.type == "inh"){
        events <- list(data = rbind(data.frame(var = c("MLuAF"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "nasal"){
        events <- list(data = rbind(data.frame(var = c("MUA"),  time = admin.time, 
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
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)
  
  
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male blood
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Gustafsson (2022) Inhalation male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "inh"
   
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Gustafsson (2022) Inhalation male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_gus_INH_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_INH_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_gus_INH_Mblood <- unlist(preds_gus_INH_Mblood) /1000 #convert ug/kg to ug/g
  
  
  obs_gus_INH_Mblood <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  score[1] <- AAFE(predictions = preds_gus_INH_Mblood, observations = obs_gus_INH_Mblood)
  
 
  ##########################
  #-------------------------
  # Gustafsson Oral male blood
  #-------------------------
  ##########################
  # Set up simulations for the 2nd case, i.e. Gustafsson (2022) Inhalation male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 0.364 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Gustafsson (2022) Oral male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_gus_OR_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_OR_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_gus_OR_Mblood <- unlist(preds_gus_OR_Mblood) /1000 #convert ug/kg to ug/g
  
  
  obs_gus_OR_Mblood <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  score[2] <- AAFE(predictions = preds_gus_OR_Mblood, observations = obs_gus_OR_Mblood)
  
  
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 3rd case, i.e. Gustafsson (2022) Inhalation male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "inh"
   
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of Gustafsson (2022) Inhalation male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CalveolarLF","Cliver","Clung", "Ckidney")
  
  preds_gus_INH_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_INH_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_gus_INH_Mtissues <- unlist(preds_gus_INH_Mtissues) /1000 #convert ug/kg to ug/g
  
  
  obs_gus_INH_Mtissues <- c( exp_data[exp_data$Tissue == "ALF", "concentration"],
                             exp_data[exp_data$Tissue == "Liver", "concentration"],
                             exp_data[exp_data$Tissue == "Lung", "concentration"], 
                             exp_data[exp_data$Tissue == "Kidney", "concentration"]) 
  
  
  
  score[3] <- AAFE(predictions = preds_gus_INH_Mtissues, observations = obs_gus_INH_Mtissues)
  
  
  ##########################
  #-------------------------
  # Gustafsson Oral male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 4th case, i.e. Gustafsson (2022) Inhalation male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 0.364 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
   
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of Gustafsson (2022) oral male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CalveolarLF","Cliver","Clung", "Ckidney")
  
  preds_gus_OR_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_OR_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_gus_OR_Mtissues <- unlist(preds_gus_OR_Mtissues) /1000 #convert ug/kg to ug/g
  
  
  obs_gus_OR_Mtissues <- c( exp_data[exp_data$Tissue == "ALF", "concentration"],
                             exp_data[exp_data$Tissue == "Liver", "concentration"],
                             exp_data[exp_data$Tissue == "Lung", "concentration"], 
                             exp_data[exp_data$Tissue == "Kidney", "concentration"]) 
  
  
  
  score[4] <- AAFE(predictions = preds_gus_OR_Mtissues, observations = obs_gus_OR_Mtissues)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single low
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 1.2 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
                
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,9,0.04)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of Hinderliter Inhalation male single low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CPlasma")
  
  preds_hind_INH_Mblood_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_low [[i]] <- solution[round(solution$time) %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Mblood_low <- unlist(preds_hind_INH_Mblood_low) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Mblood_low <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[5] <- AAFE(predictions = preds_hind_INH_Mblood_low, observations = obs_hind_INH_Mblood_low)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single medium
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 9.8 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of Hinderliter Inhalation male single medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Plasma")
  
  preds_hind_INH_Mblood_medium <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Mblood_medium <- unlist(preds_hind_INH_Mblood_medium) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Mblood_medium <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[6] <- AAFE(predictions = preds_hind_INH_Mblood_medium, observations = obs_hind_INH_Mblood_medium)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single high
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single high
  BW <- 0.25  #kg, not reported in the study
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 27 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df7=========================================================
  
  exp_data <- dataset$df7 # retrieve data of Hinderliter Inhalation male single high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Plasma")
  
  preds_hind_INH_Mblood_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Mblood_high <- unlist(preds_hind_INH_Mblood_high) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Mblood_high <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[7] <- AAFE(predictions = preds_hind_INH_Mblood_high, observations = obs_hind_INH_Mblood_high)
  

  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single low
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "F"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 1.2 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,9,0.04)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df8=========================================================
  
  exp_data <- dataset$df8 # retrieve data of Hinderliter Inhalation male single low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Plasma")
  
  preds_hind_INH_Fblood_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_low [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Fblood_low <- unlist(preds_hind_INH_Fblood_low) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Fblood_low <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[8] <- AAFE(predictions = preds_hind_INH_Fblood_low, observations = obs_hind_INH_Fblood_low)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single medium
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "F"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 9.8 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df9=========================================================
  
  exp_data <- dataset$df9 # retrieve data of Hinderliter Inhalation male single medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Plasma")
  
  preds_hind_INH_Fblood_medium <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Fblood_medium <- unlist(preds_hind_INH_Fblood_medium) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Fblood_medium <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[9] <- AAFE(predictions = preds_hind_INH_Fblood_medium, observations = obs_hind_INH_Fblood_medium)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single high
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single high
  BW <- 0.25  #kg, not reported in the study
  sex <- "F"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 27 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df10=========================================================
  
  exp_data <- dataset$df10 # retrieve data of Hinderliter Inhalation male single high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Plasma")
  
  preds_hind_INH_Fblood_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  preds_hind_INH_Fblood_high <- unlist(preds_hind_INH_Fblood_high) /1000 #convert ug/kg to ug/g
  
  
  obs_hind_INH_Fblood_high <- c(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[10] <- AAFE(predictions = preds_hind_INH_Fblood_high, observations = obs_hind_INH_Fblood_high)
  
  
  
  ########################################################################################
  # Estimate final score
  
  final_score <- mean(score, na.rm = TRUE)
  return(final_score)
  
}

################################################################################

setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")
MW <- 414.07 #g/mol

# Read data

gus_INH_Mblood <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Plasma Male rats_Inhalation.xlsx")
gus_OR_Mblood <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Plasma Male rats_Oral.xlsx")
gus_INH_Mtissues <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Tissues Male rats_Inhalation.xlsx")
gus_OR_Mtissues <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Tissues Male rats_Oral.xlsx")
hind_INH_Mblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Low_dose.xlsx")
hind_INH_Mblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Medium_dose.xlsx")
hind_INH_Mblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_High_dose.xlsx")
hind_INH_Fblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Low_dose.xlsx")
hind_INH_Fblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Medium_dose.xlsx")
hind_INH_Fblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_High_dose.xlsx")


dataset <- list("df1" = gus_INH_Mblood, "df2" = gus_OR_Mblood, "df3" = gus_INH_Mtissues, "df4" = gus_OR_Mtissues,
                "df5" = hind_INH_Mblood_low, "df6" = hind_INH_Mblood_medium, "df7" = hind_INH_Mblood_high,
                "df8" = hind_INH_Fblood_low, "df9" = hind_INH_Fblood_medium, "df10" = hind_INH_Fblood_high)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX",#"NLOPT_LN_NEWUOA","NLOPT_LN_SBPLX"
              "xtol_rel" = 1e-07,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 1, 
              "print_level" = 1)

# Create initial conditions (zero initialisation)
#Parameter names:
# Male RAFOatp_k, Male RAFOat1, Male RAFOat3, Male RAFOatp_l,Male RAFNtcp
# Female RAFOatp_k, Female RAFOat1, Female RAFOat3, Female RAFOatp_l,female RAFNtcp
#  bile_correction_factor, 11 correction factors for permeabilities

N_pars <- 4 # Number of parameters to be fitted
#fit <- log(rep(1,N_pars))
fit <- rep(1,N_pars)



# Run the optimization algorithmm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             lb	= rep(log(1e-20), N_pars),
                             ub = rep(log(1e10), N_pars),
                             opts = opts,
                             dataset = dataset)

#estimated_params <- exp(optimizer$solution)
estimated_params <- exp(optimizer$solution)


 
 
 #################################################################################
 #--------------------------------------------------------------------------------
 #                                Gustafsson Inhalation male blood
 #-------------------------------------------------------------------------------
 #################################################################################
 
 # Set up simulations for the 1st case, i.e. Gustafsson Inhalation male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
inhalation_params=estimate_BF_TV(sex)
BF = inhalation_params["BF"];TV = inhalation_params["TV"]
admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "inh"
 
 
user_input <- list('BW'=BW,
                    "admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "estimated_params" = estimated_params,
                    "sex" = sex)
 
 
 params <- create.params(user_input)
 inits <- create.inits(params)
 events <- create.events(params)
 
 
 sample_time=seq(0,48,0.2)
 solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                     y = inits, parms = params,events = events,
                                     method="lsodes",rtol = 1e-03, atol = 1e-03))
 
 preds_gus_INH_Mblood <-  solution[, c("time", "Cplasma")]
 
 
 #################################################################################
 #--------------------------------------------------------------------------------
 #                                Gustafsson Oral male blood
 #-------------------------------------------------------------------------------
 #################################################################################
 
 # Set up simulations for the 2nd case, i.e. Gustafsson Oral male blood
 BW <- 0.5125  #kg, from Gustafsson et al., 2022
 sex <- "M"
 inhalation_params=estimate_BF_TV(sex)
 BF = inhalation_params["BF"];TV = inhalation_params["TV"]
 admin.dose_per_g <- 0.364 # administered dose in mg PFOA/kg BW 
 admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
 admin.time <- 0 #time when doses are administered, in hours
 admin.type <- "oral"
  
 
 user_input <- list('BW'=BW,
                    "admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "estimated_params" = estimated_params,
                    "sex" = sex)
 
 
 params <- create.params(user_input)
 inits <- create.inits(params)
 events <- create.events(params)
 
 
 sample_time=seq(0,48,0.2)
 solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                     y = inits, parms = params,events = events,
                                     method="lsodes",rtol = 1e-03, atol = 1e-03))
 
 preds_gus_OR_Mblood <-  solution[, c("time", "Cplasma")]
 
 
 #################################################################################
 #--------------------------------------------------------------------------------
 #                                Gustafsson Inhalation male tissues
 #-------------------------------------------------------------------------------
 #################################################################################
 
 # Set up simulations for the 3rd case, i.e. Gustafsson Inhalation male tissues
 BW <- 0.5125  #kg, from Gustafsson et al., 2022
 sex <- "M"
 inhalation_params=estimate_BF_TV(sex)
 BF = inhalation_params["BF"];TV = inhalation_params["TV"]
 admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
 admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
 admin.time <- 0 #time when doses are administered, in hours
 admin.type <- "inh"
  
 
 user_input <- list('BW'=BW,
                    "admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "estimated_params" = estimated_params,
                    "sex" = sex)
 
 
 params <- create.params(user_input)
 inits <- create.inits(params)
 events <- create.events(params)
 
 
 sample_time=seq(0,48,0.2)
 solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                     y = inits, parms = params,events = events,
                                     method="lsodes",rtol = 1e-03, atol = 1e-03))
 
  preds_gus_INH_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungs", "Ckidney")]
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                                Gustafsson Oral male tissues
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 4th case, i.e. Gustafsson Inhalation male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 0.364 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time=seq(0,48,0.2)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_gus_OR_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungs", "Ckidney")]
  
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation male single low
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 1.2 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
   
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,9,0.04)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_hind_INH_Mblood_low <-  solution[, c("time", "Cplasma")]
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation male single medium
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 6th case, i.e. Hinderliter Inhalation male single medium
  BW <- 0.25  #kg, not reported in the study
  sex <- "M" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 9.8 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,30,0.1)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_hind_INH_Mblood_medium <-  solution[, c("time", "Cplasma")]
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation male single high
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 7th case, i.e. Hinderliter Inhalation male single high
  BW <- 0.25  #kg, not reported in the study
  sex <- "M"
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 27 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
   
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,30,0.1)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_hind_INH_Mblood_high <-  solution[, c("time", "Cplasma")]
  
 
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation female single low
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 8th case, i.e. Hinderliter Inhalation male single low
  BW <- 0.25  #kg, not reported in the study
  sex <- "F" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 1.2 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,9,0.04)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_hind_INH_Fblood_low <-  solution[, c("time", "Cplasma")]
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation female single medium
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 9th case, i.e. Hinderliter Inhalation female single medium
  BW <- 0.25  #kg, not reported in the study
  sex <- "F" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 9.8 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,30,0.1)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  preds_hind_INH_Fblood_medium <-  solution[, c("time", "Cplasma")]
  
  
  #################################################################################
  #--------------------------------------------------------------------------------
  #                        Hinderliter Inhalation female single high
  #-------------------------------------------------------------------------------
  #################################################################################
  
  # Set up simulations for the 10th case, i.e. Hinderliter Inhalation female single high
  BW <- 0.25  #kg, not reported in the study
  sex <- "F" 
  inhalation_params=estimate_BF_TV(sex)
  BF = inhalation_params["BF"];TV = inhalation_params["TV"]
  admin.dose_per_g <- 27 # administered dose in mg/m^3
  admin.dose <- admin.dose_per_g*6*BF*TV #ug PFOA, for 6h inhalation
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  sample_time= seq(0,30,0.1)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,events = events,
                                      method="lsodes",rtol = 1e-3, atol = 1e-3))
  
  preds_hind_INH_Fblood_high <-  solution[, c("time", "Cplasma")]
  
  
 #convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
 
 
 preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] <- preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] /1000
 preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] <- preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] /1000
 preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] <- preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] /1000
 preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] <- preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] /1000
 preds_hind_INH_Mblood_low[,2:dim(preds_hind_INH_Mblood_low)[2]] <- preds_hind_INH_Mblood_low[,2:dim(preds_hind_INH_Mblood_low)[2]] /1000
 preds_hind_INH_Mblood_medium[,2:dim(preds_hind_INH_Mblood_medium)[2]] <- preds_hind_INH_Mblood_medium[,2:dim(preds_hind_INH_Mblood_medium)[2]] /1000
 preds_hind_INH_Mblood_high[,2:dim(preds_hind_INH_Mblood_high)[2]] <- preds_hind_INH_Mblood_high[,2:dim(preds_hind_INH_Mblood_high)[2]] /1000
 preds_hind_INH_Fblood_low[,2:dim(preds_hind_INH_Fblood_low)[2]] <- preds_hind_INH_Fblood_low[,2:dim(preds_hind_INH_Fblood_low)[2]] /1000
 preds_hind_INH_Fblood_medium[,2:dim(preds_hind_INH_Fblood_medium)[2]] <- preds_hind_INH_Fblood_medium[,2:dim(preds_hind_INH_Fblood_medium)[2]] /1000
 preds_hind_INH_Fblood_high[,2:dim(preds_hind_INH_Fblood_high)[2]] <- preds_hind_INH_Fblood_high[,2:dim(preds_hind_INH_Fblood_high)[2]] /1000
 
 
 # ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 
 
 # Function that creates a plot given a compartment name and the respective predictions and observations
 create.plots <- function(predictions, observations, compartment){  
   #Colours of observations and predictions
   cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
   
   ggplot(data = predictions)+
     geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                           color = '"predictions"'),  size=1.5,alpha = 0.7) +
     geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                              color='"Observations"'), size=4)+
     labs(title = rlang::expr(!!compartment), 
          y = expression("PFOA concentration (" * mu* "g/g tissue)" ),
          x = "Time (hours)")+
     theme(plot.title = element_text(hjust = 0.5))+
     scale_color_manual("", values=cls,
                        guide = guide_legend(override.aes =
                                               list(shape = c(16,NA),
                                                    linetype = c(0,1))))+
     theme_light() + 
     theme(legend.position=c(1,1), 
           legend.justification=c(0, 1), 
           legend.key.size = unit(1.5, 'cm'),  
           legend.title = element_text(size=14),
           axis.title=element_text(size=14),
           legend.text = element_text(size=14)
     )
   
 }
 
 

 # Convert Gustafsson Inhalation male blood from long to wide format using reshape
 experiment_inh_1 <- reshape(gus_INH_Mblood[c("Tissue" ,"Time_hours", 
                                          "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_1) <- c("Time",unique(gus_OR_Mblood$Tissue))
 
 # Convert Gustafsson Oral male blood from long to wide format using reshape
 experiment_inh_2 <- reshape(gus_OR_Mblood[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_2) <- c("Time",unique(gus_OR_Mblood$Tissue))
 
 
 # Convert Gustafsson Inhalation male tissues from long to wide format using reshape
 experiment_inh_3 <- reshape(gus_INH_Mtissues[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_3) <- c("Time",unique(gus_INH_Mtissues$Tissue))
 
 
 # Convert Gustafsson Oral male tissues from long to wide format using reshape
 experiment_inh_4 <- reshape(gus_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_4) <- c("Time",unique(gus_OR_Mtissues$Tissue))
 
 
 # Convert Hinderliter Inhalation male single low from long to wide format using reshape
 experiment_inh_5 <- reshape(hind_INH_Mblood_low[c("Tissue" ,"Time_hours", 
                                               "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_5) <- c("Time",unique(hind_INH_Mblood_low$Tissue))
 
 
 # Convert Hinderliter Inhalation male single medium from long to wide format using reshape
 experiment_inh_6 <- reshape(hind_INH_Mblood_medium[c("Tissue" ,"Time_hours", 
                                                   "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_6) <- c("Time",unique(hind_INH_Mblood_medium$Tissue))
 
 
 # Convert Hinderliter Inhalation male single high from long to wide format using reshape
 experiment_inh_7 <- reshape(hind_INH_Mblood_high[c("Tissue" ,"Time_hours", 
                                                      "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_7) <- c("Time",unique(hind_INH_Mblood_high$Tissue))
 
 
 # Convert Hinderliter Inhalation female single low from long to wide format using reshape
 experiment_inh_8 <- reshape(hind_INH_Fblood_low[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_8) <- c("Time",unique(hind_INH_Fblood_low$Tissue))
 
 
 # Convert Hinderliter Inhalation female single medium from long to wide format using reshape
 experiment_inh_9 <- reshape(hind_INH_Fblood_medium[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_9) <- c("Time",unique(hind_INH_Fblood_medium$Tissue))
 
 
 # Convert Hinderliter Inhalation female single high from long to wide format using reshape
 experiment_inh_10 <- reshape(hind_INH_Fblood_high[c("Tissue" ,"Time_hours", 
                                                      "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
 colnames(experiment_inh_10) <- c("Time",unique(hind_INH_Fblood_high$Tissue))
 
 
 # Put the experiments in a list
 experiments_inh <- list(experiment_inh_1 = experiment_inh_1, experiment_inh_2 = experiment_inh_2,
                         experiment_inh_3 = experiment_inh_3, experiment_inh_4 = experiment_inh_4,
                         experiment_inh_5=experiment_inh_5, experiment_inh_6=experiment_inh_6,
                         experiment_inh_7=experiment_inh_7, experiment_inh_8=experiment_inh_8,
                         experiment_inh_9=experiment_inh_9, experiment_inh_10=experiment_inh_10)
 

  # Rename predictions so that they share the same name as the names of the experimental dataframe
 
 
 colnames(preds_gus_INH_Mblood) <- c ("Time", "Plasma")
 colnames(preds_gus_OR_Mblood) <- c ("Time", "Plasma")
 colnames(preds_gus_INH_Mtissues) <- c ("Time", "ALF", "Liver", "Lung", "Kidney")
 colnames(preds_gus_OR_Mtissues) <- c ("Time", "ALF", "Liver", "Lung", "Kidney")
 colnames(preds_hind_INH_Mblood_low) <- c ("Time", "Plasma")
 colnames(preds_hind_INH_Mblood_medium) <- c ("Time", "Plasma")
 colnames(preds_hind_INH_Mblood_high) <- c ("Time", "Plasma")
 colnames(preds_hind_INH_Fblood_low) <- c ("Time", "Plasma")
 colnames(preds_hind_INH_Fblood_medium) <- c ("Time", "Plasma")
 colnames(preds_hind_INH_Fblood_high) <- c ("Time", "Plasma")
 
 
 # Create a list containing the corresponding predictions
 simulations <- list(predictions1 = preds_gus_INH_Mblood, predictions2 = preds_gus_OR_Mblood,
                     predictions3 = preds_gus_INH_Mtissues, predictions4 = preds_gus_OR_Mtissues,
                     predictions5 = preds_hind_INH_Mblood_low, predictions6 = preds_hind_INH_Mblood_medium,
                     predictions7 = preds_hind_INH_Mblood_high, predictions8 = preds_hind_INH_Fblood_low,
                     predictions9 = preds_hind_INH_Fblood_medium, predictions10 = preds_hind_INH_Fblood_high)
                    
 
 # Iterate over all existing experiments and create the accompanying plots
 for(i in 1:length(experiments_inh)){
   # Retrieve the corresponding observations and simulations
   observations <- experiments_inh[[i]]
   predictions <- simulations[[i]]
   # Extract the compartment names
   compartments <- names(predictions)[2:length(predictions)]
   
   # Use lapply to iterate over the column names and create plots
   plots <- lapply(compartments, function(compartment) {
     create.plots(predictions, observations, compartment )
   })
   if(length(compartments) == 1){
     final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                                common.legend = TRUE, legend = "right"))
     
   }else{
     final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                                common.legend = TRUE, legend = "right"))
   }
  
   
   plot.margin=unit(c(0,0,0,0), "pt")
   
   
   # Save the plot with dynamically adjusted dimensions
   ggsave(paste0("experiment_inh", i,".png"), plot = final_plot,
          device = 'png', dpi = 300,
          width = 13,
          height = 10,
          units = "in")
 }

 