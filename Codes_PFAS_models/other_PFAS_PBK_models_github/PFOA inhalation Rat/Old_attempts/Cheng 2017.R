library(deSolve)
getwd()


create.params <- function(user.input){
  with(as.list(user.input),{
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3), assuming the density of tissue is 1 g/mL.
    
    #======Table S1=======#    
    
    PVB <- 54e-6 #13.5 mL/0.244 kg=0.055 mL/kg~55e-6 L/kg
    VB <- PVB * BW #blood volume (kg=L)
    PVplasma <- 31.2e-6 #Brown et al. 1997-Table 5
    Vplasma <- PVplasma * BW #plasma volume (kg=L)
    PVK <- 7.3e-6 #Brown et al. 1997-Table 5
    VK <- PVK * BW #kidney volume (kg=L)
    PVKB <- 0.16 #volume of residual blood in the organ relative to the volume of the organ (Brown et al. 1997)
    VKB <- PVKB * PVK * BW #kidney blood volume (kg=L)
    PVKF <- 0.13 #0.0683 mL/0.539 g = 0.1267 L/kg (Larson et al. 1984~interstitial volume/kidney weight)
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume (kg=L)
    VKT <- VK - VKF #kidney tissue volume (kg=L)
    VFil <- 2.5e-7 #renal filtrate volume (L)
    PVL <- 3.66e-5 #Brown et al. 1997-Table 5
    VL <- PVL * BW #liver volume (kg=L)
    PVLB <- 0.21 #volume of residual blood in the organ relative to the volume of the organ (Brown et al. 1997)
    VLB <- PVLB * PVL * BW #liver blood volume
    PVLF <- 0.049 #4.9% of the total parenchymal volume (Blouin et al. 1977)
    VLF <- PVLF * PVL* BW #liver interstitial fluid volume (kg=L)
    VLT <- VL - VLF #liver tissue volume (kg=L)
    PVbile <- 0.004 #0.4% liver tissue volume (Blouin et al. 1977)
    Vbile <- PVbile * PVL * BW #bile volume
    PVG <- 2.69e-5 #Brown et al. 1997-Table 5
    VG <- PVG * BW #gut volume (kg=L)
    PVGB <- 0.034 #33.9 uL/g tissue (Everett et al. 1956)
    VGB <- PVGB * PVG * BW #gut blood volume (kg=L)
    PVGF <- 0.28 #Barratt & walser 1969
    VGF <- PVGF * PVG * BW #gut interstitial fluid volume (kg=L)
    VGT <- VG - VGF #gut tissue volume (kg=L)
    PVGL <- 4.5e-5 #4.5% BW
    VGL <- PVGL * BW #gut lumen volume (kg=L)
    PVM <- 40.43e-5 #Brown et al. 1997-Table 5
    VM <- PVM * BW #muscle volume (kg=L)
    PVMB <- 0.04 #Brown et al. 1997-Table 30
    VMB <- PVMB * PVM * BW #muscle blood volume (kg=L)
    PVMF <- 0.054
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume
    VMT <- VM - VMF #muscle tissue volume
    PVA <- 7e-5 #Brown et al. 1997-Table 12-Depot weight for BW 0.244 kg
    VA <- PVA * BW #adipose volume (kg=L)
    PVAB <- 0.02 #Brown et al. 1997-Table 30
    VAB <- PVAB * PVA * BW #% adipose blood volume (kg=L)
    PVAF <- 0.174 #Bushnell et al.1998
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume (kg=L)
    VAT <- VA - VAF #adipose tissue volume (kg=L)
    PVR <- 1 / 1e3 - PVB - PVK - PVL - PVG - PVM - PVA
    VR <- PVR * BW #volume of the rest of the body (kg=L)
    PVRB <- 0.036 #Calculated on the weighted average of blood volume of the “rest of body"-Brown et al. 1997
    VRB <- PVRB * PVR * BW #volume of the blood of the rest of body (kg=L)
    PVRF <- 0.18 #Based on data availability, it was assumed to be the weighted average of brain, heart and spleen fluids-Biewald & Billmeier 1978 
    VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body (kg=L)
    VRT <- VR - VRF #tissue volume of the rest of body (kg=L)
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight
    #or weight of corresponding tissue (PAi, m^2/g) and surface area (m^2)
    
    PAK <- 350e-4 
    AK <- PAK * VK * 10^6 #kidney surface area (m^2)
    PAKG <- 68.90e-4 
    AKG <- PAKG * VK * 10^6 #the surface area of glomerular capillary (m^2)
    PAL <- 250e-4 
    AL <- PAL * VL * 10^6 #liver surface area (m^2)
    PAG <- 100e-4 
    AG <- PAG * VG * 10^6 #gut surface area (m^2)
    PAGL <- 4.14 #"The surface area of gut lumen would be 4.14 m2/kg" !!!!???
    AGL <- PAGL * BW #gut lumen surface area (m^2)
    PAM <- 70e-4 
    AM <- PAM * VM * 10^6 #muscle surface area (m^2)
    PAA <- 70e-4
    AA <- PAA * VA * 10^6 #adipose surface area (m^2)
    PAR <- 100e-4
    AR <- PAR * VR * 10^6 #surface area of rest of body (m^2)
    
    #Effective permeability (Peff, in m/s) for blood (B), liver(L), kidney(K),
    #gut(G),adipose(A), muscle(M), rest of body(R)
    
    PeffB <- 4.98e-8
    PeffK <- 4.38e-8
    PeffL <- 5.15e-8
    PeffG <- 2.65e-8
    PeffA <- 2.65e-8
    PeffM <- 2.65e-8
    PeffR <- 2.65e-8
    
    #Steady-state cell-water concentration ratios(CRss, unitless) for gut, liver, and kidney
    
    CRssG <- 3.75
    CRssL <- 7.28
    CRssK <- 6.19
    
    #Blood flow rates (QBi, in m^3/s) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 0.235/60*1e-3 * BW^0.75 #Qc = 0.235×BW^0.75 L/min, where the unit of BW is kg-> m^3/s
    PQBK <- 14.1/100 #percentage 
    QBK <- PQBK * Qcardiac #percentage of cardiac output
    PQBG <- 15.1/100 #percentage 
    QBG <- PQBG * Qcardiac #percentage of cardiac output
    PQBL <- 2.4/100 #percentage 
    QBL <- (PQBL+PQBG) * Qcardiac #percentage of cardiac output
    PQBM <- 27.8/100 #percentage 
    QBM <- PQBM * Qcardiac #percentage of cardiac output
    PQBA <- 7/100 #percentage 
    QBA <- PQBA * Qcardiac #percentage of cardiac output
    PQBR = 1 - PQBK - PQBG - PQBL - PQBM - PQBA
    QBR = PQBR * Qcardiac #percentage of cardiac output
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration
    #rate (GFR), in m^3/s
    
    Qfeces <- 5.63*1e-6/(24*3600) #mL per day -> m^3/s
    PQbile <- 90 #90 mL/d/kg BW
    Qbile <- PQbile * BW * 1e-6/(24*3600) #mL per day/kg -> m^3/s
    PQurine <- 200
    Qurine <- PQurine * BW * 1e-6/(24*3600) #mL per day/kg -> m^3/s
    PQGFR <- 10.74
    QGFR <- PQGFR * BW * 1e-6/60 #mL per min/kg -> m^3/s
    
    #======Table S2=======#  
    #Albumin concentration in blood and interstitial fluid compartments(mol/m^3)
    
    CalbB <- 281e-3*7.8 #486 μmol/L in 13.7 mL blood, C1*v1=C2*V2, v2=7.8 mL plasma, n=7.8 binding sites (Table S3)
    CalbKF <- 243e-3*7.8 #n=7.8 binding sites (Table S3)
    CalbLF <- 243e-3*7.8 #n=7.8 binding sites (Table S3)
    CalbGF <- 146e-3*7.8 #n=7.8 binding sites (Table S3)
    CalbMF <- 146e-3*7.8 #n=7.8 binding sites (Table S3)
    CalbAF <- 73e-3*7.8 #n=7.8 binding sites (Table S3)
    CalbRF <- 73e-3*7.8 #n=7.8 binding sites (Table S3)
    
    #Alpha2mu-globulin concentration in kidney tissue (mol/m^3)
    
    Ca2uKT <- 110e-3
    
    #LFABP concentration in kidney and liver tissue (mol/m^3)
    
    CL_fabpKT <- (2.65e-3)*3 #n=3 binding sites (Table S3)
    CL_fabpKT1 <- CL_fabpKT/3 
    CL_fabpKT2 <- CL_fabpKT/3
    CL_fabpKT3 <- CL_fabpKT/3
    CL_fabpLT <- (133e-3)*3 #n=3 binding sites (Table S3)
    CL_fabpLT1 <- CL_fabpLT/3
    CL_fabpLT2 <- CL_fabpLT/3
    CL_fabpLT3 <- CL_fabpLT/3
    
    #======Table S2=======#    
    #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
    #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details
    
    Ka <- 3.1
    KL_fabp1 <- 120
    KL_fabp2 <- 40.0
    KL_fabp3 <- 19.0
    Ka2u <- 0.5
    
    #Individual rate constants for association and dissociation(s^-1 and m^3/mol*s)
    #Note kon/koff=Keq
    
    koff <- 0.01 #assume koff is 0.01/s
    kon <- koff * Ka #Ka
    kL_fabpon1 <- koff * KL_fabp1
    kL_fabpon2 <- koff * KL_fabp2
    kL_fabpon3 <- koff * KL_fabp3
    kK_fabpon <- koff * Ka2u
    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    kBKF <- ((1/QBK) + 1/(PeffB * AK))^(-1)
    kBF <- PeffB * AKG
    kKFKT <- PeffK * AK
    n <- 5 #enlargement factor of apical membrane of proximal tubule
    
    kFKT <- PeffK * AK * n
    kBLF <- ((1/QBL) + 1/(PeffB * AL))^(-1)
    kLFLT <- PeffL * AL
    kBGF <- ((1/QBG) + 1/(PeffB * AG))^(-1)
    kGFGT <- PeffG * AG
    kGLGT <- PeffG * AGL
    kBMF <- ((1/QBM) + 1/(PeffB * AM))^(-1)
    kMFMT <- PeffM * AM
    kBAF <- ((1/QBA) + 1/(PeffB * AA))^(-1)
    kAFAT <- PeffA * AA
    kBRF <- ((1/QBR) + 1/(PeffB *AR))^(-1)
    kRFRT <- PeffR * AR
    kbileLT <- PeffL * AL
    
    #First-order rate constants (s^-1)
    
    bBKF <- kBKF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bKFB <- kBKF/VKF
    bKFKT <- kKFKT/VKF
    bKTKF <- kKFKT/VKT
    bFKT <- kFKT/VFil
    bKTF <- kFKT/(VKT * CRssK)
    bBF <- QGFR/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bFB <- kBF/VFil
    bBLF <- kBLF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bLFB <- kBLF/VLF
    bLFLT <- kLFLT/VLF
    bLTLF <- kLFLT/VLT
    bbileLT <- kbileLT/Vbile
    bLTbile <- kbileLT/(VLT * CRssL)
    bBGF <- kBGF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bGFB <- kBGF/VGF
    bGFGT <- kGFGT/VGF
    bGTGF <- kGFGT/VGT
    bGLGT <- kGLGT/VGL
    bGTGL <- kGLGT/(VGT * CRssG)
    bBMF <- kBMF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bMFB <- kBMF/VMF
    bMFMT <- kMFMT/VMF
    bMTMF <- kMFMT/VMT
    bBAF <- kBAF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bAFB <- kBAF/VAF
    bAFAT <- kAFAT/VAF
    bATAF <- kAFAT/VAT
    bBRF <- kBRF/(VB+VLB+VKB+VGB+VMB+VAB+VRB)
    bRFB <- kBRF/VRF
    bRFRT <- kRFRT/VRF
    bRTRF <- kRFRT/VRT
    
    #First-order rate constants (s^-1) for protein-mediated transport, see section S3-3 for details
    
    Pbclear <- 2.76e-7
    bclear <- Pbclear * AK/VKF 
    Pbreab <- 1.18e-7
    breab <- n * Pbreab * AK/VFil 
    Pbabs <- 1.78e-7
    babs <- Pbabs * AL/VLF 
    Pbefflux <- 1.38e-7
    befflux <- Pbefflux * AK/VKT 
    
    
    #Conversion between mass and concentration for protein content of tissues
    
    MalbB <- CalbB * (VB+VLB+VKB+VGB+VMB+VAB+VRB)
    MalbKF <- CalbKF * VKF
    ML_fabpKT1 <- CL_fabpKT1 * VKT
    ML_fabpKT2 <- CL_fabpKT2 * VKT
    ML_fabpKT3 <- CL_fabpKT3 * VKT
    MK_fabpKT <- Ca2uKT * VKT
    MalbLF <- CalbLF * VLF
    ML_fabpLT1 <- CL_fabpLT1 * VLT
    ML_fabpLT2 <- CL_fabpLT2 * VLT
    ML_fabpLT3 <- CL_fabpLT3 * VLT
    MalbGF <- CalbGF * VGF
    MalbMF <- CalbMF * VMF
    MalbAF <- CalbAF * VAF
    MalbRF <- CalbRF * VRF
    
    
    
    
    return(list('admin.dose'=admin.dose,'PVB'=PVB, 'VB'=VB, 'PVplsma'=PVplasma, 'Vplasma'=Vplasma, 'PVK'=PVK, 'VK'=VK, 'PVKB'=PVKB, 'VKB'=VKB, 'PVKF'=PVKF, 'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'PVL'=PVL, 'VL'=VL, 'PVLB'=PVLB, 'VLB'=VLB, 'PVLF'=PVLF, 'VLF'=VLF, 'VLT'=VLT, 'PVbile'=PVbile, 'Vbile'=Vbile, 'PVG'=PVG, 'VG'=VG, 'PVGB'=PVGB, 'VGB'=VGB,
                'PVGF'=PVGF,'VGF'=VGF, 'VGT'=VGT, 'PVGL'=PVGL, 'VGL'=VGL, 'PVM'=PVM, 'VM'=VM, 'PVMB'=PVMB, 'VMB'=VMB, 'PVMF'=PVMF, 'VMF'=VMF, 'VMT'=VMT,
                'PVA'=PVA, 'VA'=VA, 'PVAB'=PVAB, 'VAB'=VAB, 'PVAF'=PVAF, 'VAF'=VAF, 'VAT'=VAT, 'PVR'=PVR, 'VR'=VR, 'PVRB'=PVRB, 'VRB'=VRB, 'PVRF'=PVRF, 'VRF'=VRF, 'VRT'=VRT,
                'PAK'=PAK, 'AK'=AK, 'PAKG'=PAKG, 'AKG'=AKG, 'PAL'=PAL, 'AL'=AL, 'PAG'=PAG, 'AG'=AG, 'PAGL'=PAGL, 'AGL'=AGL, 'PAM'=PAM, 'AM'=AM, 'PAA'=PAA, 'AA'=AA, 'PAR'=PAR, 'AR'=AR,
                'PeffB'=PeffB, 'PeffK'=PeffK, 'PeffL'=PeffL, 'PeffG'=PeffG, 'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'CRssG'=CRssG, 'CRssL'=CRssL, 'CRssK'=CRssK, 
                'Qcardiac'=Qcardiac, 'PQBK'=PQBK, 'QBK'=QBK, 'PQBG'=PQBG, 'QBG'=QBG, 'PQBL'=PQBL, 'QBL'=QBL, 'PQBM'=PQBM, 'QBM'=QBM, 'PQBA'=PQBA, 'QBA'=QBA, 'PQBR'=PQBR, 'QBR'=QBR,
                'Qfeces'=Qfeces, 'PQbile'=PQbile, 'Qbile'=Qbile, 'PQurine'=PQurine, 'PQGFR'=PQGFR, 'QGFR'=QGFR,'Qurine'=Qurine, 'QGFR'=QGFR,
                'CalbB'=CalbB, 'CalbKF'=CalbKF, 'CalbGF'=CalbGF, 'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF,'Ca2uKT'=Ca2uKT,
                'CL_fabpKT'=CL_fabpKT, 'CL_fabpKT1'=CL_fabpKT1, 'CL_fabpKT2'=CL_fabpKT2, 'CL_fabpKT3'=CL_fabpKT3, 'CL_fabpLT'=CL_fabpLT, 'CL_fabpLT1'=CL_fabpLT1, 'CL_fabpLT2'=CL_fabpLT2, 'CL_fabpLT3'=CL_fabpLT3,
                'Ka'=Ka, 'KL_fabp1'=KL_fabp1, 'KL_fabp2'=KL_fabp2, 'KL_fabp3'=KL_fabp3, 'Ka2u'=Ka2u, 'koff'=koff, 'kon'=kon, 'kL_fabpon1'=kL_fabpon1, 'kL_fabpon2'=kL_fabpon2, 'kL_fabpon3'=kL_fabpon3, 'kK_fabpon'=kK_fabpon,
                'kBKF'=kBKF, 'kBF'=kBF, 'kKFKT'=kKFKT, 'kFKT'=kFKT, 'kBLF'=kBLF, 'kLFLT'=kLFLT, 'kBGF'=kBGF, 'kGFGT'=kGFGT, 'kGLGT'=kGLGT, 'kBMF'=kBMF,
                'kMFMT'=kMFMT, 'kBAF'=kBAF, 'kAFAT'=kAFAT, 'kBRF'=kBRF, 'kRFRT'=kRFRT,'kbileLT'=kbileLT,
                'bBKF'=bBKF, 'bKFB'=bKFB, 'bKFKT'=bKFKT, 'bKTKF'=bKTKF, 'bFKT'=bFKT, 'bKTF'=bKTF, 'bBF'=bBF, 'bFB'=bFB, 'bBLF'=bBLF, 'bLFB'=bLFB,
                'bLFLT'=bLFLT, 'bLTLF'=bLTLF, 'bbileLT'=bbileLT, 'bBAF'=bBAF, 'bAFB'=bAFB,
                'bAFAT'=bAFAT, 'bATAF'=bATAF, 'bBRF'=bBRF, 'bRFB'=bRFB, 'bRFRT'=bRFRT, 'bRTRF'=bRTRF, 'bLTbile'=bLTbile, 'bBGF'=bBGF, 'bGFB'=bGFB, 'bGFGT'=bGFGT,
                'bGTGF'=bGTGF, 'bGLGT'=bGLGT, 'bGTGL'=bGTGL, 'bBMF'=bBMF, 'bMFB'=bMFB, 'bMFMT'=bMFMT, 'bMTMF'=bMTMF, 
                'Pbclear'=Pbclear, 'bclear'=bclear, 'Pbreab'=Pbreab, 'breab'=breab, 'Pbabs'=Pbabs, 'babs'=babs, 'Pbefflux'=Pbefflux, 'befflux'=befflux,
                'MalbB'=MalbB, 'MalbKF'=MalbKF, 'ML_fabpKT1'=ML_fabpKT1, 'ML_fabpKT2'=ML_fabpKT2, 'ML_fabpKT3'=ML_fabpKT3, 'MK_fabpKT'=MK_fabpKT, 'MalbLF'=MalbLF, 'ML_fabpLT1'=ML_fabpLT1,
                'ML_fabpLT2'=ML_fabpLT2, 'ML_fabpLT3'=ML_fabpLT3, 'MalbGF'=MalbGF, 'MalbMF'=MalbMF, 'MalbAF'=MalbAF, 'MalbRF'=MalbRF
                
                
                
                
    ))
    
    
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    #==========================Protein binding=============================    
    
    #MBb: Mass of PFOA in blood bound to albumin, MBf: Mass of PFOA in blood not bound to proteins
    dMalbB = koff * MBb - kon * MalbB * MBf / (VB + VLB + VKB + VGB + VMB + VAB + VRB)
    CalbB <- MalbB / (VB + VLB + VKB + VGB + VMB + VAB + VRB)
    
    #MKFb: Mass of PFOA in interstitial fluid of kidney bound to albumin, MKFf: Mass of PFOA in blood not bound to proteins
    dMalbKF = koff * MKFb - kon * MalbKF * MKFf / VKF
    CalbKF <- MalbKF / VKF
    
    #MKTb1: Mass of PFOA in kidney tissue bound to LFABP1, MKTf: Mass of PFOA in kidney tissue not bound to proteins
    dML_fabpKT1 = koff * MKTb1 - kL_fabpon1 * ML_fabpKT1 * MKTf / VKT
    CL_fabpKT1 <- ML_fabpKT1 / VKT
    
    #MKTb2: Mass of PFOA in kidney tissue bound to LFABP2, MKTf: Mass of PFOA in kidney tissue not bound to proteins
    dML_fabpKT2 = koff * MKTb2 - kL_fabpon2 * ML_fabpKT2 * MKTf / VKT
    CL_fabpKT2 <- ML_fabpKT2 / VKT
    
    #MKTb3: Mass of PFOA in kidney tissue bound to LFABP3, MKTf: Mass of PFOA in kidney tissue not bound to proteins
    dML_fabpKT3 = koff * MKTb3 - kL_fabpon3 * ML_fabpKT3 * MKTf / VKT
    CL_fabpKT3 <- ML_fabpKT3 / VKT
    
    #MKTa2b: Mass of PFOA in kidney tissue bound to alpha2mu-globulin, MKTf: Mass of PFOA in kidney tissue not bound to proteins
    dMK_fabpKT = koff * MKTa2b - kK_fabpon * MK_fabpKT * MKTf / VKT
    Ca2uKT <- MK_fabpKT / VKT
    
    #MLFb: Mass of PFOA in interstitial fluid of liver bound to albumin, MLFf: Mass of PFOA in interstitial fluid of liver not bound to proteins
    dMalbLF = koff * MLFb - kon * MalbLF * MLFf / VLF
    CalbLF <- MalbLF / VLF
    
    #MLTb1: Mass of PFOA in liver tissue bound to LFABP1, MLTf: Mass of PFOA in liver tissue not bound to proteins
    dML_fabpLT1 = koff * MLTb1 - kL_fabpon1 * ML_fabpLT1 * MLTf / VLT
    CL_fabpLT1 <- ML_fabpLT1 / VLT
    
    #MLTb2: Mass of PFOA in liver tissue bound to LFABP2, MLTf: Mass of PFOA in liver tissue not bound to proteins
    dML_fabpLT2 = koff * MLTb2 - kL_fabpon2 * ML_fabpLT2 * MLTf / VLT
    CL_fabpLT2 <- ML_fabpLT2 / VLT
    
    #MLTb3: Mass of PFOA in liver tissue bound to LFABP3, MLTf: Mass of PFOA in liver tissue not bound to proteins
    dML_fabpLT3 = koff * MLTb3 - kL_fabpon3 * ML_fabpLT3 * MLTf / VLT
    CL_fabpLT3 <- ML_fabpLT3 / VLT
    
    #MGFb: Mass of PFOA in interstitial fluid of gut bound to albumin, MGFf: Mass of PFOA in interstitial fluid of gut not bound to proteins
    dMalbGF = koff * MGFb - kon * MalbGF * MGFf / VGF
    CalbGF <- MalbGF / VGF
    
    #MMFb: Mass of PFOA in interstitial fluid of muscle bound to albumin, MMFf: Mass of PFOA in interstitial fluid of muscle not bound to proteins
    dMalbMF = koff * MMFb - kon * MalbMF * MMFf / VMF
    CalbMF <- MalbMF / VMF
    
    #MAFb: Mass of PFOA in interstitial fluid of adipose bound to albumin, MAFf: Mass of PFOA in interstitial fluid of adipose not bound to proteins
    dMalbAF = koff * MAFb - kon * MalbAF * MAFf / VAF
    CalbAF <- MalbAF / VAF
    
    #MRFb: Mass of PFOA in interstitial fluid of rest of body bound to albumin, MAFf: Mass of PFOA in interstitial fluid of rest of body not bound to proteins
    dMalbRF = koff * MRFb - kon * MalbRF * MRFf / VRF
    CalbRF <- MalbRF / VRF
    
    
    #first-order rate constants for passive diffusion and active transport between subcompartments
    
    bBon <- CalbB * kon
    bBoff <- koff
    
    bKFon <- CalbKF * kon
    bKFoff <- koff
    
    bKTon1 <- CL_fabpKT1 * kL_fabpon1
    bKTon2 <- CL_fabpKT2 * kL_fabpon2
    bKTon3 <- CL_fabpKT3 * kL_fabpon3
    bKToff <- koff
    
    bKTa2on <- Ca2uKT * kK_fabpon
    bKTa2off <- koff
    
    bLFon <- CalbLF * kon
    bLFoff <- koff
    
    bLFon1 <- CL_fabpLT1 * kL_fabpon1
    bLFon2 <- CL_fabpLT2 * kL_fabpon2
    bLFon3 <- CL_fabpLT3 * kL_fabpon3
    bLFoff <- koff
    
    bGFon <- CalbGF * kon
    bGFoff <- koff
    
    bMFon <- CalbMF * kon
    bMFoff <- koff
    
    bAFon <- CalbAF * kon
    bAFoff <- koff
    
    bRFon <- CalbRF * kon
    bRFoff <- koff
    
    #====================PFOA mass balance at each tissue or fluid compartment==============================     
    
    #Blood subcompartment
    
    dMBf = (bKFB*MKFf + bLFB*MLFf + bGFB*MGFf + bMFB*MMFf + bAFB*MAFf + bRFB*MRFf) - (bBKF*MBf + bBLF*MBf + bBGF*MBf + bBMF*MBf + bBAF*MBf + bBRF*MBf) + bFB*MFilf - bBF*MBf + bBoff*MBb - bBon*MBf
    
    dMBb = bBon*MBf - bBoff*MBb #PFOA in blood bound to albumin
    
    #==================================================  Interstitial fluid
    
    #Kidney interstitial fluid subcompartment (bclear = bKoat1 + bKoat3)
    
    dMKFf = bBKF*MBf - bKFB*MKFf + bKTKF*MKTf - bKFKT*MKFf + befflux*MKTf - bclear*MKFf + bKFoff*MKFb - bKFon*MKFf
    
    dMKFb = bKFon*MKFf - bKFoff*MKFb #Mass of PFOA in interstitial fluid of kidney bound to albumin
    
    #Liver interstitial fluid subcompartment 
    
    dMLFf = bBLF*MBf - bLFB*MLFf + bLTLF*MLTf - bLFLT*MLFf - babs*MLFf + bLFoff*MLFb - bLFon*MLFf
    
    dMLFb = bLFon*MLFf - bLFoff*MLFb #Mass of PFOA in interstitial fluid of liver bound to albumin
    
    #Gut interstitial fluid subcompartment 
    
    dMGFf = bBGF*MBf - bGFB*MGFf + bGTGF*MGTf - bGFGT*MGFf + bGFoff*MGFb - bGFon*MGFf
    
    dMGFb = bGFon*MGFf - bGFoff*MGFb #Mass of PFOA in interstitial fluid of gut bound to albumin
    
    #Muscle interstitial fluid subcompartment 
    
    dMMFf = bBMF*MBf - bMFB*MMFf + bMTMF*MMTf - bMFMT*MMFf + bMFoff*MMFb - bMFon*MMFf
    
    dMMFb = bMFon*MMFf - bMFoff*MMFb #Mass of PFOA in interstitial fluid of muscle bound to albumin
    
    #Adipose interstitial fluid subcompartment 
    
    dMAFf = bBAF*MBf - bAFB*MAFf + bATAF*MATf - bAFAT*MAFf + bAFoff*MAFb - bAFon*MAFf
    
    dMAFb = bAFon*MAFf - bAFoff*MAFb #Mass of PFOA in interstitial fluid of adipose bound to albumin
    
    #Rest of body interstitial fluid subcompartment 
    
    dMRFf = bBRF*MBf - bRFB*MRFf + bRTRF*MRTf - bRFRT*MRFf + bRFoff*MRFb - bRFon*MRFf
    
    dMRFb = bRFon*MRFf - bRFoff*MRFb #Mass of PFOA in interstitial fluid of rest of body bound to albumin
    
    #==================================================  Tissue
    
    #Adipose tissue subcompartment 
    
    dMATf = bAFAT*MAFf - bATAF*MATf
    
    #Muscle tissue subcompartment 
    
    dMMTf = bMFMT*MMFf - bMTMF*MMTf
    
    #Rest of body tissue subcompartment 
    
    dMRTf = bRFRT*MRFf - bRTRF*MRTf
    
    #Kidney tissue subcompartment
    
    dMKTf = bKFKT*MKFf - bKTKF*MKTf + bclear*MKFf - befflux*MKTf + bFKT*MFilf - bKTF*MKTf + breab*MFilf - (bKTon1*MKTf - bKToff*MKTb1) - (bKTon2*MKTf - bKToff*MKTb2) - (bKTon3*MKTf - bKToff*MKTb3) - (bKTa2on*MKTf - bKTa2off*MKTa2b)
    
    dMKTb1 = bKTon1*MKTf - bKToff*MKTb1 #PFOA in kidney tissue bound to LFABP1
    dMKTb2 = bKTon2*MKTf - bKToff*MKTb2 #PFOA in kidney tissue bound to LFABP2
    dMKTb3 = bKTon3*MKTf - bKToff*MKTb3 #PFOA in kidney tissue bound to LFABP3
    dMKTa2b = bKTa2on*MKTf - bKTa2off*MKTa2b #PFOA in kidney tissue bound to alpha2mu-globulin
    
    dMFilf = bBF*MBf - bFB*MFilf + bKTF*MKTf - bFKT*MFilf - breab*MFilf - (Qurine/VFil)*MFilf
    dMurine = (Qurine/VFil)*MFilf
    
    #Liver tissue subcompartment
    
    dMLTf = bLFLT*MLFf - bLTLF*MLTf + babs*MLFf + bbileLT*Mbilef - bLTbile*MLTf - (bLFon1*MLTf - bLFoff*MLTb1) - (bLFon2*MLTf - bLFoff*MLTb2) - (bLFon3*MLTf - bLFoff*MLTb3)
    
    dMLTb1 = bLFon1*MLTf - bLFoff*MLTb1 #PFOA in liver tissue bound to LFABP1
    dMLTb2 = bLFon2*MLTf - bLFoff*MLTb2 #PFOA in liver tissue bound to LFABP2
    dMLTb3 = bLFon3*MLTf - bLFoff*MLTb3 #PFOA in liver tissue bound to LFABP3
    
    dMbilef = bLTbile*MLTf - bbileLT*Mbilef - (Qbile/Vbile)*Mbilef
    #dMbile = (Qbile/Vbile)*Mbilef 
    
    #Gut tissue subcompartment
    
    dMGTf = bGFGT*MGFf - bGTGF*MGTf + bGLGT*MGLf - bGTGL*MGTf
    dMGLf = bGTGL*MGTf - bGLGT*MGLf + (Qbile/Vbile)*Mbilef - (Qfeces/VGL)*MGLf
    dMfeces = (Qfeces/VGL)*MGLf
    
    Cblood <- (MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VB /Vplasma * 10^6
    Ckidney <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VKB+MKFf+MKFb+MKFf+MKTb1+MKTb2+MKTb3+MKTa2b) / (VKB+VKT+VKF) * 10^6
    Cliver <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VLB+MLFf+MLFb+MLTf+MLTb1+MLTb2+MLTb3) / (VLB+VLT+VLF) * 10^6
    Cgut <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VGB+MGFf+MGFb+MGTf) / (VGB+VGT+VGF) * 10^6
    Cmuscle <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VMB+MMFf+MMFb+MMTf) / (VMB+VMT+VMF) * 10^6
    Cadipose <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VAB+MAFf+MAFb+MATf) / (VAB+VAT+VAF) * 10^6
    Crest <- ((MBf + MBb) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VRB+MRFf+MRFb+MRTf) / (VRB+VRT+VRF) * 10^6
    Cfeces <- MGLf / VGL * 10^6
    Cbile <- Mbilef / Vbile * 10^6
    Curine <- MFilf / VFil * 10^6
    
    return(list(c('dMBf'=dMBf, 'dMBb'=dMBb, 'dMKFf'=dMKFf, 'dMKFb'=dMKFb,'dMLFf'=dMLFf, 'dMLFb'=dMLFb, 'dMGFf'=dMGFf, 'dMGFb'=dMGFb, 'dMMFf'=dMMFf, 'dMMFb'=dMMFb, 'dMAFf'=dMAFf, 'dMAFb'=dMAFb,
                  'dMRFf'=dMRFf, 'dMRFb'=dMRFb, 'dMATf'=dMATf, 'dMMTf'=dMMTf, 'dMRTf'=dMRTf, 'dMKTf'=dMKTf, 'dMKTb1'=dMKTb1, 'dMKTb2'=dMKTb2, 'dMKTb3'=dMKTb3, 'dMKTa2b'=dMKTa2b, 'dMFilf'=dMFilf, 'dMurine'=dMurine,
                  'dMLTf'=dMLTf, 'dMLTb1'=dMLTb1, 'dMLTb2'=dMLTb2, 'dMLTb3'=dMLTb3, 'dMbilef'=dMbilef, #'dMbile'=dMbile, 
                  'dMGTf'=dMGTf, 'dMGLf'=dMGLf, 'dMfeces'=dMfeces, 'dMalbB'=dMalbB, 'dMalbKF'=dMalbKF, 'dML_fabpKT1'=dML_fabpKT1,
                  'dML_fabpKT2'=dML_fabpKT2, 'dML_fabpKT3'=dML_fabpKT3, 'dMK_fabpKT'=dMK_fabpKT, 'dMalbLF'=dMalbLF, 'dML_fabpLT1'=dML_fabpLT1,
                  'dML_fabpLT2'=dML_fabpLT2, 'dML_fabpLT3'=dML_fabpLT3, 'dMalbGF'=dMalbGF, 'dMalbMF'=dMalbMF, 'dMalbAF'=dMalbAF, 'dMalbRF'=dMalbRF),   
                  'Cblood'=Cblood, 'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cgut'=Cgut, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose, 'Crest'=Crest,
                  'Cfeces'=Cfeces, 'Cbile'=Cbile, 'Curine'=Curine 
               
                
                ))
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBf <- admin.dose *BW
    MBb <- 0; MKFf <- 0; MKFb <- 0; MLFf <- 0; MLFb <- 0; MGFf <- 0; MGFb <- 0; MMFf <- 0; MMFb <- 0; MAFf <- 0; MAFb <- 0; MRFf <- 0; MRFb <- 0;
    MATf <- 0; MMTf <- 0; MRTf <- 0; MKTf <- 0; MKTb1 <- 0; MKTb2 <- 0; MKTb3 <- 0; MKTa2b <-0; 
    MFilf <- 0; Murine <- 0; MLTf <- 0; MLTb1 <- 0; MLTb2 <- 0; MLTb3 <- 0; Mbilef <- 0; Mbile <- 0; MGTf <- 0;  MGLf <- 0; Mfeces <-0
    MalbB <- MalbB; MalbKF <- MalbKF; ML_fabpKT1 <- ML_fabpKT1; ML_fabpKT2 <- ML_fabpKT2; ML_fabpKT3 <- ML_fabpKT3; MK_fabpKT <- MK_fabpKT;
    MalbLF <- MalbLF; ML_fabpLT1 <- ML_fabpLT1; ML_fabpLT2 <- ML_fabpLT2; ML_fabpLT3 <- ML_fabpLT3; MalbGF <- MalbGF; MalbMF <- MalbMF; MalbAF <- MalbAF; MalbRF <- MalbRF
    
    
    
    
    
    return(c('MBf'=MBf, 'MBb'=MBb, 'MKFf'=MKFf, 'MKFb'=MKFb,'MLFf'=MLFf, 'MLFb'=MLFb, 'MGFf'=MGFf, 'MGFb'=MGFb, 'MMFf'=MMFf, 'MMFb'=MMFb, 'MAFf'=MAFf, 'MAFb'=MAFb,
             'MRFf'=MRFf, 'MRFb'=MRFb, 'MATf'=MATf, 'MMTf'=MMTf, 'MRTf'=MRTf, 'MKTf'=MKTf, 'MKTb1'=MKTb1, 'MKTb2'=MKTb2, 'MKTb3'=MKTb3, 'MKTa2b'=MKTa2b, 'MFilf'=MFilf, 'Murine'=Murine,
             'MLTf'=MLTf, 'MLTb1'=MLTb1, 'MLTb2'=MLTb2, 'MLTb3'=MLTb3, 'Mbilef'=Mbilef, #'Mbile'=Mbile, 
             'MGTf'=MGTf, 'MGLf'=MGLf, 'Mfeces'=Mfeces, 'MalbB'=MalbB, 'MalbKF'=MalbKF, 'ML_fabpKT1'=ML_fabpKT1,
             'ML_fabpKT2'=ML_fabpKT2, 'ML_fabpKT3'=ML_fabpKT3, 'MK_fabpKT'=MK_fabpKT, 'MalbLF'=MalbLF, 'ML_fabpLT1'=ML_fabpLT1,
             'ML_fabpLT2'=ML_fabpLT2, 'ML_fabpLT3'=ML_fabpLT3, 'MalbGF'=MalbGF, 'MalbMF'=MalbMF, 'MalbAF'=MalbAF, 'MalbRF'=MalbRF
    ))
    
    
  })
}


################################################################################

BW <- 0.244  # body weight (kg)
admin.dose <- c(1 * 1e-6) # administered dose in kg PFOA/kg BW
admin.time <- c(0) # time when doses are administered, in days
user_input <- list('BW'=BW,
                   "admin.dose"=admin.dose)

params <- create.params(user_input)
inits <- create.inits(params)

# 1 mg/kg IV
# time_points <- unlist(as.vector( read.delim('time_points2.txt', header = FALSE))) * (24*3600)
# extra_time_points <- seq(0, 1*24*3600, 5)
# 
# sample_time <- sort(unique(c(time_points,extra_time_points)))

#sample_time=seq(0,22*24*3600, length.out= 19008000) # / (24*3600)
# sample_time=seq(0,22*24*3600,10)

#sample_time = seq(0, 300*60, 60)

#<<<<<<< HEAD
sample_time=seq(0,22*24*3600,10)
#=======
#>>>>>>> 477e3df2b7fb740a0303318016b4a282d4dd98f7
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
solution$time = solution$time/ (24*3600)

rowSums(solution[,2:33])


# predictions <- solution[solution$time %in% (time_points / (24*3600)),]
#write.csv(predictions, file='predictions_05.csv')
######################################################################################


# Plot with ggplot2
library (ggplot2)
compartment <- c('Cblood')
color_codes <- scales::hue_pal()(length(compartment))

plot <- ggplot(data = solution)+
  geom_line( aes(x = time, y = Cblood, color='Cblood'), size=1.3)+
  scale_y_log10()+       
  scale_x_continuous()+
 
  labs(title = 'Predicted values',
       y = 'Concentration (ng/g)', x = "Time (days)")+
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25),
        legend.text=element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
  
  scale_color_manual("Compartment", values=color_codes)+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))
print(plot)
