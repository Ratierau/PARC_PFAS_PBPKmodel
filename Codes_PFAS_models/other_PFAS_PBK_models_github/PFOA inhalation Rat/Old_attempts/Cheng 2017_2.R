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

    CalbB <- 281e-3 #*7.8 #n=7.8 binding sites 
    CalbKF <- 243e-3 #*7.8 #n=7.8 binding sites 
    CalbLF <- 243e-3 #*7.8 #n=7.8 binding sites 
    CalbGF <- 146e-3 #*7.8 #n=7.8 binding sites 
    CalbMF <- 146e-3 #*7.8 #n=7.8 binding sites 
    CalbAF <- 73e-3 #*7.8 #n=7.8 binding sites 
    CalbRF <- 73e-3 #*7.8 #n=7.8 binding sites 
    
    #Alpha2mu-globulin concentration in kidney tissue (mol/m^3)

    Ca2uKT <- 110e-3

    #LFABP concentration in kidney and liver tissue (mol/m^3)

    CLfabpKT <- 2.65e-3 #*3 #n=3 binding sites (Table S3)
    CLfabpLT <- 133e-3 #*3 #n=3 binding sites (Table S3)
  

    #======Table S2=======#
    #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
    #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details

    Ka <-  24.18    # 3.1*7.8 m3/mol multiplying by number of binding sites
    KLfabp <- 135  # 45.0*3 geo_mean of 3 binding affinity, may try normal mean
    Ka2u <- 0.5

    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    kBKF <- ((1/QBK) + 1/(PeffB * AK))^(-1)
    kBF <- PeffB * AKG #AKG-> the surface area of glomerular capillary
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

    # MalbB <- CalbB * (VB+VLB+VKB+VGB+VMB+VAB+VRB)
    # MalbKF <- CalbKF * VKF
    # MalbLF <- CalbLF * VLF
    # MalbGF <- CalbGF * VGF
    # MalbMF <- CalbMF * VMF
    # MalbAF <- CalbAF * VAF
    # MalbRF <- CalbRF * VRF
    
    
    
    
    return(list('admin.dose'=admin.dose,'PVB'=PVB, 'VB'=VB, 'PVplsma'=PVplasma, 'Vplasma'=Vplasma, 'PVK'=PVK, 'VK'=VK, 'PVKB'=PVKB, 'VKB'=VKB, 'PVKF'=PVKF, 'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'PVL'=PVL, 'VL'=VL, 'PVLB'=PVLB, 'VLB'=VLB, 'PVLF'=PVLF, 'VLF'=VLF, 'VLT'=VLT, 'PVbile'=PVbile, 'Vbile'=Vbile, 'PVG'=PVG, 'VG'=VG, 'PVGB'=PVGB, 'VGB'=VGB,
                'PVGF'=PVGF,'VGF'=VGF, 'VGT'=VGT, 'PVGL'=PVGL, 'VGL'=VGL, 'PVM'=PVM, 'VM'=VM, 'PVMB'=PVMB, 'VMB'=VMB, 'PVMF'=PVMF, 'VMF'=VMF, 'VMT'=VMT,
                'PVA'=PVA, 'VA'=VA, 'PVAB'=PVAB, 'VAB'=VAB, 'PVAF'=PVAF, 'VAF'=VAF, 'VAT'=VAT, 'PVR'=PVR, 'VR'=VR, 'PVRB'=PVRB, 'VRB'=VRB, 'PVRF'=PVRF, 'VRF'=VRF, 'VRT'=VRT,
                'PAK'=PAK, 'AK'=AK, 'PAKG'=PAKG, 'AKG'=AKG, 'PAL'=PAL, 'AL'=AL, 'PAG'=PAG, 'AG'=AG, 'PAGL'=PAGL, 'AGL'=AGL, 'PAM'=PAM, 'AM'=AM, 'PAA'=PAA, 'AA'=AA, 'PAR'=PAR, 'AR'=AR,
                
                'PeffB'=PeffB, 'PeffK'=PeffK, 'PeffL'=PeffL, 'PeffG'=PeffG, 'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'CRssG'=CRssG, 'CRssL'=CRssL, 'CRssK'=CRssK, 
                
                'Qcardiac'=Qcardiac, 'PQBK'=PQBK, 'QBK'=QBK, 'PQBG'=PQBG, 'QBG'=QBG, 'PQBL'=PQBL, 'QBL'=QBL, 'PQBM'=PQBM, 'QBM'=QBM, 'PQBA'=PQBA, 'QBA'=QBA, 'PQBR'=PQBR, 'QBR'=QBR,
                'Qfeces'=Qfeces, 'PQbile'=PQbile, 'Qbile'=Qbile, 'PQurine'=PQurine, 'PQGFR'=PQGFR, 'QGFR'=QGFR,'Qurine'=Qurine, 'QGFR'=QGFR,
                
                'CalbB'=CalbB, 'CalbKF'=CalbKF, 'CalbLF'=CalbLF, 'CalbGF'=CalbGF, 'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF,'Ca2uKT'=Ca2uKT,
                'CLfabpKT'=CLfabpKT,'CLfabpLT'=CLfabpLT,
                
                'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,'kBKF'=kBKF, 'kBF'=kBF, 'kKFKT'=kKFKT, 'kFKT'=kFKT, 'kBLF'=kBLF, 'kLFLT'=kLFLT, 'kBGF'=kBGF, 'kGFGT'=kGFGT, 'kGLGT'=kGLGT, 'kBMF'=kBMF,
                'kMFMT'=kMFMT, 'kBAF'=kBAF, 'kAFAT'=kAFAT, 'kBRF'=kBRF, 'kRFRT'=kRFRT,'kbileLT'=kbileLT,
                
                'bBKF'=bBKF, 'bKFB'=bKFB, 'bKFKT'=bKFKT, 'bKTKF'=bKTKF, 'bFKT'=bFKT, 'bKTF'=bKTF, 'bBF'=bBF, 'bFB'=bFB, 'bBLF'=bBLF, 'bLFB'=bLFB,
                'bLFLT'=bLFLT, 'bLTLF'=bLTLF, 'bbileLT'=bbileLT, 'bBAF'=bBAF, 'bAFB'=bAFB,
                'bAFAT'=bAFAT, 'bATAF'=bATAF, 'bBRF'=bBRF, 'bRFB'=bRFB, 'bRFRT'=bRFRT, 'bRTRF'=bRTRF, 'bLTbile'=bLTbile, 'bBGF'=bBGF, 'bGFB'=bGFB, 'bGFGT'=bGFGT,
                'bGTGF'=bGTGF, 'bGLGT'=bGLGT, 'bGTGL'=bGTGL, 'bBMF'=bBMF, 'bMFB'=bMFB, 'bMFMT'=bMFMT, 'bMTMF'=bMTMF, 
                'Pbclear'=Pbclear, 'bclear'=bclear, 'Pbreab'=Pbreab, 'breab'=breab, 'Pbabs'=Pbabs, 'babs'=babs, 'Pbefflux'=Pbefflux, 'befflux'=befflux
                
                
                
                
                
    ))
    
    
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
   
    #====================PFOA mass balance at each tissue or fluid compartment==============================     
    
    #Blood subcompartment
    
    dMB = (bKFB*MKFf + bLFB*MLFf + bGFB*MGFf + bMFB*MMFf + bAFB*MAFf + bRFB*MRFf) - (bBKF*MBf + bBLF*MBf + bBGF*MBf + bBMF*MBf + bBAF*MBf + bBRF*MBf) + bFB*MFil - bBF*MBf 
    
  
      #==================================================  Interstitial fluid
    
    #Kidney interstitial fluid subcompartment (bclear = bKoat1 + bKoat3)
    
    dMKF = bBKF*MBf - bKFB*MKFf + bKTKF*MKTf - bKFKT*MKFf + befflux*MKTf - bclear*MKFf 
    
    
    #Liver interstitial fluid subcompartment 
    
    dMLF = bBLF*MBf - bLFB*MLFf + bLTLF*MLTf - bLFLT*MLFf - babs*MLFf 
    
    
    #Gut interstitial fluid subcompartment 
    
    dMGF = bBGF*MBf - bGFB*MGFf + bGTGF*MGT - bGFGT*MGFf 
    
    
    #Muscle interstitial fluid subcompartment 
    
    dMMF = bBMF*MBf - bMFB*MMFf + bMTMF*MMT - bMFMT*MMFf 
    
    
    #Adipose interstitial fluid subcompartment 
    
    dMAF = bBAF*MBf - bAFB*MAFf + bATAF*MAT - bAFAT*MAFf 
    
    
    #Rest of body interstitial fluid subcompartment 
    
    dMRF = bBRF*MBf - bRFB*MRFf + bRTRF*MRT - bRFRT*MRFf 
    
    
    #==================================================  Tissue
    
    #Lung Tissue subcompartment
    
    #dLNf = bALFLN*MALFf - bLNALF*MLNf + bBLN*MBf - bLNB*MLNf
    
    #Adipose tissue subcompartment 
    
    dMAT = bAFAT*MAFf - bATAF*MAT
    
    #Muscle tissue subcompartment 
    
    dMMT = bMFMT*MMFf - bMTMF*MMT
    
    #Rest of body tissue subcompartment 
    
    dMRT = bRFRT*MRFf - bRTRF*MRT
    
    #Kidney tissue subcompartment
    
    dMKT = bKFKT*MKFf - bKTKF*MKTf + bclear*MKFf - befflux*MKTf + bFKT*MFil - bKTF*MKTf + breab*MFil 
    
    
    dMFil = bBF*MBf - bFB*MFil + bKTF*MKTf - bFKT*MFil - breab*MFil - (Qurine/VFil)*MFil
    dMurine = (Qurine/VFil)*MFil
    
    #Liver tissue subcompartment
    
    dMLT = bLFLT*MLFf - bLTLF*MLTf + babs*MLFf + bbileLT*Mbile - bLTbile*MLTf 
    
   
    dMbile = bLTbile*MLTf - bbileLT*Mbile - (Qbile/Vbile)*Mbile
  
    #Gut tissue subcompartment
    
    dMGT = bGFGT*MGFf - bGTGF*MGT + bGLGT*MGL - bGTGL*MGT
    dMGL = bGTGL*MGT - bGLGT*MGL + (Qbile/Vbile)*Mbile - (Qfeces/VGL)*MGL
    dMfeces = (Qfeces/VGL)*MGL
    
    #Mfree calculation using the expression of free fraction ff
    
    MBf = MB * 1.0 / (1.0 + CalbB * Ka)
    MKFf = MKF * 1.0 / (1.0 + CalbKF * Ka)
    MLFf = MLF * 1.0 / (1.0 + CalbLF * Ka)
    MGFf = MGF * 1.0 / (1.0 + CalbGF * Ka)
    MMFf = MMF * 1.0 / (1.0 + CalbMF * Ka)
    MAFf = MAF * 1.0 / (1.0 + CalbAF * Ka)
    MRFf = MRF * 1.0 / (1.0 + CalbRF * Ka)
    MKTf = MKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp)
    MLTf = MLT * 1.0 / (1.0 + CLfabpLT * KLfabp)
    
   
    Cblood <- MB /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VB /Vplasma * 10^6
    Ckidney <- (MB /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VKB+MKF+MKT) / (VKB+VKT+VKF) * 10^6
    Cliver <- ((MB) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VLB+MLF+MLT) / (VLB+VLT+VLF) * 10^6
    Cgut <- ((MB) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VGB+MGF+MGT) / (VGB+VGT+VGF) * 10^6
    Cmuscle <- ((MB) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VMB+MMF+MMT) / (VMB+VMT+VMF) * 10^6
    Cadipose <- ((MB) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VAB+MAF+MAT) / (VAB+VAT+VAF) * 10^6
    Crest <- ((MB) /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VRB+MRF+MRT) / (VRB+VRT+VRF) * 10^6
    Cfeces <- MGL / VGL * 10^6
    Cbile <- Mbile / Vbile * 10^6
    Curine <- MFil / VFil * 10^6
    
    return(list(c('dMB'=dMB,'dMKF'=dMKF,'dMLF'=dMLF,'dMGF'=dMGF,'dMMF'=dMMF,'dMAF'=dMAF, 'dMRF'=dMRF,  
                  'dMAT'=dMAT, 'dMMT'=dMMT, 'dMRT'=dMRT, 'dMKT'=dMKT, 'dMFil'=dMFil, 'dMurine'=dMurine,
                  'dMLT'=dMLT, 'dMbile'=dMbile, 'dMGT'=dMGT, 'dMGL'=dMGL, 'dMfeces'=dMfeces,
                  'MBf'=MBf, 'MKFf'=MKFf, 'MLFf'=MLFf, 'MGFf'=MGFf, 'MMFf'=MMFf, 'MAFf'=MAFf, 'MRFf'=MRFf, 'MKTf'=MKTf, 'MLTf'=MLTf),
                  
                #'CalbB'=CalbB, 'CalbB'=CalbB, 'CalbLF'=CalbLF, 'CalbGF'=CalbGF, 'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF
                  
                  'Cblood'=Cblood, 'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cgut'=Cgut, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose, 'Crest'=Crest,
                  'Cfeces'=Cfeces, 'Cbile'=Cbile, 'Curine'=Curine
               
                
                ))
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MB <- admin.dose *BW; MKF <- 0; MLF <- 0; MGF <- 0; MMF <- 0; MAF <- 0; MRF <- 0; 
    MAT <- 0; MMT <- 0; MRT <- 0; MKT <- 0; MFil <- 0; Murine <- 0;
    MLT <- 0; Mbile <- 0; MGT <- 0; MGL <- 0; Mfeces <-0;
    MBf <- 0; MKFf <-0; MLFf <-0; MGFf <-0; MMFf <-0; MAFf <-0; MRFf <-0; MKTf <-0; MLTf <-0
   
     
    
    return(c('MB'=MB,'MKF'=MKF,'MLF'=MLF, 'MGF'=MGF, 'MMF'=MMF, 'MAF'=MAF, 'MRF'=MRF, 
             'MAT'=MAT, 'MMT'=MMT, 'MRT'=MRT, 'MKT'=MKT, 'MFil'=MFil, 'Murine'=Murine,
             'MLT'=MLT, 'Mbile'=Mbile, 'MGT'=MGT, 'MGL'=MGL, 'Mfeces'=Mfeces,
             'MBf'=MBf, 'MKFf'=MKFf, 'MLFf'=MLFf, 'MGFf'=MGFf, 'MMFf'=MMFf, 'MAFf'=MAFf, 'MRFf'=MRFf, 'MKTf'=MKTf, 'MLTf'=MLTf
            
             
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


sample_time=seq(0,100,1)
#=======

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
solution$time = solution$time/ (24*3600)

rowSums(solution[,2:19])


predictions <- solution[solution$time %in% (time_points / (24*3600)),]
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