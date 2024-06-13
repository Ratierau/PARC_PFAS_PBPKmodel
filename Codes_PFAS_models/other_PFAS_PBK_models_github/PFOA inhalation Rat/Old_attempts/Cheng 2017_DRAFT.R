PVB <- 54e-6 #13.5 mL/0.244 kg=0.055 mL/kg~55e-3 mL/g kg=L)
VB <- PVB * BW #blood volume g=mL
PVplasma <- 31.2e-6 
Vplasma <- PVplasma * BW #plasma volume g=mL
PVK <- 7.3e-6 
VK <- PVK * BW #kidney volume g=mL
PVKB <- 0.16 
VKB <- PVKB * PVK * BW #kidney blood volume g=mL
PVKF <- 0.13 
VKF <- PVKF * PVK * BW #kidney interstitial fluid volume g=mL
VKT <- VK - VKF #kidney tissue volume g=mL
VFil <- 2.5e-7 #renal filtrate volume mL=g
PVL <- 3.66e-5 
VL <- PVL * BW #liver volume g=mL
PVLB <- 0.21 
VLB <- PVLB * PVL * BW #liver blood volume g=mL
PVLF <- 0.049 
VLF <- PVLF * PVL* BW #liver interstitial fluid volume g=mL
VLT <- VL - VLF #liver tissue volume g=mL
PVbile <- 0.004 
Vbile <- PVbile * PVL * BW #bile volume g=mL
PVG <- 2.69e-5 
VG <- PVG * BW #gut volume g=mL
PVGB <- 0.034
VGB <- PVGB * PVG * BW #gut blood volume g=mL
PVGF <- 0.28
VGF <- PVGF * PVG * BW #gut interstitial fluid volume g=mL
VGT <- VG - VGF #gut tissue volume g=mL
PVGL <- 4.5e-5 
VGL <- PVGL * BW #gut lumen volume g=mL
PVM <- 40.43e-5 
VM <- PVM * BW #muscle volume g=mL
PVMB <- 0.04 
VMB <- PVMB * PVM * BW #muscle blood volume g=mL
PVMF <- 0.054
VMF <- PVMF * PVM * BW #muscle interstitial fluid volume g=mL
VMT <- VM - VMF #muscle tissue volume g=mL
PVA <- 7e-2 
VA <- PVA * BW #adipose volume g=mL
PVAB <- 0.02 
VAB <- PVAB * PVA * BW #% adipose blood volume g=mL
PVAF <- 0.174
VAF <- PVAF * PVA * BW #adipose interstitial fluid volume g=mL
VAT <- VA - VAF #adipose tissue volume g=mL
PVR <- 1/1e3 - PVB - PVK - PVL - PVG - PVM - PVA
VR <- PVR * BW #volume of the rest of the body g=mL
PVRB <- 0.036 
VRB <- PVRB * PVR * BW #volume of the blood of the rest of body g=mL
PVRF <- 0.18  
VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body g=mL
VRT <- VR - VRF #tissue volume of the rest of body g=mL


CBf=MBf/(VB+VLB+VKB+VGB+VMB+VAB+VRB)

#Kidney
CKFf=MKFf/VKF; CKTf=MLTf/VKT; CFil=MFil/VFil
#Liver
CLFf=MLFf/VLF; CLTf=MLTf/VLT;Cbile=Mbile/Vbile
#Gut
CGFf=MGFf/vGF; CGT=MGT/VGT; CGL=MGL/VGL
#Muscle
CMFf=MMFf/VMF; CMT=MMT/VMT
#Adipose
CAFf=MAFf/VAF; CAT=MAT/VAT
#Rest of Body
CRFr=MRFr/VRF; CRT=MRT/VRT



#Blood subcompartment

dCB = kBKF*(CKFf-CBf) + kBLF*(CLFf-CBf) + kBGF*(CGFf-MBf) + kBMF*(CMFf-CBf) + kBAF*(CAFf-CBf) + kBRF*(CRFf-CBf) + kBF*(CFil-CBf) 

#Kidney

#interstitial fluid subcompartment
dCKF = kBKF*(CBf-CKFf) + kKFKT*(CKTf-CKFf) + (VmK_Osta*CKTf/KmK_Osta+CKTf) - (VmK_Oat1*CKFf/KmK_Oat1+CKFf) - (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
#Kidney proximal tubule cells subcompartment
dCKT = kKFKT*(CKFf-CKTf) + kFKT*(CFil-CKTf) + (VmK_Oatp*CFil/KmK_Oatp+CFil) - (VmK_Osta*CKTf/KmK_Osta+CKTf) + (VmK_Oat1*CKFf/KmK_Oat1+CKFf) + (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
dCFil = kBF*(CBf-CFil) + kFKT*(CKTf-CFil) - (VmK_Oatp*CFil/KmK_Oatp+CFil) - (Qurine/VFil)*CFil
dCurine = (Qurine/VFil)*CFil

#Liver

#interstitial fluid subcompartment 
dCLF = kBLF*(CBf-CLFf) + kLFLT*(CLTf-CLFf) - (VmL_Oatp*CLFf/KmL_Oatp+CLFf) - (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf)
#Liver tissue subcompartment
dCLT = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/KmL_Oatp+CLFf) + (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf) + kbileLT*(Cbile-CLTf)
dCbile = kbileLT*(CLTf-Cbile) - (Qbile/Vbile)*Cbile


#Gut

#interstitial fluid subcompartment 
dCGF = kBGF*(CBf-CGFf) + kGFGT*(CGT-CGFf) 
#Gut tissue subcompartment
dCGT = kGFGT*(CGFf-CGT) + kGLGT*(CGL-CGT)
dCGL = kGLGT*(CGT-CGL) + (Qbile/Vbile)*Cbile - (Qfeces/VGL)*CGL
dCfeces = (Qfeces/VGL)*CGL

#Muscle

#interstitial fluid subcompartment 
dCMF = kBMF*(CBf-CMFf) + kMFMT*(CMT-CMFf) 
#Muscle tissue subcompartment 
dCMT = kMFMT*(CMFf-CMT)

#Adipose

#interstitial fluid subcompartment 
dCAF = kBAF*(CBf-CAFf) + kAFAT*(CAT-CAFf) 
#Adipose tissue subcompartment 
dCAT = kAFAT*(CAFf-CAT)

#Rest of body

#interstitial fluid subcompartment 
dCRF = kBRF*(CBf-CRFf) + kRFRT*(CRT-CRFf) 
#Rest of body tissue subcompartment 
dCRT = kRFRT*(CRFf-CRT)

Cblood <- CBf * VB /Vplasma * 1e-9
Ckidney <- CBf * (VKB+MKFf+MKTf) / (VKB+VKT+VKF) * 1e-9
Cliver <- CBf * (VLB+MLFf+MLTf) / (VLB+VLT+VLF) * 1e-9
Cgut <- CBf * (VGB+MGFf+MGT) / (VGB+VGT+VGF) * 1e-9
Cmuscle <- CBf * (VMB+MMFf+MMT) / (VMB+VMT+VMF) * 1e-9
Cadipose <- CBf * (VAB+MAFf+MAT) / (VAB+VAT+VAF) * 1e-9
Crest <- CBf * (VRB+MRFf+MRT) / (VRB+VRT+VRF) * 1e-9
Cfeces <- CGL * 1e-9
Cbile <- Cbile * 1e-9
Curine <- CFil * 1e-9




#Blood subcompartment

dMB = kBKF*(MKFf-MBf) + kBLF*(MLFf-MBf) + kBGF*(MGFf-MBf) + kBMF*(MMFf-MBf) + kBAF*(MAFf-MBf) + kBRF*(MRFf-MBf) + kBF*(MFil-MBf) 

#Kidney

#interstitial fluid subcompartment
dMKF = kBKF*(MBf-MKFf) + kKFKT*(MKTf-MKFf) + (VmK_Osta*MKTf/KmK_Osta+MKTf) - (VmK_Oat1*MKFf/KmK_Oat1+MKFf) - (VmK_Oat3*MKFf/KmK_Oat3+MKFf)
#Kidney proximal tubule cells subcompartment
dMKT = kKFKT*(MKFf-MKTf) + kFKT*(MFil-MKTf) + (VmK_Oatp*MFil/KmK_Oatp+MFil) - (VmK_Osta*MKTf/KmK_Osta+MKTf) + (VmK_Oat1*MKFf/KmK_Oat1+MKFf) + (VmK_Oat3*MKFf/KmK_Oat3+MKFf)
dMFil = kBF*(MBf-MFil) + kFKT*(MKTf-MFil) - (VmK_Oatp*MFil/KmK_Oatp+MFil) - (Qurine/VFil)*MFil
dMurine = (Qurine/VFil)*MFil

#Liver

#interstitial fluid subcompartment 
dMLF = kBLF*(MBf-MLFf) + kLFLT*(MLTf-MLFf) - (VmL_Oatp*MLFf/KmL_Oatp+MLFf) - (VmL_Ntcp*MLFf/KmL_Ntcp+MLFf)
#Liver tissue subcompartment
dMLT = kLFLT*(MLFf-MLTf) + (VmL_Oatp*MLFf/KmL_Oatp+MLFf) + (VmL_Ntcp*MLFf/KmL_Ntcp+MLFf) + kbileLT*(Mbile-MLTf)
dMbile = kbileLT*(MLTf-Mbile) - (Qbile/Vbile)*Mbile


#Gut

#interstitial fluid subcompartment 
dMGF = kBGF*(MBf-MGFf) + kGFGT*(MGT-MGFf) 
#Gut tissue subcompartment
dMGT = kGFGT*(MGFf-MGT) + kGLGT*(MGL-MGT)
dMGL = kGLGT*(MGT-MGL) + (Qbile/Vbile)*Mbile - (Qfeces/VGL)*MGL
dMfeces = (Qfeces/VGL)*MGL

#Muscle

#interstitial fluid subcompartment 
dMMF = kBMF*(MBf-MMFf) + kMFMT*(MMT-MMFf) 
#Muscle tissue subcompartment 
dMMT = kMFMT*(MMFf-MMT)

#Adipose

#interstitial fluid subcompartment 
dMAF = kBAF*(MBf-MAFf) + kAFAT*(MAT-MAFf) 
#Adipose tissue subcompartment 
dMAT = kAFAT*(MAFf-MAT)

#Rest of body

#interstitial fluid subcompartment 
dMRF = kBRF*(MBf-MRFf) + kRFRT*(MRT-MRFf) 
#Rest of body tissue subcompartment 
dMRT = kRFRT*(MRFf-MRT)



MBf = MB * 1.0 / (1.0 + CalbB * Ka)
MKFf = MKF * 1.0 / (1.0 + CalbKF * Ka)
MLFf = MLF * 1.0 / (1.0 + CalbLF * Ka)
MGFf = MGF * 1.0 / (1.0 + CalbGF * Ka)
MMFf = MMF * 1.0 / (1.0 + CalbMF * Ka)
MAFf = MAF * 1.0 / (1.0 + CalbAF * Ka)
MRFf = MRF * 1.0 / (1.0 + CalbRF * Ka)
MKTf = MKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp)
MLTf = MLT * 1.0 / (1.0 + CLfabpLT * KLfabp)


Cblood <- MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VB /Vplasma * 10^9
Ckidney <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VKB+MKFf+MKTf) / (VKB+VKT+VKF) * 10^9
Cliver <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VLB+MLFf+MLTf) / (VLB+VLT+VLF) * 10^9
Cgut <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VGB+MGFf+MGT) / (VGB+VGT+VGF) * 10^9
Cmuscle <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VMB+MMFf+MMT) / (VMB+VMT+VMF) * 10^9
Cadipose <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VAB+MAFf+MAT) / (VAB+VAT+VAF) * 10^9
Crest <- (MBf /(VB+VLB+VKB+VGB+VMB+VAB+VRB) * VRB+MRFf+MRT) / (VRB+VRT+VRF) * 10^9
Cfeces <- MGL / VGL * 10^9
Cbile <- Mbile / Vbile * 10^9
Curine <- MFil / VFil * 10^9

return(list(c('dMB'=dMB,'dMKF'=dMKF,'dMLF'=dMLF,'dMGF'=dMGF,'dMMF'=dMMF,'dMAF'=dMAF, 'dMRF'=dMRF,  
              'dMAT'=dMAT, 'dMMT'=dMMT, 'dMRT'=dMRT, 'dMKT'=dMKT, 'dMFil'=dMFil, 'dMurine'=dMurine,
              'dMLT'=dMLT, 'dMbile'=dMbile, 'dMGT'=dMGT, 'dMGL'=dMGL, 'dMfeces'=dMfeces,
              'MBf'=MBf, 'MKFf'=MKFf, 'MLFf'=MLFf, 'MGFf'=MGFf, 'MMFf'=MMFf, 'MAFf'=MAFf, 'MRFf'=MRFf, 'MKTf'=MKTf, 'MLTf'=MLTf),
            
            #'CalbB'=CalbB, 'CalbB'=CalbB, 'CalbLF'=CalbLF, 'CalbGF'=CalbGF, 'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF
            
            'Cblood'=Cblood, 'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cgut'=Cgut, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose, 'Crest'=Crest,
            'Cfeces'=Cfeces, 'Cbile'=Cbile, 'Curine'=Curine