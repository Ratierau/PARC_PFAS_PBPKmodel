#Epididymal fat
PVF <- 2.84/268 # https://doi.org/10.1039/c4fo00903g --> Figures 1,2
VF <- PVF * BW
PVFB <- NA
VFB <-PVFB * PVF * BW #volume of the blood of testis kg=L
PVFF <- NA
VFF <- PVFF * PVF * BW #testis interstitial fluid volume kg=L
VFT <- VF - VFF #testis tissue volume kg=L
PQBF <- 0.29/100 #https://doi.org/10.1152/ajpregu.1987.253.2.R228 Qcard=0.235*(0.335^0.75)*60 (L/h) and PQBF = (0.304*60/1000)/Qcard *100
QBF <- PQBF * Qcardiac #L/h


#Epididymal fat
epifat_cells = (5.326/1.517+5.425/3.082+6.935/5.489)*10^6/3 # cells/g tissue https://doi.org/10.1172/JCI105894 --> Tables I,II
Cmedium_F = 1*MW# Kimura et al. 2017 1uM umol/L -->  ug/L
epifat_protein <- muscle_protein#NA#5034
epifat_protein_total <- epifat_protein * (1000*VFT)
ClFFT <- ClFFT_unscaled * epifat_protein_total#uL/min
kFFFT <-  (60*ClFFT)/1e06 #L/h

#Epididymal fat

#blood subcompartment
dMBF = QBF*CBfart - (QBF-QBF/500)*CRBf - PeffF*AF*(CFBf-CFFf) - CFBf*QBF/500
#interstitial fluid subcompartment 
dMFF = CFBf*QBF/500 - CFFf*QBF/500 + PeffF*AF*(CFBf-CFFf) - kFFFT*(CFFf -CFT) 
#Epididymal fat tissue subcompartment 
dMFT = kFFFT*(CFFf -CFT)