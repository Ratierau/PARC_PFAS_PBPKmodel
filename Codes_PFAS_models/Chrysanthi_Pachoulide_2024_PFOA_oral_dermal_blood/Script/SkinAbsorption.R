
### As coded by Husoy ###
tchng = 50*365*24  # Duration of exposure (h); 50 years; turn dose on and off
Tinput = 24  # durration of dose (h)
tinterval = 24 # the interval the dose should be repeated (h). Do not have to be similar to Tinput

Skinarea = 972  # Exposed area on skin (cm^2); Changes from 5 since 2018 EFSA
SkinTarea = 9.1*(BW*1000)^0.666  # Total area of the skin (cm^2)
Skinthickness = 0.1 # Skin thickness (cm)
VSk = (Skinarea*Skinthickness)/1000  # Skin volume (L)
AbsPFOA =  0.016 #0.00048     # Changed to the absorption measured by Abraham and Monien 2022, of 1.6% of applied dose from sunscreen. 

PSk = 0.1  # Plasma/skin partition coefficient

QSkC = 0.058 # Fraction cardiac output going to skin
QCP = QC*(1-Htc)  # Cardiac output adjusted for plasma flow (L/h)

QSk <- QSkC*QCP*(Skinarea/SkinTarea) #ifelse(Dermconc>0.0,QSkC*QCP*(Skinarea/SkinTarea),0.0) # plasma flow to the skin
#      actualflowtoskin*fractionalareaexposed

Input2 = 0.0 # amount of PFOA entering the body via dermal exposure
Dermconc = 0 #as.numeric(SumExpPFOA_LB_val[i,14])
Dermdose = Dermconc*BW*AbsPFOA     # Internal dose from dermal absorption (Ug/day)

Input2 <- Dermdose/Tinput #*(t %% tinterval<Tinput)

## Skin compartment
RSk <- QSk*(CA*Free-CSk*FreeSk) + Input2*DoseOn # Rate of PFOA amount change in skin
#      fromtoplasma               changingtimeabsorption


### Westerhout ###
QSkC = 0.058  # Fraction cardiac output going to skin #same as Husoy

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


### My calculation ###

# Scaling in vitro to in vivo sking absorption

#To calculate skin area
 SA = 0.007184 * H^0.725 * W^0.425 #SA skin area in m2, H height in m, W weight in kg

Papp = 
Peff.human = (10^(0.4926*log10(ZEN.dat$Papp)- 0.1454));                           #human Peff(cm^-4/s)
Peff.rat = Peff.human/11.04;                                                      #rat Peff(cm^-4/s)
Peff = Peff.rat/10000*3600;                                                       #conversion to cm/h
ka = Peff*2/ZEN.dat$R;    