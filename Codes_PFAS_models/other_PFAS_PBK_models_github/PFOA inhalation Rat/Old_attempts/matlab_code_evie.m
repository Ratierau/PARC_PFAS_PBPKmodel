%% Diffusion-Limited PBPK for PFOA in Male Rats
%% Code written by Weixiao Cheng, University of Pittsburgh
%% Last Modified: 26 July, 2017
% The purpose of this model is to predict the toxicokinetics and tissue  
% distributioin of PFOA in male rats without the need to fit experimental   % data.
% This version implements Monte Carlo uncertainty analysis. Parameter
% values are therefore chosen from distributions (normal, lognormal, 
% or uniform, see SI Section S4 for details) over 10,000 iterations.
% The code below is parameterized to replicate the conditions of the 1 mg/kg 
% IV dose experiment from Kemper (33).
% For the code for other dose experiments, the paramters including body
% weight(BW), Dose, and sampling time (seconds) need to be modified to
% corresponding values based on different studies. Moreover, the initial
% condition is different between IV and oral dose treatment. Namely, for IV
% dose, the initial concentration of PFOA in the blood compartment is Dose
% times BW (i.e., in line 311 , y1 = Dose.*BW.*ones([1,Q])), with PFOA 
% concentration in other compartments being zero; for oral dose, the initial 
% concentration of PFOA in the gut lumen compartment is Dose time BW
% (i.e.,in line 340, y26 = Dose.*BW.*ones([1,Q])), with PFOA in other  24 % compartments being zero.
clc;
clear all;
close all;
Q = 1; % number of iterations for Monte Carlo analysis
BW = 0.244; % body weight (kg)
Dose = 1e-6; % the unit is kg PFOA/kg BW   
% Volume of tissue i as percentage of body weight (PVi, unitless) and 
% volume (Vi, m^3), assuming the density of tissue is 1e3 kg/m^3.

PVB = 54e-6;
VB = PVB.*BW; % blood volume
PVplasma = 31.2e-6;
Vplasma = PVplasma.*BW; % plasma volume
PVK = 0.73/100/1e3;
VK = PVK.*BW; % kidney volume
PVKB = 0.16;
VKB = PVKB.*PVK.*BW; % kidney blood volume
PVKF = 0.13;
VKF = PVKF.*PVK.*BW; % kidney interstitial fluid volume
VKT = VK - VKF; % kidney tissue volume
VFil = 0.25/1e6; % renal filtrate volume
PVL = 3.66/100/1e3;
VL = PVL.*BW; % liver volume
PVLB = 0.21;
VLB = PVLB.*PVL.*BW; % liver blood volume
PVLF = 0.049;
VLF = PVLF.*PVL.*BW; % liver interstitial fluid volume
VLT = VL - VLF; % liver tissue volume
PVbile = 0.004;
Vbile = PVbile.*PVL.*BW; % bile volume
PVG = 2.69/100/1e3;
VG = PVG.*BW; % gut volume
PVGB = 0.034;
VGB = PVGB.*PVG.*BW; % gut blood volume
PVGF = 0.28;
VGF = PVGF.*PVG.*BW; % gut interstitial fluid volume
VGT = VG - VGF; % gut tissue volume
PVGL = 4.5/100/1e3;
VGL = PVGL.*BW; % gut lumen volume
PVM = 40.43/100/1e3;
VM = PVM.*BW; % muscle volume
PVMB = 0.04;
VMB = PVMB.*PVM.*BW; % muscle blood volume
PVMF = 0.054;
VMF = PVMF.*PVM.*BW; % muscle interstitial fluid volume
VMT = VM - VMF; % muscle tissue volume
PVA = 7/100/1e3;
VA = PVA.*BW; % adipose volume
PVAB = 0.02;
VAB = PVAB.*PVA.*BW; % adipose blood volume
PVAF = 0.174;
VAF = PVAF.*PVA.*BW; % adipose interstitial fluid volume
VAT = VA - VAF; % adipose tissue volume
PVR = 1/1e3 - PVB - PVK - PVL - PVG - PVM - PVA;
VR = PVR.*BW; % volume of the rest of body
PVRB = 0.036; 
VRB = PVRB.*PVR.*BW; % volume of the blood of the rest of body
PVRF = 0.18; 
VRF = PVRF.*PVR.*BW; % interstitial fluid volume of the rest of body  
VRT = VR - VRF; % tissue volume of the rest of body

% Capillary surface area for each tissue (Ai) as percentage of body weight  
% or weight of corresponding tissue (PAi, unitless) and surface area (m^2). 
PAK = 350e-4; 
AK = PAK.*VK.*10^6; % kidney surface area
PAKG = 6890/1e6; 
AKG = PAKG.*VK.*10^6; % the surface area of glomerular capillary
PAL = 250e-4;
AL = PAL.*VL.*10^6; % liver surface area
PAG = 100e-4;
AG = PAG.*VG.*10^6; % gut surface area
PAGL = 4.14;
AGL = PAGL.*BW; % gut lumen surface area
PAM = 70e-4;
AM = PAM.*VM.*10^6; % muscle surface area
PAA = 70e-4;
AA = PAA.*VA.*10^6; % adipose surface area
PAR = 100e-4;
AR = PAR.*VR.*10^6; % surface area of rest of body

% Effective permeability (Peff, in m/s) for blood (B), liver(L), kidney(K), 109 % gut(G),adipose(A), muscle(M), rest of body(R).
PeffB = 4.98e-8;
PeffK = 4.38e-8;
PeffL = 5.15e-8;
PeffG = 2.65e-8;
PeffA = 2.65e-8;
PeffM = 2.65e-8;
PeffR = 2.65e-8;

% Steady-state cell-water concentration ratios(CRss) for gut, liver, and 
% kidney.
CRssG = 3.75;
CRssL = 7.28;
CRssK = 6.19;

% Blood flow rates (QBi, in m^3/s) to different tissues (i=L, K, G, A, M, R) 126 % as a percentage of cardiac output (Qcardiac), which itself is a function 127 % of body weight (BW).
Qcardiac = 0.235/60*1e-3*BW.^0.75;
PQBK = 14.1/100;
QBK = PQBK.*Qcardiac;
PQBG = 15.1/100;
QBG = PQBG.*Qcardiac; 
PQBL = 2.4/100;
QBL = (PQBL+PQBG).*Qcardiac; 
PQBM = 27.8/100;
QBM = PQBM.*Qcardiac; 
PQBA = 7/100;
QBA = PQBA.*Qcardiac; 
PQBR = 1 - PQBK - PQBG - PQBL - PQBM - PQBA;
QBR = PQBR.*Qcardiac;

% Flow rate of fluids including feces, bile, urine and glomerular filtration 
% rate (GFR), in m^3/s.
Qfeces = 5.63*1e-6/(24*3600);
PQbile = 90;
Qbile = PQbile.*BW.*1e-6/(24*3600);
PQurine = 200;
Qurine = PQurine.*BW.*1e-6/(24*3600);
PQGFR = 10.74;
QGFR = PQGFR.*BW.*1e-6/60;

% Albumin concentration in blood and interstitial fluid compartments(mol/m^3).
CalbB = 281e-3*7.8;
CalbKF = 243e-3*7.8;
CalbLF = 243e-3*7.8;
CalbGF = 146e-3*7.8;
CalbMF = 146e-3*7.8;
CalbAF = 73e-3*7.8;
CalbRF = 73e-3*7.8;

% Alpha2mu-globulin concentration in kidney tissue (mol/m^3). 
Ca2uKT = 110e-3;
% LFABP concentration in kidney and liver tissue (mol/m^3).
CL_fabpKT = 2.65e-3*3;
CL_fabpKT1 = CL_fabpKT./3;
CL_fabpKT2 = CL_fabpKT./3;
CL_fabpKT3 = CL_fabpKT./3;
CL_fabpLT = 133e-3*3;
CL_fabpLT1 = CL_fabpLT./3;
CL_fabpLT2 = CL_fabpLT./3;
CL_fabpLT3 = CL_fabpLT./3;

% Equilibrium association constant (m^3/mol) for albumin(Ka), LFABP(KL_fabp), 
% and alpha2mu-globulin(Ka2u). See SI section S2-2 for details.
Ka = 3.1;
KL_fabp1 = 120; 
KL_fabp2 = 40.0; 
KL_fabp3 = 19.0; 
Ka2u = 0.5;

% Individual rate constants for association and dissociation(s^-1 and m^3/mol*s). 
% Note kon/koff=Keq.
koff = 0.01; % assume koff is 0.01/s
kon = koff.*Ka;
kL_fabpon1 = koff.*KL_fabp1;
kL_fabpon2 = koff.*KL_fabp2;
kL_fabpon3 = koff.*KL_fabp3;
kK_fabpon = koff.*Ka2u;

% Overall mass transfer coefficients between subcompartments and passive 
% diffusion rate constants. See SI section S3-1 for details.
kBKF = ((1./QBK) + 1./(PeffB.*AK)).^(-1);
kBF =  PeffB.*AKG;
kKFKT = PeffK.*AK;
n = 5; % enlargement factor of apical membrane of proximal tubule 
kFKT = PeffK.*AK.*n;
kBLF = ((1./QBL) + 1./(PeffB.*AL)).^(-1);
kLFLT = PeffL.*AL;
kBGF = ((1./QBG) + 1./(PeffB.*AG)).^(-1);
kGFGT = PeffG.*AG;
kGLGT = PeffG.*AGL;
kBMF = ((1./QBM) + 1./(PeffB.*AM)).^(-1);
kMFMT = PeffM.*AM;
kBAF = ((1./QBA) + 1./(PeffB.*AA)).^(-1);
kAFAT = PeffA.*AA;
kBRF = ((1./QBR) + 1./(PeffB.*AR)).^(-1);
kRFRT = PeffR.*AR;
kbileLT = PeffL.*AL;

% First-order rate constants (s^-1).
bBKF = kBKF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bKFB = kBKF./VKF;
bKFKT = kKFKT./VKF;
bKTKF = kKFKT./VKT;
bFKT = kFKT./VFil; 
bKTF = kFKT./(VKT.*CRssK);
bBF = QGFR./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bFB = kBF./VFil;
bBLF = kBLF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bLFB = kBLF./VLF;
bLFLT = kLFLT./VLF;
bLTLF = kLFLT./VLT;
bbileLT = kbileLT./Vbile;
bLTbile = kbileLT./(VLT.*CRssL);
bBGF = kBGF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bGFB = kBGF./VGF;
bGFGT = kGFGT./VGF;
bGTGF = kGFGT./VGT;
bGLGT = kGLGT./VGL;
bGTGL = kGLGT./(VGT.*CRssG);
bBMF = kBMF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bMFB = kBMF./VMF;
bMFMT = kMFMT./VMF;
bMTMF = kMFMT./VMT;
bBAF = kBAF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bAFB = kBAF./VAF;
bAFAT = kAFAT./VAF;
bATAF = kAFAT./VAT;
bBRF = kBRF./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
bRFB = kBRF./VRF;
bRFRT = kRFRT./VRF;
bRTRF = kRFRT./VRT;

% First-order rate constants (s^-1) for protein-mediated transport, see 
% section S3-3 for details.
Pbclear = 2.76e-7;
bclear = Pbclear.*AK./VKF; 
Pbreab = 1.18e-7;
breab = n.*Pbreab.*AK./VFil; 
Pbabs = 1.78e-7;
babs = Pbabs.*AL./VLF;
Pbefflux = 1.38e-7; 
befflux = Pbefflux.*AK./VKT;

% Conversion between mass and concentration for protein content of tissues.
MalbB = CalbB.*(VB+VLB+VKB+VGB+VMB+VAB+VRB);
MalbKF = CalbKF.*VKF;
ML_fabpKT1 = CL_fabpKT1.*VKT;
ML_fabpKT2 = CL_fabpKT2.*VKT;
ML_fabpKT3 = CL_fabpKT3.*VKT;
MK_fabpKT = Ca2uKT.*VKT;
MalbLF = CalbLF.*VLF;
ML_fabpLT1 = CL_fabpLT1.*VLT;
ML_fabpLT2 = CL_fabpLT2.*VLT;
ML_fabpLT3 = CL_fabpLT3.*VLT;
MalbGF = CalbGF.*VGF;
MalbMF = CalbMF.*VMF;
MalbAF = CalbAF.*VAF;
MalbRF = CalbRF.*VRF;

% Below is the numerical method used to solve mass balance equations.
seconds = 22*24*3600; % simulation time, 22 days
h = 0.5; % step size
tspan = (1:h:seconds);
steps = seconds./h;

% Initial condition for each comparment.
% Mass of PFOA in blood not bound to proteins: y1
% Mass of PFOA in blood bound to albumin: y2
% Mass of PFOA in interstitial fluid of kidney not bound to proteins: y3
% Mass of PFOA in interstitial fluid of kidney bound to albumin: y4
% Mass of PFOA in kidney tissue not bound to proteins: y5
% Mass of PFOA in kidney tissue bound to LFABP: y6, y61, and y62 (LFABP 
% has 3 binding sites)
% Mass of PFOA in kidney tissue bound to alpha2mu-globulin: y7
% Mass of PFOA in renal filtrate not bound to proteins: y8
% Mass of PFOA in interstitial fluid of liver not bound to proteins: y9
% Mass of PFOA in interstitial fluid of liver bound to albumin: y10
% Mass of PFOA in liver tissue not bound to proteins: y11
% Mass of PFOA in liver tissue bound to LFABP: y12, y121, and y122 (LFABP 
% has 3 binding sites)
% Mass of PFOA in bile not bound to proteins: y13
% Mass of PFOA in interstitial fluid of gut not bound to proteins: y14
% Mass of PFOA in interstitial fluid of gut bound to albumin: y15
% Mass of PFOA in gut tissue not bound to proteins: y16
% Mass of PFOA in interstitial fluid of muscle not bound to proteins: y17
% Mass of PFOA in interstitial fluid of muscle bound to albumin: y18
% Mass of PFOA in muscle tissue not bound to proteins: y19
% Mass of PFOA in interstitial fluid of adipose not bound to proteins: y20
% Mass of PFOA in interstitial fluid of adipose bound to albumin: y21
% Mass of PFOA in adipose tissue not bound to proteins: y22
% Mass of PFOA in interstitial fluid of the rest of body not bound to proteins: y23
% Mass of PFOA in interstitial fluid of the rest of body bound to albumin: y24
% Mass of PFOA in tissue of the rest of body not bound to proteins: y25
% Mass of PFOA in gut lumen not bound to proteins: y26
y1 = Dose.*BW.*ones([1,Q]); 
y2 = zeros([1,Q]);
y3 = zeros([1,Q]);
y4 = zeros([1,Q]);
y5 = zeros([1,Q]);
y6 = zeros([1,Q]);
y61 = zeros([1,Q]);
y62 = zeros([1,Q]);
y7 = zeros([1,Q]);
y8 = zeros([1,Q]);
y9 = zeros([1,Q]);
y10 = zeros([1,Q]);
y11 = zeros([1,Q]);
y12 = zeros([1,Q]);
y121 = zeros([1,Q]);
y122 = zeros([1,Q]);
y13 = zeros([1,Q]);
y14 = zeros([1,Q]);
y15 = zeros([1,Q]);
y16 = zeros([1,Q]);
y17 = zeros([1,Q]);
y18 = zeros([1,Q]);
y19 = zeros([1,Q]);
y20 = zeros([1,Q]);
y21 = zeros([1,Q]);
y22 = zeros([1,Q]);
y23 = zeros([1,Q]);
y24 = zeros([1,Q]);
y25 = zeros([1,Q]);
y26 = zeros([1,Q]);
for j = 1:steps

% Mass balance for available protein binding sites.
MalbB_new = MalbB+h.*(koff.*y2-kon.*MalbB.*y1./(VB+VLB+VKB+VGB+VMB+VAB+VRB));
MalbB = MalbB_new;
CalbB = MalbB./(VB+VLB+VKB+VGB+VMB+VAB+VRB);
MalbKF_new = MalbKF+h.*(koff.*y4-kon.*MalbKF.*y3./VKF);
MalbKF = MalbKF_new;
CalbKF = MalbKF./VKF;
ML_fabpKT1_new = ML_fabpKT1+h.*(koff.*y6-kL_fabpon1.*ML_fabpKT1.*y5./VKT);
ML_fabpKT1 = ML_fabpKT1_new;
CL_fabpKT1 = ML_fabpKT1./VKT;
ML_fabpKT2_new = ML_fabpKT2+h.*(koff.*y61-kL_fabpon2.*ML_fabpKT2.*y5./VKT);
ML_fabpKT2 = ML_fabpKT2_new;
CL_fabpKT2 = ML_fabpKT2./VKT;
ML_fabpKT3_new = ML_fabpKT3+h.*(koff.*y62-kL_fabpon3.*ML_fabpKT3.*y5./VKT);
ML_fabpKT3 = ML_fabpKT3_new;
CL_fabpKT3 = ML_fabpKT3./VKT;
MK_fabpKT_new = MK_fabpKT+h.*(koff.*y7-kK_fabpon.*MK_fabpKT.*y5./VKT);
MK_fabpKT = MK_fabpKT_new;
Ca2uKT = MK_fabpKT./VKT;
MalbLF_new = MalbLF+h.*(koff.*y10-kon.*MalbLF.*y9./VLF);
MalbLF = MalbLF_new;
CalbLF = MalbLF./VLF;
ML_fabpLT1_new = ML_fabpLT1 + h.*(koff.*y12 - kL_fabpon1.*ML_fabpLT1.*y11./VLT);
ML_fabpLT1 = ML_fabpLT1_new;
CL_fabpLT1 = ML_fabpLT1./VLT;
ML_fabpLT2_new = ML_fabpLT2+h.*(koff.*y121-kL_fabpon2.*ML_fabpLT2.*y11./VLT);
ML_fabpLT2 = ML_fabpLT2_new;
CL_fabpLT2 = ML_fabpLT2./VLT;
ML_fabpLT3_new = ML_fabpLT3+h.*(koff.*y122-kL_fabpon3.*ML_fabpLT3.*y11./VLT);
ML_fabpLT3 = ML_fabpLT3_new;
CL_fabpLT3 = ML_fabpLT3./VLT;
MalbGF_new = MalbGF+h.*(koff.*y15-kon.*MalbGF.*y14./VGF);
MalbGF = MalbGF_new;
CalbGF = MalbGF./VGF;
MalbMF_new = MalbMF+h.*(koff.*y18-kon.*MalbMF.*y17./VMF);
MalbMF = MalbMF_new;
CalbMF = MalbMF./VMF;
MalbAF_new = MalbAF+h.*(koff.*y21-kon.*MalbAF.*y20./VAF);
MalbAF = MalbAF_new;
CalbAF = MalbAF./VAF;
MalbRF_new = MalbRF+h.*(koff.*y24-kon.*MalbRF.*y23./VRF);
MalbRF = MalbRF_new;
CalbRF = MalbRF./VRF;
bBon = CalbB.*kon;
bBoff = koff;
bKFon = CalbKF.*kon;
bKFoff = koff;
bL_fabpKTon1 = CL_fabpKT1.*kL_fabpon1;
bL_fabpKTon2 = CL_fabpKT2.*kL_fabpon2;
bL_fabpKTon3 = CL_fabpKT3.*kL_fabpon3;
bL_fabpKToff = koff;
bK_fabpKTon = Ca2uKT.*kK_fabpon;
bK_fabpKToff = koff;
bLFon = CalbLF.*kon;
bLFoff = koff;
bL_fabpLTon1 = CL_fabpLT1.*kL_fabpon1;
bL_fabpLTon2 = CL_fabpLT2.*kL_fabpon2;
bL_fabpLTon3 = CL_fabpLT3.*kL_fabpon3;
bL_fabpLToff = koff;
bGFon = CalbGF.*kon;
bGFoff = koff;
bMFon = CalbMF.*kon;
bMFoff = koff;
bAFon = CalbAF.*kon;
bAFoff = koff;
bRFon = CalbRF.*kon;
bRFoff = koff;

% Differential mass balance for each tissue or fluid compartment.
y1_new = y1 + h.*(bKFB.*y3+bLFB.*y9-bBLF.*y1-bBKF.*y1-bBon.*y1+bBoff.*y2-bBF.*y1+bFB.*y8+bGFB.*y14-bBGF.*y1+bMFB.*y17-bBMF.*y1+bAFB.*y20-bBAF.*y1+bRFB.*y23-bBRF.*y1);
y1 = max(y1_new,0);
y2_new = y2 + h.*(bBon.*y1-bBoff.*y2);
y2 = max(y2_new,0);
y3_new = y3 + h.*(bBKF.*y1-bKFB.*y3+befflux.*y5+bKTKF.*y5-bKFKT.*y3-bclear.*y3+bKFoff.*y4-bKFon.*y3);
y3 = max(y3_new,0);
y4_new = y4 + h.*(bKFon.*y3-bKFoff.*y4);
y4 = max(y4_new,0);
y5_new = y5 + h.*(bKFKT.*y3+bFKT.*y8+breab.*y8+bclear.*y3-befflux.*y5-bKTKF.*y5-bKTF.*y5+bL_fabpKToff.*y6-bL_fabpKTon1.*y5+bL_fabpKToff.*y61-bL_fabpKTon2.*y5+bL_fabpKToff.*y62-bL_fabpKTon3.*y5+bK_fabpKToff.*y7-bK_fabpKTon.*y5);
y5 = max(y5_new,0);
y6_new = y6 + h.*(bL_fabpKTon1.*y5-bL_fabpKToff.*y6);
y6 = max(y6_new,0);
y61_new = y61 + h.*(bL_fabpKTon2.*y5-bL_fabpKToff.*y61);
y61 = max(y61_new,0);
y62_new = y62 + h.*(bL_fabpKTon3.*y5-bL_fabpKToff.*y62);
y62 = max(y62_new,0);
y7_new = y7 + h.*(bK_fabpKTon.*y5-bK_fabpKToff.*y7);
y7 = max(y7_new,0);
y8_new = y8 + h.*(bKTF.*y5-breab.*y8+bBF.*y1-bFB.*y8-bFKT.*y8-Qurine.*y8./VFil);
y8 = max(y8_new,0);
y9_new = y9 + h.*(bBLF.*y1-bLFB.*y9+bLTLF.*y11-bLFLT.*y9-babs.*y9-bLFon.*y9+bLFoff.*y10);
y9 = max(y9_new,0);
y10_new = y10 + h.*(bLFon.*y9-bLFoff.*y10);
y10 = max(y10_new,0);
y11_new = y11 + h.*(bLFLT.*y9+babs.*y9+bbileLT.*y13-bLTbile.*y11-bLTLF.*y11-bL_fabpLTon1.*y11+bL_fabpLToff.*y12-bL_fabpLTon2.*y11+bL_fabpLToff.*y121-bL_fabpLTon3.*y11+bL_fabpLToff.*y122);
y11 = max(y11_new,0);
y12_new = y12 + h.*(bL_fabpLTon1.*y11-bL_fabpLToff.*y12);
y12 = max(y12_new,0);
y121_new = y121 + h.*(bL_fabpLTon2.*y11-bL_fabpLToff.*y121);
y121 = max(y121_new,0);
y122_new = y122 + h.*(bL_fabpLTon3.*y11-bL_fabpLToff.*y122);
y122 = max(y122_new,0);
y13_new = y13 + h.*(bLTbile.*y11-bbileLT.*y13-Qbile.*y13./Vbile);
y13 = max(y13_new,0);
y14_new = y14 + h.*(bBGF.*y1-bGFB.*y14+bGTGF.*y16-bGFGT.*y14-bGFon.*y14+bGFoff.*y15);
y14 = max(y14_new,0);
y15_new = y15 + h.*(bGFon.*y14-bGFoff.*y15);
y15 = max(y15_new,0);
y16_new = y16 + h.*(bGFGT.*y14-bGTGF.*y16+bGLGT.*y26-bGTGL.*y16);
y16 = max(y16_new,0);
y17_new = y17 + h.*(bBMF.*y1-bMFB.*y17+bMTMF.*y19-bMFMT.*y17-bMFon.*y17+bMFoff.*y18);
y17 = max(y17_new,0);
y18_new = y18 + h.*(bMFon.*y17-bMFoff.*y18);
y18 = max(y18_new,0);
y19_new = y19 + h.*(bMFMT.*y17-bMTMF.*y19);
y19 = max(y19_new,0);
y20_new = y20 + h.*(bBAF.*y1-bAFB.*y20+bATAF.*y22-bAFAT.*y20-bAFon.*y20+bAFoff.*y21);
y20 = max(y20_new,0);
y21_new = y21 + h.*(bAFon.*y20-bAFoff.*y21);
y21 = max(y21_new,0);
y22_new = y22 + h.*(bAFAT.*y20-bATAF.*y22);
y22 = max(y22_new,0);
y23_new = y23 + h.*(bBRF.*y1-bRFB.*y23+bRTRF.*y25-bRFRT.*y23-bRFon.*y23+bRFoff.*y24);
y23 = max(y23_new,0);
y24_new = y24 + h.*(bRFon.*y23-bRFoff.*y24);
y24 = max(y24_new,0);
y25_new = y25 + h.*(bRFRT.*y23-bRTRF.*y25);
y25 = max(y25_new,0);   
y26_new = y26 + h.*(bGTGL.*y16-bGLGT.*y26+Qbile.*y13./Vbile-Qfeces.*y26./VGL);
y26 = max(y26_new,0);

yBfree = y1;
yBbound = y2;
yKFfree = y3;
yKFbound = y4;
yKTfree = y5;
yKTLbound1 = y6;
yKTLbound2 = y61;
yKTLbound3 = y62;
yKTKbound = y7;
yFfree = y8;
yLFfree = y9;
yLFbound = y10;
yLTfree = y11;
yLTbound1 = y12;
yLTbound2 = y121;
yLTbound3 = y122;
yBile = y13;
yGFfree = y14;
yGFbound = y15;
yGTfree = y16;
yMFfree = y17;
yMFbound = y18;
yMTfree = y19;
yAFfree = y20;
yAFbound = y21;
yATfree = y22;
yRFfree = y23;
yRFbound = y24;
yRTfree = y25;
yGLfree = y26;

i = round(1 + j./3600)
t(i,:) = j*h./(24*3600);
% Unit conversion from kg/m^3 to ng/g.

Blood_f(i,:) = yBfree;
Blood_b(i,:) = yBbound;
Kidney_Tf(i,:)= yKTfree;
Kidney_Tb1(i,:)= yKTLbound1;
Kidney_Tb2(i,:)=yKTLbound2;
Kidney_Tb3(i,:)=yKTLbound3;
Kidney_Tb(i,:)=yKTKbound;
Liver_Tf(i,:)= yLTfree;
Liver_Tb1(i,:)= yLTbound1;
Liver_Tb2(i,:)=yLTbound2;
Liver_Tb3(i,:)=yLTbound3;
Gut_Tf(i,:)= yGTfree;
Muscle_Tf(i,:)= yMTfree;
Adipose_Tf(i,:)= yATfree;
Rest_Tf(i,:) = yRTfree;
Bile(i,:) = yBile;
Kidney_Ff(i,:)= yKFfree;
Kidney_Fb(i,:)= yKFbound;
Liver_Ff(i,:)= yLFfree;
Liver_Fb(i,:)= yLFbound;
Gut_Ff(i,:)= yGFfree;
Gut_Fb(i,:)= yGFbound;
Gut_Lf(i,:)= yGLfree;
Muscle_Ff(i,:)= yMFfree;
Muscle_Fb(i,:)= yMFbound;
Adipose_Ff(i,:)= yAFfree;
Adipose_Fb(i,:)= yAFbound;
Rest_Ff(i,:)= yRFfree;
Rest_Fb(i,:)= yRFbound;
Ff(i,:)= yFfree;

C_Blood(i,:) = (yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VB./Vplasma.*10^6; 
C_Kidney(i,:)=((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VKB+yKFfree-yKFbound+yKTfree+yKTLbound1+yKTLbound2+yKTLbound3+yKTKbound)./(VKB+VKT+VKF).*10^6;
C_Liver(i,:)=((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VLB+yLFfree+yLFbound+yLTfree+yLTbound1+yLTbound2+yLTbound3)./(VLB+VLT+VLF).*10^6; 
C_Gut(i,:)=((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VGB+yGFfree+yGFbound+ yGTfree)./(VGB+VGT+VGF).*10^6;
C_Muscle(i,:)=((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VMB+yMFfree+yMFbound+yMTfree)./(VMB+VMT+VMF).*10^6; 
C_Adipose(i,:) = ((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VAB+yAFfree+yAFbound+yATfree)./(VAB+VAT+VAF).*10^6; 
C_Rest(i,:) = ((yBfree+yBbound)./(VB+VLB+VKB+VGB+VMB+VAB+VRB).*VRB+yRFfree+yRFbound+yRTfree)./(VRB+VRT+VRF).*10^6;
C_Feces(i,:) = yGLfree./VGL.*10^6;
C_Bile(i,:) = yBile./Vbile.*10^6;
C_Urine(i,:) = yFfree./VFil.*10^6;
   
end

% Sensitivity analysis using Pearson ranked correlation between sampled
% parameters and predicted PFOA concentrations, shown here for blood, 
% sampled on Day 12.
rho1 = corr((Blood(4115,:))',[BW',Dose',PVB',PVplasma',PVK',PVKB',...
PVKF',VFil',PVL',PVLB',PVLF',PVbile',PVG',PVGB',PVGF',PVGL',...
PVM',PVMB',PVMF',PVA',PVAB',PVAF',PVRB',PVRF',PQBK',...
PQGFR',PQurine',PQBL',PQbile',Qfeces',PQBG',PQBM',PQBA',...
PAK',n',PAKG',PAL',PAG',PAGL',PAM',PAA',PAR',PeffB',PeffK',...
PeffL',PeffG',PeffM',PeffA',PeffR',CRssG',CRssL',CRssK',...
CalbB',CalbKF',CL_fabpKT',Ca2uKT',CalbLF',CL_fabpLT',CalbGF',...
CalbMF',CalbAF',CalbRF',Ka',KL_fabp1',KL_fabp2',KL_fabp3',Ka2u',...
koff',Pbclear',Pbreab',Pbefflux',Pbabs']); 
