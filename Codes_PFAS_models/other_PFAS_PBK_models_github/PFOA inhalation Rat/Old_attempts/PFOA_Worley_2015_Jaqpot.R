# The code replicates the Loccisano et al.(2012) PBK model
# and simulates IV or Oral PFOA distribution for female or male rats of constant or varying
# body weight. The Bernstein et al. (2021) code was used for replicating the
# model and the parameter set was expanded to also describe female rats
# (the original code only described male rats). 
# We did not include simultaneous IV and oral because there are no such experiments

library(deSolve)

#=========================
#1. Parameters of the model
#=========================

create.params  <- function(user_input){
  with( as.list(user_input),{
    # User input: BW(kg), sex(F/M)
    #Physiological Parameters
    #BW = 0.3	#bodyweight (kg)
    MKC = 0.0084	#fraction mass of kidney (percent of BW); Brown 1997
    MLC = 0.034	#fraction mass of liver (percent of BW); Brown 1997
    
    #Cardiac Output and Bloodflow (as fraction of cardiac output)
    QCC = 14.0 #cardiac output in L/h/kg^0.75; Brown 1997		 
    QLC = 0.183 #fraction blood flow to liver; Brown 1997	
    QKC = 0.141 #fraction blood flow to kidney; Brown 1997.
    Htc = 0.46 #hematocrit for the rat; Davies 1993
    
    #Tissue Volumes 
    VplasC = 0.0312 #fraction vol. of plasma (L/kg BW); Davies 1993
    VLC = 0.035 #fraction vol. of liver (L/kg BW); Brown 1997
    VKC = 0.0084 #fraction vol. of kidney (L/kg BW); Brown 1997
    VfilC = 8.4e-4	#fraction vol. of filtrate (L/kg BW)
    VPTCC = 1.35e-4 #vol. of proximal tubule cells (L/g kidney) (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
    
    #Chemical Specific Parameters
    MW = 414.07	#PFOA molecular mass (g/mol)
   
    #Kidney Transport Parameters
    Vmax_baso_invitro = 393.45 #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of rOAT1 and rOAT3 from Nakagawa, 2007
    Km_baso = 27.2 #Km of basolateral transporter (mg/L) Average of rOAT1 and rOAT3 from Nakagawa et. al, 2007
    Vmax_apical_invitro = 9300 #Vmax of apical transporter (pmol/mg protein/min); invitro value for Oatp1a1 from Weaver, 2010
    Km_apical = 52.3#Km of apical transporter (mg/L), in vitro value for Oatp1a1 from Weaver, 2010.
    protein = 2.0e-6	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
     
    #Partition Coefficients (from rat tissue data, Kudo et al, 2007)
    PL = 2.2 #liver:blood
    PK = 1.05 #kidney:blood
    PR = 0.11 #rest of body:blood
    
    #rate constants
    kdif = 0.001	#diffusion rate from proximal tubule cells (L/h)
    kabsc = 2.12	#rate of absorption of chemical from small intestine to liver (1/(h*BW^0.25))(fit to data)
    kunabsc = 7.06e-5	#rate of unabsorbed dose to appear in feces (1/(h*BW^0.25))(fit to data)
    GEC = 0.54 #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013
    k0C = 1.0	#rate of uptake from the stomach into the liver (1/(h*BW^0.25)) (fit to data)
    
    keffluxc = 2.49 #rate of clearance of PFOA from proximal tubule cells into blood (1/(h*BW^0.25))
    kbilec = 0.004 #biliary elimination rate ((male); liver to feces storage (1/(h*BW^0.25)) (fit to data)
    kurinec = 1.6 #rate of urine elimination from urine storage (male) (1/(h*BW^0.25))(fit to data)
                                                       
    if(sex == "M"){
      Free = 0.006 #free fraction in plasma (male)
      GFRC = 62.1	#glomerular filtration rate (L/hr/kg kidney) (male); Corley, 2005
      RAFbaso = 3.99	#relative activity factor, basolateral transporters (male) (fit to data)	
      RAFapi = 4.07	#relative activity factor, apical transporters (male) (fit to data)	
    }else if(sex == "F"){
      Free = 0.09 #free fraction in plasma (female)
      GFRC = 41.04	#glomerular filtration rate (L/hr/kg kidney) (female); Corley, 2005
      RAFapi = 0.001356	#relative activity factor, apical transporters (female) (fit to data)
      RAFbaso = 0.01356	#relative activity factor, basolateral transporters (female) (fit to data)
    }
    
    
    #Scaled Parameters
    #Cardiac output and blood flows
    QC = QCC*(BW^0.75)*(1-Htc)	#cardiac output in L/h; adjusted for plasma
    QK = (QKC*QC)	#plasma flow to kidney (L/h)
    QL = (QLC*QC)	#plasma flow to liver (L/h)
    QR = QC - QK - QL 	#plasma flow to rest of body (L/h)
    QBal = QC - (QK + QL + QR) #Balance check of blood flows; should equal zero
    
    #Tissue Volumes
    MK = MKC*BW*1000	#mass of the kidney (g)
    PTC = MKC*6e7	#number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
    VPTC = MK*VPTCC	#volume of proximal tubule cells (L)	
    VPlas = VplasC*BW 	#volume of plasma (L) 
    VK = VKC*BW 	#volume of kidney (L)
    VKb = VK*0.16	#volume of blood in the kidney (L); fraction blood volume of kidney (0.16) from Brown, 1997
    Vfil = VfilC*BW	#volume of filtrate (L)
    VL = VLC*BW	#volume of liver (L)
    VR = (0.93*BW) - VPlas - VPTC - Vfil - VL	#volume of remaining tissue (L); 
    VBal = (0.93*BW) - (VR + VL + VPTC + Vfil + VPlas)	#Balance check of tissue volumes; should equal zero 
    ML = MLC*BW*1000	#mass of the liver (g)
    
    #Kidney Parameters
    MK = MKC*BW*1000	#mass of the kidney (g)
    PTC = MKC*6e7	#number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
    VPTC = MK*VPTCC	#volume of proximal tubule cells (L)	
    MPTC = VPTC*1000 #mass of the proximal tubule cells (g) (assuming density 1 kg/L)	
    Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000)#Vmax of basolateral transporters (mg/h/kg BW)
    Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000) #Vmax of basolateral transporters (mg/h/kg BW)
    Vmax_baso = Vmax_basoC*BW^0.75	#(mg/h)
    Vmax_apical = Vmax_apicalC*BW^0.75	#(mg/h)
    kbile = kbilec*BW^(-0.25)	#biliary elimination; liver to feces storage (/h)
    kurine = kurinec*BW^(-0.25)	#urinary elimination, from filtrate (/h)
    kefflux = keffluxc*BW^(-0.25)	#efflux clearance rate, from PTC to blood (/h)
    GFR = GFRC*(MK/1000)	#glomerular filtration rate, scaled to mass of kidney(in kg)(L/h)
    
    #GI Tract Parameters
    kabs = kabsc*BW^(-0.25)	#rate of absorption of chemical from small intestine to liver (/h)
    kunabs = kunabsc*BW^(-0.25)	#rate of unabsorbed dose to appear in feces (/h)
    GE = GEC/BW^0.25	#gastric emptying time (/h)
    k0 = k0C/BW^0.25 	#rate of uptake from the stomach into the liver (/h)
    
    return(list( "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR, 
                 "VPlas" = VPlas,
                 "VKb" = VKb, "Vfil" = Vfil, "VL" = VL, "VR" = VR, "ML" = ML,
                 "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
                 'kdif' = kdif,
                 "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
                 "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, "GE" = GE, "k0" = k0))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    "AR" = 0; "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "AST" = 0;
    "AabsST" = 0; "ASI" = 0; "AabsSI" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    
    return(c("AR" = AR, "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
             "ACl" = ACl, "Aefflux" = Aefflux,
             "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
             "Aurine" = Aurine, "AST" =AST,
             "AabsST" = AabsST, "ASI" = ASI, "AabsSI" = AabsSI, "Afeces" = Afeces, 
             "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free))
  })
}

#===================
#3. Events function
#===================

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
        events <- list(data = rbind(data.frame(var = c("Aplas_free"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("AST"),  time = admin.time, 
                                               value = c( admin.dose), method = c("add")) ))
      }
    }
    return(events)
  })
}

#==================
#4. Custom function 
#==================
custom.func <- function(){
  return()
}

#==============
#5. ODEs System
#==============

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    
    CR = AR/VR #concentration in rest of body (mg/L)
    CVR = CR/PR	#concentration in venous blood leaving the rest of the body (mg/L)
    CKb = AKb/VKb	#concentration in kidney blodd (mg/L) 
    CVK = CKb #/PK	#concentration in venous blood leaving kidney (mg/L)
    CPTC = APTC/VPTC	#concentration in PTC (mg/L)
    Cfil = Afil/Vfil	#concentration in filtrate (mg/L)
    CL = AL/VL	#concentration in the liver (mg/L)
    CVL = CL/PL	#concentration in the venous blood leaving the liver (mg/L)
    CA_free = Aplas_free/VPlas		#concentration in plasma (mg)
    CA = CA_free/Free	#concentration of total PFOA in plasma (mg/L)
    
    # Rest of Body (Tis)
    dAR = QR*(CA-CVR)*Free	#rate of change in rest of body (mg/h)
   
    #Kidney 
    #Kidney Blood (Kb)
    dAdif = kdif*(CKb - CPTC)	#rate of diffusion from into the PTC (mg/hr)
    dA_baso = (Vmax_baso*CKb)/(Km_baso + CKb)	
    dAKb = QK*(CA-CVK)*Free - CA*GFR*Free - dAdif - dA_baso #rate of change in kidney blood (mg/h).
    dACl = CA*GFR*Free	#rate of clearance via glormerular filtration (mg/h)
   
    #Proximal Tubule Cells (PTC)
    dAefflux = kefflux*APTC
    dA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil)
    dAPTC =  dAdif + dA_apical + dA_baso - dAefflux #rate of change in PTC(mg/h)
    
    #Filtrate (Fil)
    dAfil = CA*GFR*Free - dA_apical - Afil*kurine	#rate of change in filtrate (mg/h)
    
    #Urinary elimination 
    dAurine = kurine*Afil	#rate of change in urine (mg/h)

    #GI Tract (Absorption site of oral dose)
    #Stomach
    dAST=  - k0*AST - GE*AST	#rate of change in the stomach (mg/h)
    dAabsST = k0*AST	#rate of absorption in the stomach (mg/h)

    #Small Intestine
    dASI = GE*AST - kabs*ASI - kunabs*ASI	#rate of change in the small intestine (mg/hr)
    dAabsSI = kabs*ASI	#rate of absorption inthe small intestine (mg/hr)

    total_oral_uptake = AabsSI + AabsST	#total oral uptake in the GI tract (mg) 
    
    #Feces compartment
    dAfeces = kbile*AL + kunabs*ASI #rate of change in the feces compartment (mg/h)

    #Liver
    dAL = QL*(CA-CVL)*Free - kbile*AL + kabs*ASI + k0*AST #rate of change in the liver (mg/h)
    dAbile = kbile*AL
    amount_per_gram_liver = (AL/ML)*1000	#amount of PFOA in liver per gram liver (ug/g)
    
    #Plasma compartment
    dAplas_free = (QR*CVR*Free) + (QK*CVK*Free) + (QL*CVL*Free) - 
             (QC*CA*Free) + dAefflux  #rate of change in the plasma (mg/h) 
    
    #Mass Balance Check
    Atissue = Aplas_free + AR + AKb + Afil + APTC + AL + AST + ASI 	#sum of mass in all compartments (mg)
    Aloss = Aurine + Afeces #sum of mass lost through urinary and fecal excretion (mg)
    Atotal = Atissue + Aloss 	#total mass; should equal total dose
    
    list(c("dAR" = dAR, "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
           "dACl" = dACl, "dAefflux" = dAefflux,
           "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
           "dAurine" = dAurine, "dAST" =dAST,
           "dAabsST" = dAabsST, "dASI" = dASI, "dAabsSI" = dAabsSI, "dAfeces" = dAfeces, 
           "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free), 
         "total_oral_uptake" = total_oral_uptake, "amount_per_gram_liver" = amount_per_gram_liver,
         "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal)
    
  })
}

#=============
#6. User input 
#=============
# Parameters for reproducing example of sheet "MKudo1BW" in "PFAS_template_parameters_PFOA.xlsx" of Bernstein et al.2021
sex <- "M" # rat sex, values: M/F
BW <- 0.29# rat body weight in kg
admin.type <-  "iv" # administration type values: iv/oral
admin.dose <- 0.041 * BW # administered dose in mg
admin.time <- 0.001 # time when doses are administered, in hours

user_input <- list( "admin.type" = admin.type,
                    "admin.dose" = admin.dose, 
                    "admin.time" = admin.time,
                    "BW"=BW,  "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,3,0.1)  #hours

solution <-  ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                 events = events, method="bdf",rtol = 1e-05, atol = 1e-05)
print(tail(solution))

#====================
#7. Upload on Jaqpot 
#===================

# # Subset of features to be displayed on the user interface
# predicted.feats <- c("AR", "Adif" , "A_baso", "AKb", "ACl", "Aefflux" ,
#                      "A_apical" , "APTC" , "Afil" , "Aurine" , "AST",
#                      "AabsST", "ASI", "AabsSI", "Afeces", 
#                      "AL" , "Abile" , "Aplas_free")
# # Log in Jaqpot server
# jaqpotr::login.cred()
# 
# # Deploy the model on the Jaqpot server to create a web service
# jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
#                      create.params = create.params,  create.inits = create.inits,
#                      create.events = create.events, custom.func = custom.func, 
#                      method = "bdf",url = "https://api.jaqpot.org/jaqpot/")