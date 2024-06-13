ga_fitness <- function(chromosome) 
{ 
  setwd("C:/Users/user/Documents/GitHub/PBK_Grouping/PFOA")
    
  
  library(deSolve)
  library(ggplot2)
  
  #=========================
  #1. Parameters of the model
  #=========================
  create.params  <- function(user_input){
    with( as.list(user_input),{
      # Initialise a vector for correction factors
      CF <- rep(NA, length(group))
      # Initialise a vector for storing fitted parameters. The first value is 1, i.e. no change
      parameter_values <- c(1,rep(NA, N_pars))
      # Retrieve the fitted parameter information
      for (i in 1:N_pars){
        parameter_values[i+1] <- exp(fitted_pars[i])
      }
       
      #Pass grouping information from GA to the correction factors
      for (i in 1:10){
        CF[i] <- parameter_values[group[i]]
      }
      # User input: BW(kg), sex(F/M)
      #Physiological Parameters
      MKC = 0.0073	#fraction mass of kidney (percent of BW); Brown 1997
      MLC = 0.0366	#fraction mass of liver (percent of BW); Brown 1997
         
      #Cardiac Output and Bloodflow (as fraction of cardiac output)
      QCC = 14.0 #cardiac output in L/h/kg^0.75; Brown 1997		 
      QL_hepatic_arteryC = 0.021 #fraction blood flow to liver; Brown 1997	
      QKC = 0.141 #fraction blood flow to kidney; Brown 1997.
      QgonadsC = 0.500/53  	#fraction blood flow to gonads; from doi:10.1124/dmd.107.015644
      QintestineC = 0.451/2.58	#fraction blood flow to intestine; from doi:10.1007/bf02353860
      QspleenC = 0.63/74	#fraction blood flow to spleen; davies 1993
      QheartC = 0.051	#fraction blood flow to heart; Brown 1997
      QlungC = 	1#fraction blood flow to lung; Brown 1997
      QbrainC = 0.02	#fraction blood flow to brain; Brown 1997
      QstomachC = 0.068 /2.58	#fraction blood flow to stomach; from doi:10.1007/bf02353860
      Htc = 0.46 #hematocrit for the rat; Davies 1993
      
      #Tissue Volumes 
      VplasC = 0.0312 #fraction vol. of plasma (L/kg BW); Davies 1993
      VliverC = 0.0366 #fraction vol. of liver (L/kg BW); Brown 1997
      VkidneyC = 0.0073 #fraction vol. of kidney (L/kg BW); Brown 1997
      VgonadsC = 2.50/250	#fraction vol. of gonads (L/kg BW); from doi:10.1124/dmd.107.015644
      VintestineC = 0.0140+0.0084	#fraction vol. of tointestine (L/kg BW); Brown 1997
      VspleenC = 0.0020 #fraction vol. of spleen (L/kg BW); Brown 1997
      VheartC = 0.0033	#fraction vol. of heart (L/kg BW); Brown 1997
      VlungC = 	0.0050#fraction vol. of lung (L/kg BW); Brown 1997
      VbrainC = 0.0057	#fraction vol. of brain (L/kg BW); Brown 1997
      VstomachC = 0.0046	#fraction vol. of stomach (L/kg BW); Brown 1997
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
      Pliver = 2.2 #liver:blood
      Pkidney = 1.05 #kidney:blood
      Prest = 0.030/0.254 #rest of body:blood
      Pintestine = 0.017/0.254 
      Pgonads = 0.018/0.254 
      Pspleen = 0.015/0.254 
      Pheart = 0.029/0.254 
      Plung = 0.039/0.254
      Pbrain = 0.003/0.254
      Pstomach = 0.009/0.254 
      
      #rate constants
      kdif = 0.001	#diffusion rate from proximal tubule cells (L/h)
      kabsc = 2.12	#rate of absorption of chemical from small intestine to liver (1/(h*BW^0.25))(fit to data)
      kunabsc = 7.06e-5	#rate of unabsorbed dose to appear in feces (1/(h*BW^0.25))(fit to data)
      GEC = 0.54 #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013
      k0C = 1.0	#rate of uptake from the stomach into the liver (1/(h*BW^0.25)) (fit to data)
      
      keffluxc = 2.49 #rate of clearance of PFOA from proximal tubule cells into blood (1/(h*BW^0.25))
      kbilec = 0.004 #biliary elimination rate ((male); liver to feces storage (1/(h*BW^0.25)) (fit to data)
      kurinec = 1.6 #rate of urine elimination from urine storage (male) (1/(h*BW^0.25))(fit to data)
      Free = 0.09#0.006 #free fraction in plasma (male)
      
      if(sex == "M"){
        GFRC = 62.1	#glomerular filtration rate (L/hr/kg kidney) (male); Corley, 2005
        RAFbaso = 4.07	#relative activity factor, basolateral transporters (male) (fit to data)	
        RAFapi = 35	#relative activity factor, apical transporters (male) (fit to data)	
      }else if(sex == "F"){
        GFRC = 41.04	#glomerular filtration rate (L/hr/kg kidney) (female); Corley, 2005
        RAFapi = 0.001356	#relative activity factor, apical transporters (female) (fit to data)
        RAFbaso = 0.01356	#relative activity factor, basolateral transporters (female) (fit to data)
      }
      
      
      #Scaled Parameters
      #Cardiac output and blood flows
      QC <- QCC*(BW^0.75)*(1-Htc)	#cardiac output in L/h; adjusted for plasma
      Qlung <- QC
      QK <- (QKC*QC)	#plasma flow to kidney (L/h)
      QL_hepatic_artery <- (QL_hepatic_arteryC*QC)	#plasma flow to liver (L/h)
      Qspleen <- (QspleenC*QC)	#plasma flow to spleen (L/h)
      Qgonads <- (QgonadsC*QC)	#plasma flow to gonads (L/h)
      Qintestine <- (QintestineC*QC)	#plasma flow to intestine (L/h)
      Qstomach <- (QstomachC*QC)	#plasma flow to stomach (L/h)
      Qheart <- (QheartC*QC)	#plasma flow to heart (L/h)
      Qbrain <- (QbrainC*QC)	#plasma flow to liver (L/h)
      Qrest <- QC - QK - QL_hepatic_artery - Qspleen - Qgonads - Qintestine	- Qstomach -
               Qheart - Qbrain#plasma flow to rest of body (L/h)
      QBal <- QC - (QK + QL_hepatic_artery + Qspleen + Qgonads + Qintestine	+ Qstomach +
                      Qheart + Qbrain+Qrest) #Balance check of blood flows; should equal zero
      
      #Tissue Volumes
      MK = MKC*BW*1000	#mass of the kidney (g)
      PTC = MKC*6e7	#number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
      VPTC = MK*VPTCC	#volume of proximal tubule cells (L)	
      VPlas = VplasC*BW 	#volume of plasma (L) 
      V_venous_blood <- BW*11.3/250 	#volume of venous plasma (L); from doi:10.1007/bf02353860
      V_arterial_blood <- BW*5.6/250	#volume of arterial plasma (L); from doi:10.1007/bf02353860
      V_blood <-   V_venous_blood + V_arterial_blood
      Vven_plasma = VPlas*V_venous_blood/V_blood	#volume of venous plasma (L)
      Vart_plasma = VPlas*V_arterial_blood/V_blood  	#volume of arterial plasma (L)
      Vgonads = VgonadsC*BW 	#volume of gonads (L)
      Vspleen = VspleenC*BW 	#volume of spleen (L)
      Vheart = VheartC*BW 	#volume of heart (L)
      Vstomach = VstomachC*BW 	#volume of stomach (L)
      Vintestine = VintestineC*BW 	#volume of intestine (L)
      Vlung = VlungC*BW 	#volume of lung (L)
      Vbrain = VbrainC*BW 	#volume of brain (L)
      Vkidney = VkidneyC*BW 	#volume of kidney (L)
      Vliver = VliverC*BW	#volume of liver (L)
      
      Vkidneyb <- Vkidney*0.16	#volume of blood in the kidney (L); fraction blood volume of kidney (0.16) from Brown, 1997
      Vbrainb <- Vbrain*0.21	#volume of blood in the brain (L); fraction blood volume of brain (0.21) from Brown, 1997
      Vliverb <- Vliver*0.03	#volume of blood in the liver (L); fraction blood volume of liver (0.03) from Brown, 1997
      Vfil = VfilC*BW	#volume of filtrate (L)
      Vrest = BW - V_venous_blood -V_arterial_blood -Vfil  - Vliver - Vkidney - Vbrain - Vlung-
           Vintestine - Vstomach - Vheart - Vspleen -  Vgonads #volume of remaining tissue (L); 
      VBal = BW - (Vrest + Vliver + VPTC + Vfil + VPlas)	#Balance check of tissue volumes; should equal zero 
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
      
      return(list( "QC" = QC, "QK" = QK, "QL_hepatic_artery" = QL_hepatic_artery, "Qrest" = Qrest, 
                   "Qgonads" = Qgonads, "Qstomach" = Qstomach, "Qintestine" = Qintestine,
                   "Qbrain" = Qbrain, "Qlung" = Qlung,  "Qheart" = Qheart,
                   "Qspleen" = Qspleen,
                   
                   "VPlas" = VPlas, "Vart_plasma" = Vart_plasma, "Vven_plasma" = Vven_plasma,
                   "Vkidneyb" = Vkidneyb, "Vfil" = Vfil, "Vliver" = Vliver, "Vrest" = Vrest, 
                   "Vgonads" = Vgonads, "Vstomach" = Vstomach, "Vintestine" = Vintestine,
                   "Vbrain" = Vbrain, "Vlung" = Vlung,  "Vheart" = Vheart,
                   "Vspleen" = Vspleen, "Vkidney" = Vkidney,"Vbrainb" = Vbrainb, 
                   "Vliverb" = Vliverb,
                   
                   "GFR" = GFR, "VPTC" = VPTC,"Km_baso" = Km_baso, "Km_apical" = Km_apical,
                   "Vmax_apical" = Vmax_apical, "kbile" = kbile, "kurine" = kurine, 
                   "kunabs" = kunabs, "GE" = GE,
                   
                   "Pliver" = Pliver*CF[1],  "Prest" = Prest*CF[2], 
                   "Pintestine" = Pintestine*CF[3], "Pgonads" = Pgonads*CF[4],
                   "Pspleen" = Pspleen*CF[5], "Pheart" = Pheart*CF[6],
                   "Plung" = Plung*CF[7], "Pbrain" = Pbrain*CF[8], "Pstomach" =Pstomach*CF[9], 
                   "Pkidney" = Pkidney*CF[10],
                   
                   "Vmax_baso" = Vmax_baso*CF[11], 
                   'kdif' = kdif*CF[12], 
                  "kefflux" = kefflux*CF[13],
                   "kabs" = kabs*CF[14],  "k0" = k0*CF[15],
                   "Free" = Free*CF[16], 
                   
                  
                   "admin.type" = admin.type,
                   "admin.time" = admin.time, "admin.dose" = admin.dose))
      
    })
  }
  
  #===============================================
  #2. Function to create initial values for ODEs 
  #===============================================
  
  create.inits <- function(parameters){
    with( as.list(parameters),{
      "Arest" <- 0; "Adif" <- 0; "A_baso" <- 0; "Akidney" <- 0;
      "ACl" <- 0; "Aefflux" <- 0;
      "A_apical" <- 0; "APTC" <- 0; "Afil" <- 0;
      "Aurine" <- 0; "Astomach_lumen" <- 0;
      "Astomach" <- 0; "Aintestine_lumen" <- 0; "Aintestine" <- 0; "Afeces" <- 0;
      "Aliver" <- 0; "Abile" <- 0;
      "Agonads" <- 0; "Aspleen" <- 0; "Aheart" <- 0;
      "Alung" <- 0; "Abrain" <- 0; "Aven_free" <- 0; "Aart_free" <- 0;
      
      return(c("Arest" = Arest, "Agonads" = Agonads, "Aspleen" = Aspleen, "Aheart" = Aheart,
               "Alung" = Alung, "Abrain" = Abrain, "Adif" = Adif, "A_baso" = A_baso, 
               "Akidney" = Akidney,
               "ACl" = ACl, "Aefflux" = Aefflux,
               "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
               "Aurine" = Aurine, "Astomach_lumen" =Astomach_lumen,
               "Astomach" = Astomach, "Aintestine_lumen" = Aintestine_lumen, "Aintestine" = Aintestine,
               "Afeces" = Afeces, 
               "Aliver" = Aliver, "Abile" = Abile, "Aven_free" = Aven_free, "Aart_free" = Aart_free))
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
          events <- list(data = rbind(data.frame(var = c("Aven_free"),  time = admin.time, 
                                                 value = admin.dose, method = c("add")) ))
        }else if (admin.type == "oral"){
          events <- list(data = rbind(data.frame(var = c("Astomach_lumen"),  time = admin.time, 
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
      Crest = Arest/Vrest #concentration in rest of body (mg/L)
      CVrest = Crest/Prest	#concentration in venous blood leaving the rest of the body (mg/L)
      Ckidney = Akidney/Vkidney	#concentration in kidney blodd (mg/L) 
      CVkidney = Ckidney /Pkidney	#concentration in venous blood leaving kidney (mg/L)
      CPTC = APTC/VPTC	#concentration in PTC (mg/L)
      Cfil = Afil/Vfil	#concentration in filtrate (mg/L)
      Cliver = Aliver/Vliver	#concentration in the liver (mg/L)
      CVliver = Cliver/Pliver	#concentration in the venous blood leaving the liver (mg/L)
      Cgonads = Agonads/Vgonads	#concentration in the gonads (mg/L)
      CVgonads = Cgonads/Pgonads	#concentration in the venous blood leaving the gonads (mg/L)
      Cspleen = Aspleen/Vspleen	#concentration in the spleen (mg/L)
      CVspleen = Cspleen/Pspleen	#concentration in the venous blood leaving the spleen (mg/L)
      Cheart = Aheart/Vheart	#concentration in the heart (mg/L)
      CVheart = Cheart/Pheart	#concentration in the venous blood leaving the heart (mg/L)
      Clung = Alung/Vlung	#concentration in the lung (mg/L)
      CVlung = Clung/Plung	#concentration in the venous blood leaving the lung (mg/L)
      Cbrain = Abrain/Vbrain	#concentration in the brain (mg/L)
      CVbrain = Cbrain/Pbrain	#concentration in the venous blood leaving the brain (mg/L)
      Cintestine = Aintestine/Vintestine	#concentration in the liver (mg/L)
      CVintestine = Cintestine/Pintestine	#concentration in the venous blood leaving the liver (mg/L)
      Cstomach = Astomach/Vstomach	#concentration in the liver (mg/L)
      CVstomach = Cstomach/Pstomach	#concentration in the venous blood leaving the liver (mg/L)
    
      Cart_free = Aart_free/Vart_plasma		#concentration in arterial  plasma (mg)
      Cart = Cart_free/Free	#concentration of total PFOA in arterial plasma (mg/L)
      Cven_free = Aven_free/Vven_plasma		#concentration in venous plasma (mg)
      Cven = Cven_free/Free	#concentration of total PFOA in venous plasma (mg/L)
      
      # Rest of Body (Tis)
      dArest = Qrest*(Cart-CVrest)*Free	#rate of change in rest of body (mg/h)
      dAgonads = Qgonads*(Cart-CVgonads)*Free	#rate of change in rest of body (mg/h)
      dAspleen = Qspleen*(Cart-CVspleen)*Free	#rate of change in rest of body (mg/h)
      dAheart = Qheart*(Cart-CVheart)*Free	#rate of change in rest of body (mg/h)
      dAlung = Qlung*(Cven-CVlung)*Free	#rate of change in rest of body (mg/h)
      dAbrain = Qbrain*(Cart-CVbrain)*Free	#rate of change in rest of body (mg/h)
      
      #Kidney 
      #Kidney Blood (Kb)
      dAdif <- kdif*(Cven_free - CPTC)	#rate of diffusion from into the PTC (mg/hr)
      dA_baso <- (Vmax_baso*Cven_free)/(Km_baso + Cven_free)	
      dAkidney <- QK*(Cart-CVkidney)*Free  #rate of change in kidney blood (mg/h).
      dACl <- Cven*GFR*Free	#rate of clearance via glormerular filtration (mg/h)
     
      #Proximal Tubule Cells (PTC)
      dAefflux <- kefflux*APTC
      dA_apical <- (Vmax_apical*Cfil)/(Km_apical + Cfil)
      dAPTC <-  dAdif + dA_apical + dA_baso - dAefflux #rate of change in PTC(mg/h)
      
      #Filtrate (Fil)
      dAfil = Cart*GFR*Free - dA_apical - Afil*kurine	#rate of change in filtrate (mg/h)
      
      #Urinary elimination 
      dAurine = kurine*Afil	#rate of change in urine (mg/h)
  
      #GI Tract (Absorption site of oral dose)
      #Stomach
      dAstomach_lumen=  - k0*Astomach_lumen - GE*Astomach_lumen 	#rate of change in the stomach lumen (mg/h)
      dAstomach = k0*Astomach_lumen + Qstomach*(Cart-CVstomach)*Free	#rate change in the stomach (mg/h)
  
      #Small Intestine
      dAintestine_lumen = GE*Astomach_lumen - kabs*Aintestine_lumen - kunabs*Aintestine_lumen	#rate of change in the  intestine lumen (mg/hr)
      dAintestine = kabs*Aintestine_lumen +  Qintestine*(Cart-CVintestine)*Free	#rate change in the  intestine (mg/hr)
  
      #Feces compartment
      dAfeces = kbile*Aliver + kunabs*Aintestine_lumen #rate of change in the feces compartment (mg/h)
  
      #Liver
      dAliver = QL_hepatic_artery*Cart*Free - kbile*Aliver + kabs*Aintestine_lumen + k0*Astomach_lumen +
        Qspleen*CVspleen*Free +Qstomach*CVstomach*Free+ Qintestine*CVintestine*Free-
        (QL_hepatic_artery+Qspleen+Qstomach+Qintestine)*CVliver*Free#rate of change in the liver (mg/h)
      dAbile = kbile*Aliver  
  
      #Venous Plasma compartment
      dAven_free = Qrest*CVrest*Free + Qgonads*CVgonads*Free +  Qheart*CVheart*Free + Qbrain*CVbrain*Free +
                 (QK*CVkidney*Free) + ((QL_hepatic_artery+Qspleen+Qstomach+Qintestine) *CVliver*Free) - 
               (Qlung*Cven*Free) + dAefflux- dA_baso - dAdif  #rate of change in the plasma (mg/h) 
      
      #Arterial Plasma compartment
      dAart_free =  Qlung*CVlung*Free - Cart*GFR*Free- Cart*Free*(Qrest+Qgonads+Qspleen+Qheart+
                                                      Qbrain+QK+Qstomach+Qintestine+QL_hepatic_artery)
      
      #Mass Balance Check
      Atissue = Aart_free +Aven_free+ Arest + Akidney + Afil + APTC + Aliver + 
        Astomach + Astomach_lumen + Aintestine+ Aintestine_lumen+
        Agonads + Aspleen + Aheart + Alung + Abrain #sum of mass in all compartments (mg)
      Aloss = Aurine + Afeces #sum of mass lost through urinary and fecal excretion (mg)
      Atotal = Atissue + Aloss 	#total mass; should equal total dose
      
      list(c("dArest" = dArest, "dAgonads" = dAgonads, "dAspleen" = dAspleen, "dAheart" = dAheart,
             "dAlung" = dAlung, "dAbrain" = dAbrain, 
             "dAdif" = dAdif, "dA_baso" = dA_baso, "dAkidney" = dAkidney,
             "dACl" = dACl, "dAefflux" = dAefflux,
             "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
             "dAurine" = dAurine, "dAstomach_lumen" =dAstomach_lumen,
             "dAstomach" = dAstomach, "dAintestine_lumen" = dAintestine_lumen,
             "dAintestine" = dAintestine, "dAfeces" = dAfeces, 
             "dAliver" = dAliver, "dAbile" = dAbile, "dAven_free" = dAven_free,
             "dAart_free" = dAart_free), 
           "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal, "Crest" =Crest,
           "CVrest" = CVrest,
           "CVkidney" = CVkidney, "CPTC" = CPTC,
           "Cfil" = Cfil,  "CVliver" = CVliver, "Cart_free" = Cart_free,
           "Cart" = Cart, 
           "Cplasma" = Aven_free/VPlas/Free,
           "Cliver" = (Aliver + CVliver*Vliverb)/Vliver,
           "Ckidneys" = (APTC+Akidney+CVkidney*Vkidneyb)/Vkidney,
           "Cbrain" = (Abrain + CVbrain*Vbrainb)/Vbrain)
      
    })
  }
  
  #======================
  #3. Objective function  
  #======================

  obj.func <- function(fitted_pars, group, serum_male, serum_indices_male, 
                       tissue_male,tissue_indices_male, serum_female, serum_indices_female, 
                       tissue_female,tissue_indices_female, inits, N_pars){
    #Subroutine for returning the goodness-of-fit on the serum data, male or female
    solve_odes <- function(serum, admin.type, admin.dose, indices, 
                           BW, sex, fitted_pars, group, N_pars ){
    # Calculate PBK parameters
    parameters <- create.params( list("BW" = BW  , sex = sex, 
                                      "admin.type" = admin.type,
                                      "admin.time" = 0.01, "admin.dose" = admin.dose,
                                      "fitted_pars" = fitted_pars, "group" = group,
                                       "N_pars" = N_pars))
    events <- create.events(parameters)
    # Structure the in silico time vector in a way that the sampling time points are included
    if(sex == "F"){
      sample_time <- c(0,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5, 1e-4, 1e-3,
                       seq(0.01,0.08, 0.001), 0.083, 0.1, 0.2, 0.25, 
                       seq(0.3,  0.9, 0.1), seq(1,23,1),
                       seq(24,196 , 4))
    }else if(sex == "M"){
      sample_time <- c(0,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5, 1e-4, 1e-3,
                       seq(0.01,0.08, 0.001), 0.083, 0.1, 0.2, 0.25, 
                       seq(0.3,  0.9, 0.1), seq(1,23,1),
                       seq(24, 192, 4), seq(196, 516,16), 528, seq(540, 860,16),
                       864, seq(880, 1200,32))
    }else(
      stop(" Provide a valid sex. Chose between 'F' or 'M'")
    )
    solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                        y = inits, parms = parameters, events = events,
                                        method="lsodes",rtol = 1e-4, atol = 1e-4))
    
    concentrations <- data.frame("time" = solution$time, "Cplasma" = solution$Cplasma)
    experimental_time_points <- serum$Time[indices[1]:indices[2]]
    concentrations <- concentrations[concentrations$time %in% experimental_time_points, "Cplasma"]
    
    observed <- list("Cplasma" = serum$Mass[indices[1]:indices[2]])
    predicted <- list("Cplasma" = concentrations)
    #Calculate goodness-of-fit
    discrepancy <- SODI(observed, predicted)
    return(list("discrepancy" = discrepancy,"solution" = solution))
    }
    # Weight from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en based on SD rats 8 weeks old
    BW_male <- 0.3#kg
    BW_female <- 0.2#kg
    #######################
    # Goodness-of-fit on male serum
    # IV, 6mg/kg
    discrepancy_iv_6 <- solve_odes(serum = serum_male, admin.type = unique(serum_male$Type)[1],
                                    admin.dose = BW_male*unique(serum_male$Dose)[1], 
                                   indices = c(1,serum_indices_male[1]-1),
                                   BW = BW_male, sex = "M", fitted_pars = fitted_pars,
                                   group = group, N_pars = N_pars)$discrepancy
    # Oral, 6mg/kg
    discrepancy_oral_6 <- solve_odes(serum = serum_male,admin.type = unique(serum_male$Type)[2],
                                    admin.dose = BW_male*unique(serum_male$Dose)[1], 
                                    indices = c(serum_indices_male[1],serum_indices_male[2]-1),
                                                BW = BW_male, sex = "M", fitted_pars = fitted_pars,
                                                group = group, N_pars = N_pars)$discrepancy
    # Oral, 12mg/kg
    solution_oral_12 <- solve_odes(serum = serum_male,admin.type = unique(serum_male$Type)[2],
                                      admin.dose = BW_male*unique(serum_male$Dose)[2], 
                                   indices = c(serum_indices_male[2],serum_indices_male[3]-1),
                                   BW = BW_male, sex = "M", fitted_pars = fitted_pars,
                                   group = group, N_pars = N_pars)
    discrepancy_oral_12 <- solution_oral_12$discrepancy
    
    # Oral, 48mg/kg
    discrepancy_oral_48 <- solve_odes(serum = serum_male,admin.type = unique(serum_male$Type)[2],
                                       admin.dose = BW_male*unique(serum_male$Dose)[3], 
                                      indices = c(serum_indices_male[3],dim(serum_male)[1]),
                                      BW = BW_male, sex = "M", fitted_pars = fitted_pars,
                                      group = group, N_pars = N_pars)$discrepancy
    
    # #######################
    # # Goodness-of-fit on female serum
    # # IV, 40mg/kg
    # discrepancy_iv_40 <- solve_odes(serum = serum_female, admin.type = unique(serum_female$Type)[1],
    #                                admin.dose = BW_female*unique(serum_female$Dose)[1], 
    #                                indices = c(1,serum_indices_female[1]-1),
    #                                BW = BW_female, sex = "F", fitted_pars = fitted_pars,
    #                                group = group, N_pars = N_pars)$discrepancy
    # # Oral, 40mg/kg
    # discrepancy_oral_40 <- solve_odes(serum = serum_female, admin.type = unique(serum_female$Type)[2],
    #                                  admin.dose = BW_female*unique(serum_female$Dose)[1], 
    #                                  indices = c(serum_indices_female[1],serum_indices_female[2]-1),
    #                                  BW = BW_female, sex = "F", fitted_pars = fitted_pars,
    #                                  group = group, N_pars = N_pars)$discrepancy
    # # Oral, 80mg/kg
    # solution_oral_80 <- solve_odes(serum = serum_female, admin.type = unique(serum_female$Type)[2],
    #                                admin.dose = BW_female*unique(serum_female$Dose)[2], 
    #                                indices = c(serum_indices_female[2],serum_indices_female[3]-1),
    #                                BW = BW_female, sex = "F", fitted_pars = fitted_pars,
    #                                group = group, N_pars = N_pars)
    # discrepancy_oral_80 <- solution_oral_80$discrepancy
    # 
    # 
    # # Oral, 320mg/kg
    # discrepancy_oral_320 <- solve_odes(serum = serum_female, admin.type = unique(serum_female$Type)[2],
    #                                   admin.dose = BW_female*unique(serum_female$Dose)[3], 
    #                                   indices = c(serum_indices_female[3],dim(serum_female)[1]),
    #                                   BW = BW_female, sex = "F", fitted_pars = fitted_pars,
    #                                   group = group, N_pars = N_pars)$discrepancy
    # 
    #######################
    # Estimate the goodness-of-fit on the male tissues
    concentrations <- data.frame("time" = solution_oral_12$solution$time, 
                                 "Ckidneys" = solution_oral_12$solution$Ckidneys,
                                 "Cliver" =solution_oral_12$solution$Cliver,
                                 "Cbrain" = solution_oral_12$solution$Cbrain)
    
    experimental_time_points_liver <- tissue_male$Time[1:(tissue_indices_male[1]-1)]
    concentration_liver <- concentrations[concentrations$time %in% experimental_time_points_liver, "Cliver"]
    
    experimental_time_points_kidneys <- tissue_male$Time[tissue_indices_male[1]:(tissue_indices_male[2]-1)]
    concentration_kidneys <- concentrations[concentrations$time %in% experimental_time_points_kidneys, "Ckidneys"]
    
    experimental_time_points_brain<- tissue_male$Time[tissue_indices_male[2]:dim(tissue_male)[1]]
    concentration_brain <- concentrations[concentrations$time %in% experimental_time_points_brain, "Ckidneys"]
    
    observed <- list("Cliver" = tissue_male[tissue_male$Tissue == "Liver","Mass"],
                     "Ckidneys" = tissue_male[tissue_male$Tissue == "Kidneys","Mass"],
                     "Ckidneys" = tissue_male[tissue_male$Tissue == "Brain","Mass"])
            
    predicted <- list("Cliver" = concentration_liver, "Ckidneys" = concentration_kidneys,
                      "Cbrain" = concentration_brain)
    #Calculate goodness-of-fit
    discrepancy_tissues_male <- SODI(observed, predicted)
    
    #######################
    # # Estimate the goodness-of-fit on the female tissues
    # concentrations <- data.frame("time" = solution_oral_12$solution$time, 
    #                              "Ckidneys" = solution_oral_12$solution$Ckidneys,
    #                              "Cliver" =solution_oral_12$solution$Cliver,
    #                              "Cbrain" = solution_oral_12$solution$Cbrain)
    # 
    # experimental_time_points_liver <- tissue_female$Time[1:(tissue_indices_female[1]-1)]
    # concentration_liver <- concentrations[concentrations$time %in% experimental_time_points_liver, "Cliver"]
    # 
    # experimental_time_points_kidneys <- tissue_female$Time[tissue_indices_female[1]:(tissue_indices_female[2]-1)]
    # concentration_kidneys <- concentrations[concentrations$time %in% experimental_time_points_kidneys, "Ckidneys"]
    # 
    # experimental_time_points_brain<- tissue_female$Time[tissue_indices_female[2]:dim(tissue_female)[1]]
    # concentration_brain <- concentrations[concentrations$time %in% experimental_time_points_brain, "Ckidneys"]
    # 
    # observed <- list("Cliver" = tissue_female[tissue_female$Tissue == "Liver","Mass"],
    #                  "Ckidneys" = tissue_female[tissue_female$Tissue == "Kidneys","Mass"],
    #                  "Ckidneys" = tissue_female[tissue_female$Tissue == "Brain","Mass"])
    # 
    # predicted <- list("Cliver" = concentration_liver, "Ckidneys" = concentration_kidneys,
    #                   "Cbrain" = concentration_brain)
    # #Calculate goodness-of-fit
    # discrepancy_tissues_female <- SODI(observed, predicted)
    # 
    
    # Calculate total discrepancy
    total_discrepancy <- discrepancy_iv_6 + discrepancy_oral_6 + discrepancy_oral_12 + 
                         discrepancy_oral_48+ discrepancy_tissues_male
                         #  discrepancy_iv_40 +discrepancy_oral_40 + 
                         # discrepancy_oral_80 + discrepancy_oral_320+
                         #   + discrepancy_tissues_female
    return(total_discrepancy)
  }

    # SODI function the returns the SODI index described in Tsiros et al.2024
  # predictions: list of vectors containing the predicted data
  # names of the compartments
  
  SODI <- function(observed, predicted, comp.names =NULL){
    # Check if the user provided the correct input format
    if (!is.list(observed) || !is.list(predicted)){
      stop(" The observations and predictions must be lists")
    }
    # Check if the user provided equal length lists
    if (length(observed) != length(predicted)){
      stop(" The observations and predictions must have the same compartments")
    }
    Ncomp <- length(observed) # Number of compartments
    I <- rep(NA, Ncomp) # Compartment discrepancy index
    N_obs <- rep(NA, Ncomp) #Number of observations per compartment
    #loop over the compartments
    for (i in 1:Ncomp){
      Et <- 0 #relative error with observations
      St <- 0  #relative error with simulations
      N <- length(observed[[i]]) # number of observations for compartment i
      # Check if observations and predictions have equal length
      if(N != length(predicted[[i]])){
        stop(paste0("Compartment ",i," had different length in the observations and predictions"))
      }
      N_obs[i] <- N # populate the N_obs vector
      for (j in 1:N){
        # sum of relative squared errors (error = observed - predicted)
        Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
        St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
      }
      
      # root mean of the square of observed values
      RMEt <- sqrt(Et/N)
      # root mean of the square of simulated values
      RMSt <- sqrt( St/N)
      
      I[i] <- (RMEt + RMSt)/2   
    }
    # Total number of observations
    Ntot <- sum(N_obs)
    # Initialise the consolidated discrepancy index
    Ic <-0
    for (i in 1:Ncomp){
      # Give weight to compartments with more observations (more information)
      Ic <- Ic +  I[i]* N_obs[i]/Ntot
    }
    # Name the list of compartment discrepancy indices
    if ( !is.null(comp.names)){
      names(I) <- comp.names
    }else if (!is.null(names(observed))){
      names(I) <- names(observed)
    } else if (!is.null(names(predicted)) && is.null(comp.names) ){
      names(I) <- names(predicted)
    }
    return(Ic)
    #return(list(Total_index = Ic, Compartment_index= I))
  }
  #==============================
  #5. Decode chromosomes  
  #==============================
  # Function for decoding the GA output. Simply, we take the floor of the continuous number
  decode_ga_real <- function(real_num){ 
    CF <- rep(NA, length(real_num))
    # Grouping of correctin factors
    for (i in 1:length(CF)){
      CF[i] <-  floor(real_num[i])
    }
    
    return(CF)
  }
  
  #===============
  # Load data  
  #===============
  MW = 414.07	#PFOA molecular mass (g/mol)
  # Load raw data from paper Kreyling et al.2017, which are given in %ID/g tissue
  df_serum_male <- openxlsx::read.xlsx("serum_male.xlsx",  colNames = T, rowNames = F)
  df_tissue_male <- openxlsx::read.xlsx("tissue_male.xlsx", colNames = T, rowNames = F)
  df_serum_female <- openxlsx::read.xlsx("serum_female.xlsx",  colNames = T, rowNames = F)
  df_tissue_female <- openxlsx::read.xlsx("tissue_female.xlsx", colNames = T, rowNames = F)
  #Rename columns for easier handling
  colnames(df_serum_male) <- c("Time", "Mass", "Dose", "Type")
  colnames(df_tissue_male) <- c("Time", "Mass", "Dose", "Tissue")
  colnames(df_serum_female) <- c("Time", "Mass", "Dose", "Type")
  colnames(df_tissue_female) <- c("Time", "Mass", "Dose", "Tissue")
  
  #Transform microMolar to mg/L
  df_serum_male$Mass <- df_serum_male$Mass*MW/1000
  df_tissue_male$Mass <- df_tissue_male$Mass*MW/1000
  df_serum_female$Mass <- df_serum_female$Mass*MW/1000
  df_tissue_female$Mass <- df_tissue_female$Mass*MW/1000
  
  #Decode the chromosome of the genetic algorithm
  group <- decode_ga_real(chromosome)

  #Initialise optimiser to NULL for better error handling later
  optimizer <- NULL
  opts <- list( "algorithm" = "NLOPT_LN_NEWUOA",
                "xtol_rel" = 1e-07,
                "ftol_rel" = 0.0,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0 ,
                "maxeval" = 500)
  # Create initial conditions (zero initialisation)
  inits <- create.inits(list(NULL))
  N_pars <- 5 # Number of parameters to be fitted
  fit <- log(rep(1,N_pars))
  try(
    # Run the optimization algorithmm to estimate the parameter values
    optimizer <- nloptr::nloptr( x0= fit,
                                 eval_f = obj.func,
                                 lb	= rep(log(1e-4), N_pars),
                                 ub = rep(log(1000), N_pars),
                                 opts = opts,
                                 group = group,
                                 serum_male = df_serum_male,
                                 serum_indices_male = c(13,23,33),# index where dose changes
                                 tissue_male = df_tissue_male,
                                 tissue_indices_male = c(8,14),# index where tissue changes
                                 serum_female = df_serum_female,
                                 serum_indices_female = c(10,19,29),
                                 tissue_female = df_tissue_female,
                                 tissue_indices_female = c(6,11),
                                 inits = inits,
                                 N_pars = N_pars), 

    silent = TRUE
  )
  
  # If the otpimizer didn't run, set a small predefined value for the objective function value
  if(is.null(optimizer)){
    of_value <- -1000 
  }else{
    of_value <- -optimizer$objective
  }
  
  
  return(of_value)
}



##############################
#=============================
#  ***  Genetic algorithm  ***
#=============================
##############################

#=======================================================================
#                    Available tuning parameters                       
#                        (real encoding)                             
#=======================================================================

#                            /Selection/          
# gareal_lrSelection:Linear-rank selection
# gareal_nlrSelection:Nonlinear-rank selection.
# gareal_rwSelection:Proportional (roulette wheel) selection.
# gareal_tourSelection: (Unbiased) tournament selection
# gareal_lsSelection: Fitness proportional selection with fittness linear scaling.
# gareal_sigmaSelection: Fitness proportional selection with Goldberg's sigma truncation scaling
#
#                            /Crossover/                               
# gareal_spCrossover: Single-point crossover
# gareal_uCrossover: Uniform crossover
# gareal_waCrossover: Whole arithmetic crossover.
# gareal_laCrossover: Local arithmetic crossover.
# gareal_blxCrossover: Blend crossover.
#
#                           /Mutation/
# gareal_raMutation: Uniform random mutation
# gareal_nraMutation: Nonuniform random mutation.
# gareal_rsMutation: Random mutation around the solution.
setwd("C:/Users/user/Documents/GitHub/PBK_Grouping/PFOA")
N_genes <- 16 #number of parameters to be included in the grouping process
N_part <- 10 #number of partition coefficients
N_pars <- 5 # Number of parameters to be fitted
start <- Sys.time()
GA_results <- GA::ga(type = "real", fitness = ga_fitness, 
                     lower = c(rep(1,N_part), rep(5,(N_genes-N_part))),
                     upper = c(rep(4.99999,N_geneN_parts),rep(7.99999,(N_genes-N_part)) ), 
                     population = "gareal_Population",
                     selection = "gareal_lsSelection",
                     crossover = "gareal_laCrossover", 
                     mutation = "gareal_raMutation",
                     popSize =  5*parallel::detectCores(), #the population size.
                     pcrossover = 0.85, #the probability of crossover between pairs of chromosomes.
                     pmutation = 0.4, #the probability of mutation in a parent chromosome
                     elitism = 5, #the number of best fitness individuals to survive at each generation. 
                     maxiter = 200, #the maximum number of iterations to run before the GA search is halted.
                     run = 50, # the number of consecutive generations without any improvement
                     #in the best fitness value before the GA is stopped.
                     keepBest = TRUE, # best solutions at each iteration should be saved in a slot called bestSol.
                     parallel = (parallel::detectCores()),
                     monitor =plot,
                     seed = 1234)
stop <- Sys.time()
print(paste0("Time ellapsed was ", stop-start))
save.image(file = "PFOA.RData")
