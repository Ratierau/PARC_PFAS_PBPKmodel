setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/Fabrega et al.2015')
create.params <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al)
    # fractional blood flows to tissues
    QCC = 12.5		  # Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		  # Fraction cardiac output going to fat
    QLC = 0.19 #0.25 from Loccisano 2011 		  # Fraction cardiac output going to liver
    QKC = 0.175		  # Fraction cardiac output going to kiney
    #QFilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kiney blood flow)
    QGC = 0.181		  # Fraction cardiac output going to gut
    QLuC = 0.034    # Fraction cardiac output going to lungs
    QBC = 0.117    # Fraction cardiac output going to brain
    QBmC = 0.05    # Fraction cardiac output going to bone marrow
    
    #fractional tissue volumes
    VLC = 0.026		  # Fraction liver volume
    VFC = 0.02 #0.214		  # Fraction fat volume
    VKC = 0.004		  # Fraction kidney volume
    VFilC = 0.0004	# Fraction filtrate compartment volume (10% of kidney volume)
    VGC = 0.0171		# Fraction gut volume
    VPlasC = 0.0428	# Fraction plasma volume
    VLuC = 0.014	# Fraction lungs volume
    VBC = 0.021	# Fraction lungs volume
    VBmC = 0.1429	# Fraction lungs volume
    
    # Scaling parameters
    QC = QCC*BW**0.75	#Cardiac output (L/h)
    Htc = 0#0.44        #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/h)
    QF = QFC*QCP			# Plasma flow to fat (L/h)
    QK = QKC*QCP	    # Plasma flow to kiney (L/h)
    QFil = 0.2*QK   	# Plasma flow to filtrate compartment (L/h)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/h)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/h)
    QB = QBC*QCP	    # Plasma flow to brain (L/h)
    QBm = QBmC*QCP	  # Plasma flow to bone marrow (L/h)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB - QBm	 # Plasma flow to rest of the boy (L/h)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB+QBm)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kiney volume (L)
    VFil = VFilC*BW	    # Fitrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VBm = VBmC*BW			    # Bone marrow volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VBm - VPlas		# Rest of the boy volume (L) # Loccisano 2011 uses 0.84*BW
    
    Vbal = (1*BW)-(VL+VF+VK+VFil+VG+VPlas+VR+VLu+VB+VBm)         # Balance check
    
    loaded_params <- openxlsx::read.xlsx('Fabrega_2015_data.xlsx')
    keep_params <- loaded_params[loaded_params$Substance==substance,2:dim(loaded_params)[2]]
    
    Tm = keep_params$Tm #ug/h Fabrega (2015) Table 2
    Kt = keep_params$Kt #ug/L Fabrega (2015) Table 2
    Free = keep_params$Free #unitless Fabrega (2015) Table 2
    PL = keep_params$Liver # Partition Coefficient for Liver: Fabrega (2015) Table 2
    PF = 0.467 # Partition Coefficient for Fat: Fabrega (2014) Table 1
    PB = keep_params$Brain # Partition Coefficient for Brain: Fabrega (2015) Table 2
    PBm = keep_params$Bone_marrow # Partition Coefficient for Bone marrow: Fabrega (2015) Table 2
    PLu = keep_params$Lung # Partition Coefficient for Lung: Fabrega (2015) Table 2
    PK = keep_params$Kidney # Partition Coefficient for Kidney: Fabrega (2015) Table 2
    PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
    PR = 0.12 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
    kurinec = 3e-04	#urinary elimination rate constant  (/h/kg^-0.25); 
    kurine = kurinec*BW**(-0.25) # Elimination rate (1/h)
    Total_hourly_Intake = keep_params$Intake/1000/24 # ug/h

    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'QFil'=QFil, 'QG'=QG, 'QLu'=QLu, 'QB'=QB, 'QBm'=QBm, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'VFil'=VFil, 'VG'=VG, 'VLu'=VLu, 'VB'=VB, 'VBm'=VBm, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PB'=PB, 'PBm'=PBm, 'PLu'=PLu, 'PK'=PK, 
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free, 'kurine'=kurine,
                "Total_hourly_Intake"=Total_hourly_Intake))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    APlas<-0; AG<-0; AL<-0; AF<-0; ALu<-0; AB<-0; ABm<-0
    AK<-0; AFil<-0; AStor<-0; AUrine<-0
    AR<-0; Intake<-0; A_in<-0
    
    return(c("APlas"=APlas, "AG"=AG, "AL"=AL, "AF"=AF, "ALu"=ALu, "AB"=AB, "ABm"=ABm,
             "AK"=AK, "AFil"=AFil, "AStor"=AStor, "AUrine"=AUrine,
             "AR"=AR, "Intake"=Intake, "A_in"=A_in))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    ltimes <- length(exposure_times)
    lexposure <- length(Total_hourly_Intake)
    
    events <- data.frame(var = rep('Intake', ltimes), time = exposure_times,
                         value = Total_hourly_Intake, method = rep('rep',ltimes))
    return(list(data=events))
  })
}


ode.func <- function(time, inits, params){
  with(as.list(c(inits,params)),{
    
    # Units
    # C_i: Concentration in tissue i (ug/kg) #cosidering density = 1g/ml
    # Intake: hourly ingestion of PFASs (ug/h)
    
    # Concentrations (ng/g)
    CPlas <- APlas/VPlas #total concentration of chemical in plasma
    CG <- AG/VG # Concentration in Gut
    CL <- AL/VL # Concentration in Liver
    CF <- AF/VF # Concentration in Fat
    CLu <- ALu/VLu # Concentration in Lungs
    CB <- AB/VB # Concentration in Brain
    CBm <- ABm/VBm # Concentration in Bone marrow
    CK <- AK/VK # Concentration in Kidney
    CFil <- AFil/VFil # Concentration in Filtration compartment
    CR <- AR/VR  # Concentration in Rest of the body compartment
    
    # Plasma compartment (ng)
    dAPlas <- QF*Free*CF/PF + (QL+QG)*Free*CL/PL + QR*Free*CR/PR + QK*Free*CK/PK + 
      QB*Free*CB/PB + QBm*Free*CBm/PBm + QLu*Free*CLu/PLu - QCP*Free*CPlas 
    
    # Gut compartment
    dAG <- QG*Free*(CPlas - CG/PG) + Intake
    
    # Liver compartment
    dAL <- QL*Free*CPlas + QG*Free*CG/PG - (QL+QG)*Free*CL/PL
    
    # Fat compartment
    dAF <- QF*Free*(CPlas - CF/PF)
    
    # Lungs Compartment
    dALu <- QLu*Free*(CPlas - CLu/PLu)
    
    # Brain Compartment
    dAB <- QB*Free*(CPlas - CB/PB)

    # Bone Marrow Compartment
    dABm <- QBm*Free*(CPlas - CBm/PBm)
    
    # Kidney compartment
    dAK <- QK*Free*(CPlas - CK/PK) + (Tm*CFil)/(Kt+CFil)
    
    # Filtrate compartment
    dAFil = QFil*(Free*CPlas - CFil) - (Tm*CFil)/(Kt+CFil)
    
    # Storage compartment for urine
    dAStor <- QFil*CFil - kurine*AStor
    
    # Amount excrete via urine
    dAUrine <- kurine*AStor
    
    # Rest of the body
    dAR <- QR*Free*(CPlas - CR/PR)
    
    dIntake <- 0 
    
    dA_in <- Intake
    Mass_balance <- A_in - (APlas+AG+AL+AF+ALu+AB+ABm+AK+AFil+AStor+AR + AUrine)
    
    return(list(c("dAPlas"=dAPlas, "dAG"=dAG, "dAL"=dAL, "dAF"=dAF, "dALu"=dALu, 
                  "dAB"=dAB, "dABm"=dABm, "dAK"=dAK, "dAFil"=dAFil, "dAStor"=dAStor, "dAUrine"=dAUrine,
                  "dAR"=dAR, "dIntake"=dIntake, "dA_in"=dA_in),
                "CPlas"=CPlas, "CG"=CG, "CL"=CL, "CF"=CF, "CLu"=CLu, "CB"=CB, "CBm"=CBm,
                "CK"=CK, "CFil"=CFil, "CR"=CR, "Mass_Balance"=Mass_balance))
  })
}

########################################
BW <- 70 # kg
substance <- 'PFOA'
exposure <- 0
exposure_times <- 0
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "exposure"=exposure,
                   "exposure_times"= exposure_times)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
Lifetime <- 90 #years
sample_time <- seq(0,Lifetime*360*24,24)

solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events, 
                           method="lsodes",rtol = 1e-05, atol = 1e-05)) 

paste0('Simulation for ',substance)
print(tail(solution,1))

# plot(solution$time/24/360, solution$CL)
# plot(solution$time/24/360, solution$Mass_Balance)
#plot(solution$time/24/360, solution$CFil)
