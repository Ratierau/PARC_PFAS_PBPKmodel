create.params <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al)
    # fractional bloo flows to tissues
    QCC = 12.5		  # Cariac blood output (L/h/kg^0.75)
    QFC = 0.061		  # Fraction cardiac output going to fat
    QLC = 0.189 #0.25 from Loccisano 2011 		  # Fraction cariac output going to liver
    QKC = 0.177		  # Fraction cardiac output going to kiney
    QFilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kiney blood flow)
    QGC = 0.181		  # Fraction cardiac output going to gut
    QLuC = 0.034    # Fraction cardiac output going to lungs
    QBC = 0.117    # Fraction cardiac output going to brain
    
    
    #fractional tissue volumes
    VLC = 0.023		  # Fraction liver volume
    VFC = 0.02 #0.214		  # Fraction fat volume
    VKC = 0.004		  # Fraction kiney volume
    VFilC = 0.00042	# Fraction filtrate compartment volume (10% of kiney volume)
    VGC = 0.016		# Fraction gut volume
    VPlasC = 0.039	# Fraction plasma volume
    VLuC = 0.014	# Fraction lungs volume
    VBC = 0.021	# Fraction lungs volume
    
    # Scaling parameters
    QC = 312#QCC*BW**0.75	#Cardiac output (L/h)
    Htc = 0#0.44        #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/h)
    QF = QFC*QCP			# Plasma flow to fat (L/h)
    QK = QKC*QCP	    # Plasma flow to kiney (L/h)
    QFil = QFilC*QCP   	# Plasma flow to filtrate compartment (L/h)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/h)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/h)
    QB = QBC*QCP	    # Plasma flow to brain (L/h)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB	 # Plasma flow to rest of the boy (L/h)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kiney volume (L)
    VFil = VFilC*BW	    # Fitrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VPlas		# Rest of the boy volume (L) # Loccisano 2011 uses 0.84*BW
    
    Vbal = (1*BW)-(VL+VF+VK+VFil+VG+VPlas+VR+VLu+VB)         # Balance check
    
    
    if(substance == 'PFOA'){
      Tm = 147.4#6*BW^0.75 #ug/h
      Kt = 0.116 #ug/L
      Free = 0.03 # Fabrega (2014) Table 1
      PL = 1.03 # Partition Coefficient for Liver: Fabrega (2014) Table 1
      PF = 0.467; # Partition Coefficient for Fat: Fabrega (2014) Table 1
      PB = 0.17; # Partition Coefficient for Brain: Fabrega (2014) Table 1
      PLu = 1.27; # Partition Coefficient for Lungs: Fabrega (2014) Table 1
      PK = 1.17 # Partition Coefficient for Kiney: Fabrega (2014) Table 1
      PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
      PR = 0.12 # Partition Coefficient for Rest of boy: Loccisano (2011) Table 1
      kurinec = 3e-04	#urinary elimination rate constant  (/h/kg^-0.25); 
      kurine = kurinec*BW**(-0.25) # Elimination rate (1/h)
      Total_hourly_Intake = 0.218/24 # ug/h
    }else if(substance == 'PFOS'){
      Tm = 86#3.5*BW^0.75 #ug/h
      Kt = 0.0176 #ug/L
      Free = 0.03 # Fabrega (2014) Table 1
      PL = 2.67 # Partition Coefficient for Liver: Fabrega (2014) Table 1
      PF = 0.33 # Partition Coefficient for Fat: Fabrega (2014) Table 1
      PB = 0.255 # Partition Coefficient for Brain: Fabrega (2014) Table 1
      PLu = 0.155 # Partition Coefficient for Lungs: Fabrega (2014) Table 1
      PK = 1.26 # Partition Coefficient for Kiney: Fabrega (2014) Table 1
      PG = 0.57 # Partition Coefficient for Gut: Loccisano (2011) Table 1
      PR = 0.2 # Partition Coefficient for Rest of boy: Loccisano (2011) Table 1
      kurinec = 1e-3	#urinary elimination rate constant  (/h/kg^-0.25); estimated from Harada, 
      kurine = kurinec*BW**(-0.25) # Elimination rate (1/h)
      Total_hourly_Intake = 0.252/24 # ug/h
    }
    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'QFil'=QFil, 'QG'=QG, 'QLu'=QLu, 'QB'=QB, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'VFil'=VFil, 'VG'=VG, 'VLu'=VLu, 'VB'=VB, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PB'=PB, 'PLu'=PLu, 'PK'=PK, 
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free, 'kurine'=kurine,
                "Total_hourly_Intake"=Total_hourly_Intake))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    APlas<-0; AG<-0; AL<-0; AF<-0; ALu<-0; AB<-0
    AK<-0; AFil<-0; AStor<-0; AUrine<-0
    AR<-0; Intake<-0; A_in<-0
    
    return(c("APlas"=APlas, "AG"=AG, "AL"=AL, "AF"=AF, "ALu"=ALu, "AB"=AB,
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
    CK <- AK/VK # Concentration in Kiney
    CFil <- AFil/VFil # Concentration in Filtration compartment
    CR <- AR/VR  # Concentration in Rest of the boy compartment
    
    # Plasma compartment (ng)
    dAPlas <- QF*Free*CF/PF + (QL+QG)*Free*CL/PL + QR*Free*CR/PR + QK*Free*CK/PK + 
      QB*Free*CB/PB + QLu*Free*CLu/PLu - QCP*Free*CPlas 
    
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
    
    # Kiney compartment
    dAK <- QK*Free*(CPlas - CK/PK) + (Tm*CFil)/(Kt+CFil)
    
    # Filtrate compartment
    dAFil = QFil*(Free*CPlas - CFil) - (Tm*CFil)/(Kt+CFil)
    
    # Storage compartment for urine
    dAStor <- QFil*CFil - kurine*AStor
    
    # Amount excrete via urine
    dAUrine <- kurine*AStor
    
    # Rest of the boy
    dAR <- QR*Free*(CPlas - CR/PR)
    
    dIntake <- 0 
    
    dA_in <- Intake
    Mass_balance <- A_in - (APlas+AG+AL+AF+ALu+AB+AK+AFil+AStor+AR + AUrine)
    
    return(list(c("dAPlas"=dAPlas, "dAG"=dAG, "dAL"=dAL, "dAF"=dAF, "dALu"=dALu, 
                  "dAB"=dAB, "dAK"=dAK, "dAFil"=dAFil, "dAStor"=dAStor, "dAUrine"=dAUrine,
                  "dAR"=dAR, "dIntake"=dIntake, "dA_in"=dA_in),
                "CPlas"=CPlas, "CG"=CG, "CL"=CL, "CF"=CF, "CLu"=CLu, "CB"=CB,
                "CK"=CK, "CFil"=CFil, "CR"=CR, "Mass_Balance"=Mass_balance))
  })
}

########################################
BW <- 71.4 # kg
substance <- 'PFOS'
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

print(tail(solution,1))

# plot(solution$time/24/360, solution$CL)
# plot(solution$time/24/360, solution$Mass_Balance)
plot(solution$time/24/360, solution$CFil)
