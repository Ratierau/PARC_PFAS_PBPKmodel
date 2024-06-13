library(deSolve)
setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/Loccisano et al.2011/Monkey_model')

create.params <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al 1997)
    #fractional blood flows
    QCC = 15 		# Cardiac blood output (L/h/kg^0.75)
    QFC = 0.02		# Fraction cardiac output going to fat
    QLC = 0.2 		# Fraction cardiac output going to liver
    QKC = 0.123		# Fraction cardiac output going to kidney
    #QfilC = 0.0246		# Fraction cardiac output to the filtrate compartment (20% of kidney blood flow)
    QSkC = 0.048	# Fraction cardiac output going to skin
    QGC = 0.096		# Fraction cardiac output going to gut
    
    #fractional tissue volumes
    VLC = 0.026		# Fraction liver volume
    VFC = 0.05		# Fraction fat volume
    VKC = 0.004		# Fraction kidney volume
    VfilC = 0.0004	# Fraction filtrate compartment volume (10% of kidney volume)
    VGC = 0.017#0.0138		# Fraction gut volume
    VPlasC = 0.0448	# Fraction plasma volume (58% of blood)
    Htc = 0.41    # hematocrit
    
    # for dermal exposure
    # Dermal exposure
    Dermconc <- 0.0		# Dermal concentration (mg/mL)
    Dermvol <- 0.0		  # Dermal exposure volume (mL)
    Dermdose <- Dermconc*Dermvol*1000  # (ug)
    Skinarea = 5		#Exposed area on skin (cm^2)
    SkinTarea = 9.1*((BW*1000)**0.666)	  # Total area of skin (cm^2)
    Skinthickness = 0.1	                  # Skin thickness (cm)
    
    # Scaling parameters
    QC = QCC*BW**0.75	# Cardiac output (L/h)
    QCP = QC*(1-Htc)	# adjust for plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/h)
    QF = QFC*QCP			# Plasma flow to fat (L/h)
    QK = QKC*QCP	    # Plasma flow to kidney (L/h)
    Qfil = 0.2*QK		  # Plasma flow to filtrate compartment (L/h)# 20% of QK
    QG = QGC*QCP		  # Plasma flow to gut (L/h)
    QSk <- ifelse(Dermconc > 0, QSkC*QCP*(Skinarea/SkinTarea), 0.0) # plasma flow to skin
    #QSk <- QSkC*QCP*(Skinarea/SkinTarea) # plasma flow to skin
    QR = QCP - QL - QF - QK - Qfil - QG - QSk	# Plasma flow to rest of the body (L/h)
    
    Qbal = QCP - (QL+QF+QK+Qfil+QG+QSk+QR)       # balance check--better be 0
    
    VL = VLC*BW			# Liver volume (L)
    VF = VFC*BW			# Fat volume (L)
    VK = VKC*BW			# Kidney volume (L)
    Vfil = VfilC*BW			# Fitrate compartment volume (L)
    VG = VGC*BW			# Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VSk = (Skinarea*Skinthickness)/1000	                                    # Skin volume (L)
    VR = 0.84*BW - VL - VF - VK - Vfil - VG - VPlas - VSk	          # Rest of the body volume (L)
    Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+VR)               # Balance check--better be 0
    
    if(substance == 'PFOA'){
      # Chemical-specific parameters (PFOS)
      Tmc = 0.15			# Maximum resorption rate
      Kt = 0.055		# Resorption affinity 
      Free = 0.02	# Free fraction of PFOS in plasma
      PL = 2.2			# Liver/blood partition coefficient
      PF = 0.04			# Fat/blood partition coefficient
      PK = 1.05			# Kidney/blood partition coefficient
      PSk = 0.1		# Skin/blood partition coefficient
      PR = 0.12			# Rest of the body/blood partition coefficient
      PG = 0.05     # Gut/blood partition coeff. 
      kurinec = 50		# urinary elimination rate constant  (/h/kg^-0.25)# estimated from Harada, et al 2005
    } else if(substance == 'PFOS'){
      # Chemical-specific parameters (PFOA)
      Tmc = 1.3			# Maximum resorption rate
      Kt = 0.055#0.023		# Resorption affinity# same as monkey
      Free = 0.025		# Free fraction of PFOA in plasma# same as monkey
      PL = 3.72			# Liver/plasma partition coefficient
      PF = 0.14			# Fat/plasma partition coefficient
      PK = 1.05#0.8			# Kidney/plasma partition coefficient
      PSk = 0.29			# Skin/plasma partition coefficient
      PR = 0.12#0.2			# Rest of the body/plasma partition coefficient
      PG = 0.57     # Gut/blood plasma coeff.
      kurinec = 0.003		# Elimination rate (1/h)# estimated from data of Harada, et al 2005
    }
    
    kurine = kurinec*BW**(-0.25)
    Tm = Tmc*BW**0.75   #transporter maximum
    
    # Free fraction of chemical in tissues
    FreeL = Free/PL  #liver
    FreeF = Free/PF  #fat
    FreeK = Free/PK  #kidney
    FreeSk = Free/PSk #skin
    FreeR = Free/PR  #rest of tissues
    FreeG = Free/PG  #gut
    
    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'Qfil'=Qfil, 'QG'=QG, 'QSk'=QSk, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'Vfil'=Vfil, 'VG'=VG, 'VSk'=VSk, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PK'=PK, 'PSk'=PSk,
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free,
                'FreeL'=FreeL, 'FreeF'=FreeF, 'FreeK'=FreeK,
                'FreeSk'=FreeSk, 'FreeR'=FreeR, 'FreeG'=FreeG, 
                'kurine'=kurine,
                'C_water'=C_water, 
                'C_water_times'=C_water_times,
                'IV_intake_hourly'=IV_intake_hourly,
                'IV_times'=IV_times,
                'Oral_intake_hourly'=Oral_intake_hourly,
                'Oral_times'=Oral_times,
                'IV_doses'=IV_doses,
                'IV_doses_times'=IV_doses_times,
                'Drinking_rate'=Drinking_rate
    ))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    APlas<-0; AG<-0; AL<-0; AK<-0; AF<-0; Afil<-0; AStore<-0; AUrine<-0; ASk<-0;
    AR<-0; C_water<-0; IV_intake<-0; Oral_intake<-0
    return(c('APlas'=APlas, 'AG'=AG, 'AL'=AL, 'AF'=AF, 'AK'=AK, 
             'Afil'=Afil, 'AStore'=AStore, 'AUrine'=AUrine, 'ASk'=ASk, 
             'AR'=AR, 'C_water'=C_water, 'IV_intake'=IV_intake, 'Oral_intake'=Oral_intake))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    C_water_times <- C_water_times #hours
    C_water_values <- C_water #ug/L
    C_water_ltimes <- length(C_water_times)
    
    IV_intake_hourly <- IV_intake_hourly
    IV_times <- IV_times
    IV_ltimes <- length(IV_times)
    
    Oral_intake_hourly <- Oral_intake_hourly
    Oral_times <- Oral_times
    Oral_ltimes <- length(Oral_times)
    
    IV_doses <- IV_doses
    IV_doses_times<-IV_doses_times
    IV_doses_ltimes <- length(IV_doses_times)
    
    events <- data.frame(var = c(rep('C_water', C_water_ltimes), rep('IV_intake', IV_ltimes), rep('Oral_intake', Oral_ltimes), rep('APlas', IV_doses_ltimes)), 
                         time = c(C_water_times,IV_times, Oral_times, IV_doses_times),
                         value = c(C_water_values,IV_intake_hourly, Oral_intake_hourly, IV_doses),
                         method = rep('rep',C_water_ltimes+IV_ltimes+Oral_ltimes, IV_doses_ltimes))
    events <- events[order(events$time),] 
    return(list(data=events))
  })
}

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Concentrations
    CPlas <- APlas/VPlas # Concentration in plasma 
    CG <- AG/VG # Concentration in gut 
    CL <- AL/VL # Concentration in liver
    CF <- AF/VF # Concentration in fat
    CK <- AK/VK # Concentration in kidney
    Cfil <- Afil/Vfil # Concentration in filtrate
    CSk <- ASk/VSk # Concentration in skin
    CR <- AR/VR
    
    # Plasma Compartment
    dAPlas <- QF*CF*FreeF + (QL+QG)*CL*FreeL + QR*CR*FreeR + QSk*CSk*FreeSk + 
      QK*CK*FreeK - QCP*CPlas*Free + IV_intake
    
    dIV_intake <- 0
    dOral_intake <- 0
      
    # Gut compartment
    dAG <- QG*(CPlas*Free - CG*FreeG) + C_water*Drinking_rate + Oral_intake
    
    # Liver compartment
    dAL <- QL*CPlas*Free + QG*CG*FreeG - (QL+QG)*CL*FreeL 
    
    # Fat compartment
    dAF <- QF*(CPlas*Free - CF*FreeF)
    
    # Kidney compartment
    dAK <- QK*(CPlas*Free - CK*FreeK) + (Tm*Cfil)/(Kt+Cfil)
    
    # Filtrate compartment
    dAfil = Qfil*(CPlas*Free - Cfil) - (Tm*Cfil)/(Kt+Cfil)
    
    # Storage compartment
    dAStore <- Qfil*Cfil - kurine*AStore
    
    # Urine
    dAUrine <- kurine*AStore
    
    # Skin compartment
    dASk <- QSk*(CPlas*Free-CSk*FreeSk) #+ input4*DoseOn
    
    # Rest of the body
    dAR <- QR*(CPlas*Free - CR*FreeR)
    
    # Water concentration
    dC_water <- 0
    
    return(list(c('dAPlas'=dAPlas, 'dAG'=dAG, 'dAL'=dAL, 'dAF'=dAF,
                  'dAK'=dAK, 'dAfil'=dAfil, 'dAStore'=dAStore,
                  'dAUrine'=dAUrine, 'dASk'=dASk, 'dAR'=dAR, 'dC_water'=dC_water,
                  'dIV_intake'=dIV_intake, 'dOral_intake'=dOral_intake),
                'CPlas'=CPlas, 'CG'=CG, 'CL'=CL, 'CF'=CF,
                'CK'=CK, 'Cfil'=Cfil, 'CSk'=CSk, 'CR'=CR))
  })
}

################################################################################

reproduce_results <- function(user_input){
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                             events = events,
                             method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  plasma_plot <- plot(x=solution$time/24, y=solution$APlas)#, log="y")
  plasma_plot <- recordPlot(plasma_plot)
  
  return(list(plasma_plot, max(solution[,'CPlas'])))
}

BW = 3.5 # kg
substance <- 'PFOA'
C_water <- c(0) # ug/L
C_water_times <- c(0)
IV_intake_hourly <- c(0) # ug/h
IV_times <- c(0)
Oral_intake_hourly <- c(10*BW/24, 0) # mg/h
Oral_times <- c(0,26*7*24)
IV_doses <- c(0)*BW
IV_doses_times <- c(0)

Drinking_rate <- 0#BW*11/1000/24 # L/h
sample_time <- seq(0,18*30*24,1)
user_input <- list('BW'=BW, 'substance'=substance,
                   'C_water'=C_water,
                   'C_water_times'=C_water_times,
                   'IV_intake_hourly'=IV_intake_hourly,
                   'IV_times'=IV_times,
                   'Oral_intake_hourly'=Oral_intake_hourly,
                   'Oral_times'=Oral_times,
                   'IV_doses'=IV_doses, 
                   'IV_doses_times'=IV_doses_times,
                   'Drinking_rate'=Drinking_rate,
                   'sample_time'=sample_time)

fig4 = reproduce_results(user_input)
fig4
########################################################

BW = 3.5 # kg
substance <- 'PFOS'
C_water <- c(0) # ug/L
C_water_times <- c(0)
IV_intake_hourly <- c(0) # ug/h
IV_times <- c(0)
Oral_intake_hourly <- c(0.75*BW/24, 0) # mg/h
Oral_times <- c(0,6*30*24)
IV_doses <- c(0)*BW
IV_doses_times <- c(0)

Drinking_rate <- 0#BW*11/1000/24 # L/h
sample_time <- seq(0,18*30*24,1)

user_input <- list('BW'=BW, 'substance'=substance,
                   'C_water'=C_water,
                   'C_water_times'=C_water_times,
                   'IV_intake_hourly'=IV_intake_hourly,
                   'IV_times'=IV_times,
                   'Oral_intake_hourly'=Oral_intake_hourly,
                   'Oral_times'=Oral_times,
                   'IV_doses'=IV_doses, 
                   'IV_doses_times'=IV_doses_times,
                   'Drinking_rate'=Drinking_rate,
                   'sample_time'=sample_time)

fig6 = reproduce_results(user_input)
fig6
########################################################
BW = 3.5 # kg
substance <- 'PFOS'
C_water <- c(0) # mg/L
C_water_times <- c(0)
IV_intake_hourly <- c(0) # mg/h
IV_times <- c(0)
Oral_intake_hourly <- c(0) # mg/h
Oral_times <- c(0)
IV_doses <- c(2)*BW # mg 
IV_doses_times <- c(0)

Drinking_rate <- 0#BW*11/1000/24 # L/h
sample_time <- seq(0,160*24,0.1)

user_input <- list('BW'=BW, 'substance'=substance,
                   'C_water'=C_water,
                   'C_water_times'=C_water_times,
                   'IV_intake_hourly'=IV_intake_hourly,
                   'IV_times'=IV_times,
                   'Oral_intake_hourly'=Oral_intake_hourly,
                   'Oral_times'=Oral_times,
                   'IV_doses'=IV_doses, 
                   'IV_doses_times'=IV_doses_times,
                   'Drinking_rate'=Drinking_rate,
                   'sample_time'=sample_time)

fig5 = reproduce_results(user_input)
fig5

########################################################

# figure 3
BW = 3.5 # kg
substance <- 'PFOA'
C_water <- c(0) # ug/L
C_water_times <- c(0)
IV_intake_hourly <- c(0) # ug/h
IV_times <- c(0)
Oral_intake_hourly <- c(0) # mg/h
Oral_times <- c(0)
IV_doses <- c(10)*BW
IV_doses_times <- c(0)

Drinking_rate <- 0#BW*11/1000/24 # L/h
sample_time <- seq(0,100*24,1)

user_input <- list('BW'=BW, 'substance'=substance,
                   'C_water'=C_water,
                   'C_water_times'=C_water_times,
                   'IV_intake_hourly'=IV_intake_hourly,
                   'IV_times'=IV_times,
                   'Oral_intake_hourly'=Oral_intake_hourly,
                   'Oral_times'=Oral_times,
                   'IV_doses'=IV_doses, 
                   'IV_doses_times'=IV_doses_times,
                   'Drinking_rate'=Drinking_rate,
                   'sample_time'=sample_time)

fig3 = reproduce_results(user_input)
fig3
