#########################################################################################################################
#                                                                                                                       #
#  Function for computing content of PFASs in food: exposure()                                                          #
#                                                                                                                       #
#  Arguments:                                                                                                           #
#  expoHierarchyCode.x   not used in this version                                                                       #
#  expoHierarchyCode.2   expo hierarchy code EFSA                                                                       #
#  pfas                  compound code as in the data data.frame                                                        #
#  region                one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"              #
#  country               ISO 3166 alpha-2 code                                                                          #
#  year                  if given, uses trend (+ CI) for estimate exposure in that year                                 #
#  start                 lower limit of pfas data used for the estimation                                               #
#  end                   upper limit of pfas data used for the estimation                                               #
#  fit                   either "moments" for GM & GSD or "MLE" for most likelihood                                     #
#  method                either "numeric" for randomly generated distribution or "analytical" for returning GM and GSD  #
#  q                     quantiles to return if method is set to "numeric"                                              # 
#  data                  pfas food exposure data (data.frame as provided by Daria Sapunova)                             #
#  enable.regions        TRUE for replacing missing country data by the country region, FALSE otherwise                 #
#                                                                                                                       #
#########################################################################################################################

# This function computes exposure for individual food items in given country and given time (based on time trend) or over the whole period of selected data (if year=NA).
# The function uses MLE or simple 1st and 2nd moments estimation to define a lognormal distribution based on dataset of PFASs' food concentrations collected by Daria Sapunova (data_pfas.rds).
# If there is not enough data for selected country and enable.regions=TRUE, it uses broader regions (4 parts of Europe according to HBM4EU) or the whole dataset instead of missing data.

exposure<-function(expoHierarchyCode.x,expoHierarchyCode.2=substr(expoHierarchyCode.x,1,10),pfas,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),data=data_pfas,enable.regions=TRUE) {
  
  warn<-character(0)
  
  #####################
  # Selection of data #
  #####################  
  
  # At the time the data_pfas dataset only contains food expo 2 code, so the food expo x code is not used.
  if(is.na(country)) {
    selection<-data[which(data$food_cat_code==expoHierarchyCode.2&data$chemical==pfas&data$Europe_region==region&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
    if(nrow(selection)<3) {
      selection<-data[which(data$food_cat_code==expoHierarchyCode.2&data$chemical==pfas&!is.na(data$median)),]
      if(nrow(selection)<3) {
        warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in given period. No values to return."))
      } else {
        warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",region,". I am using the whole dataset instead."))
      }
    }
  } else {
    selection<-data[which(data$food_cat_code==expoHierarchyCode.2&data$chemical==pfas&grepl(country,data$country)&!is.na(data$median)),]
    
    if(nrow(selection)<3&enable.regions==TRUE) {
      region.of.country<-unique(data_pfas$Europe_region[which(data_pfas$country==country)])[1]
      selection<-data[which(data$food_cat_code==expoHierarchyCode.2&data$chemical==pfas&data$Europe_region==region.of.country&!is.na(data$median)),]
      if(nrow(selection)<3) {
        selection<-data[which(data$food_cat_code==expoHierarchyCode.2&data$chemical==pfas&!is.na(data$median)),]
        if(nrow(selection)<3) {
          warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in given period. No values to return."))
        } else {
          warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",country," neither the region of ",region.of.country,". I am using the whole dataset instead."))
        }
      } else {
        warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",country,". I am using the region of ",region.of.country," instead."))
      }
    }
    
    if(nrow(selection)<3&enable.regions==FALSE) {
      if(!is.na(country)) {
        warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",country,". No values to return."))
      } else {
        warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",region,". No values to return."))
      }
    }
  }
  
  ##########################
  # Computing distribution #
  ##########################
  
  if(!is.na(year)) {
    
    ################
    # Time related #
    ################
    
    subselection<-selection[which(!is.na(selection$median&!is.na(selection$date_midpoint))),]
    if(nrow(subselection)<3|length(unique(subselection$date_midpoint))<2) {
      warn<-c(warn,paste0("There is less than 3 timepoints or they are all in the same year. It is not possible to derive any trend."))
    } else {
      y<-log10(subselection$median)
      x<-subselection$date_midpoint
      model<-lm(y~x)
      
      results<-numeric(0)*NA
      for(quantile in q) {
        if(quantile<0.5) {
          results<-c(results,10^predict(model,newdata=data.frame("x"=year),level=+1-2*quantile,interval="confidence")[1,"lwr"])
        }
        if(quantile==0.5) {
          results<-c(results,10^predict(model,newdata=data.frame("x"=year),level=0,interval="confidence")[1,"fit"])
        }
        if(quantile>0.5) {
          results<-c(results,10^predict(model,newdata=data.frame("x"=year),level=-1+2*quantile,interval="confidence")[1,"upr"])
        }
      }
      names(results)<-q
    }
  } else {
    
    ##########
    # Static #
    ########## 
    
    subselection<-selection[which(!is.na(selection$median)),]
    if(nrow(subselection)<3) {
      warn<-c(warn,paste0("There is less than 3 points. It is not possible to derive any distribution."))
      
      if(method=="numeric") {
        results<-numeric(0)*NA
        names(results)<-q
      }
      
      if(method=="analytical") {
        results<-c(NA,NA)
        names(results)<-c("gmean","gsd")
      }
      
    } else {
      
      ## Computing distribution moments
      
      if(fit=="moments") {
        mean<-mean(log10(subselection$median))
        sd  <-  sd(log10(subselection$median))
      }
      
      if(fit=="mle") {
        fitted.distribution<-fitdist(subselection$median,"lnorm")
        mean<-fitted.distribution$estimate["meanlog"]/log(10)
        sd  <-fitted.distribution$estimate["sdlog"]  /log(10)
      }
      
      
      ## Generating results
      
      if(method=="numeric") {
        results<-numeric(0)*NA
        for(quantile in q) {
          results<-c(results,10^qnorm(quantile,mean,sd))
        }
        names(results)<-q
      }
      
      if(method=="analytical") {
        results<-c(10^mean,10^sd)
        names(results)<-c("gmean","gsd")
      }
    }
  }
  
  ############
  # Warnings #
  ############ 
  
  for(i in warn) {
    warning(i,call.=FALSE)
  }
  return(results)
}



################################################################################################################################################################
#                                                                                                                                                              #
#  Function for individual exposure including individual food frequencies: exposure.freq()                                                                     #
#                                                                                                                                                              #
#  Arguments:                                                                                                                                                  #
#  food.freq.individual   data.frame with individual weekly food frequencies (columns "expoHierarchyCode.x", "expoHierarchyCode.2", "freq")                    #
#  portions               data.frame with portion sizes in grams in form of log-normal distributions (columns "expo.level.x", "expo.level.2", "gmean", "gsd")  #
#  pfass                  vector of pfas compound codes as in the data data.frame                                                                              #
#  region                 one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"                                                    #
#  country                ISO 3166 alpha-2 code                                                                                                                #
#  year                   if given, uses trend (+ CI) for estimate exposure in that year                                                                       #
#  start                  lower limit of pfas data used for the estimation                                                                                     #
#  end                    upper limit of pfas data used for the estimation                                                                                     #
#  fit                    either "moments" for GM & GSD or "MLE" for most likelihood                                                                           #
#  method                 either "numeric" for summing randomly generated distributions or "analytical" for using Fenton-Wilkinson approximation               #
#  n                      number of generated values in case of method="numeric"                                                                               #
#  data                   pfas food exposure data (data.frame as provided by Daria Sapunova)                                                                   #
#  enable.regions         TRUE for replacing missing country data by the country region, FALSE otherwise                                                       #
#                                                                                                                                                              #
################################################################################################################################################################

# This function uses previously defined function exposure() to compute PFASs' exposure for one individual with known food frequencies.
# It uses exposure() for handling food PFASs' concentrations, portions.rds for a distribution of a typical portion size and individual food frequencies on form of a data frame

exposure.freq<-function(food.freq.individual,portions,pfass,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),n=10000,data=data_pfas,enable.regions=TRUE) {
  
  warn<-character(0)
  results<-data.frame()
  
  ffnl<-nrow(food.freq.individual)
  food.freq.individual<-aggregate(list("freq"=food.freq.individual$freq),by=list("expoHierarchyCode.x"=food.freq.individual$expoHierarchyCode.x,"expoHierarchyCode.2"=food.freq.individual$expoHierarchyCode.2),FUN="sum")
  if(ffnl!=nrow(food.freq.individual)) {
    warn<-c(warn,"In food.freq.individual, multiple rows for same food items were aggregated.")
  }
  
  for(pfas in pfass) {
    
    contrib<-data.frame("gmean"=numeric(0),"gsd"=numeric(0),weight=numeric(0))
    for(i in 1:nrow(food.freq.individual)) {
      contrib<-rbind(contrib,c(exposure(food.freq.individual$expoHierarchyCode.x[i],food.freq.individual$expoHierarchyCode.2[i],pfas,region=region,country=country,year=year,fit=fit,method="analytical",data=data,enable.regions=enable.regions),"weight"=food.freq.individual$freq[i]))
    }
    rownames(contrib)<-food.freq.individual$expoHierarchyCode.x
    colnames(contrib)<-c("gmean","gsd","weight")
    
    if(method=="numeric") {
      ## Monte-Carlo
      cs<-numeric(0)
      for(l in 1:n) {
        c<-0
        for(row in 1:nrow(contrib)) {
          if(!is.na(rowSums(contrib)[row])) {
            
            # finds appropriate portion size distribution
            portion.dist<-portions[which(portions$expo.level.x==rownames(contrib)[row]),]
            
            if(nrow(portion.dist)>0) {
              c<-c+contrib$weight[row]*rlnorm(1,log(contrib$gmean[row]),log(contrib$gsd[row]))*rlnorm(1,log(portion.dist$gmean[1]),log(portion.dist$gsd[1]))
            } else {
              warn<-c(warn,paste0("No portion sizes were found for the food item ",rownames(contrib)[row],"."))
            }
          } else  {
            warn<-c(warn,paste0("No concentration of ",pfas," was found for the food item ",rownames(contrib)[row],"."))
          }
        }
        cs<-c(cs,c)
      }
      addline<-as.data.frame(t(quantile(cs,q)))
      colnames(addline)<-q
      results<-rbind(results,addline)
    }
    
    if(method=="analytical") {
      concentrationmis<-log(contrib$gmean*contrib$weight)[which(!is.na(rowSums(contrib)))]
      concentrationsds<-log(contrib$gsd  *             1)[which(!is.na(rowSums(contrib)))]
      
      portionmis<-log(portions[rownames(contrib)[which(!is.na(rowSums(contrib)))],"gmean"])
      portionsds<-log(portions[rownames(contrib)[which(!is.na(rowSums(contrib)))],"gsd"]  )
      
      # Putting together lognormally distributed concentration and portion size.
      mis<-concentrationmis+portionmis
      sds<-sqrt(concentrationsds^2+portionsds^2)
      
      # Central moments of the new distribution based on https://arxiv.org/pdf/1502.03619, formulas (6) and (7) if there is no mutual covariance
      totalm<-sum(exp(mis)*exp((sds^2)/2),na.rm=TRUE)
      totals<-sqrt(sum(exp(2*mis+(sds^2))*(exp(sds^2)-1),na.rm=TRUE))
      
      # Deriving characteristics of the resulting log-normal distribution based on central moments
      mi<-log(totalm/sqrt((totals/totalm)^2+1))
      sd<-sqrt(log(((totals/totalm)^2+1)))
      
      # Storing as results
      addline<-data.frame("gmean"=exp(mi),"gsd"=exp(sd))
      results<-rbind(results,addline)
    }
  }
  rownames(results)<-pfass
  
  ############
  # Warnings #
  ############ 
  
  for(i in warn) {
    warning(i,call.=FALSE)
  }
  return(results)
}



################################################################################################################################################################
#                                                                                                                                                              #
#  Function for individual exposure based on food consumption in mass/day: exposure.mass()                                                                     #
#                                                                                                                                                              #
#  Arguments:                                                                                                                                                  #
#  food.mass.individual   data.frame with individual daily food consumption in grams (columns "expoHierarchyCode.x", "expoHierarchyCode.2", "mass")            #
#  pfass                  vector of pfas compound codes as in the data data.frame                                                                              #
#  region                 one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"                                                    #
#  country                ISO 3166 alpha-2 code                                                                                                                #
#  year                   if given, uses trend (+ CI) for estimate exposure in that year                                                                       #
#  start                  lower limit of pfas data used for the estimation                                                                                     #
#  end                    upper limit of pfas data used for the estimation                                                                                     #
#  fit                    either "moments" for GM & GSD or "MLE" for most likelihood                                                                           #
#  method                 either "numeric" for summing randomly generated distributions or "analytical" for using Fenton-Wilkinson approximation               #
#  n                      number of generated values in case of method="numeric"                                                                               #
#  data                   pfas food exposure data (data.frame as provided by Daria Sapunova)                                                                   #
#  enable.regions         TRUE for replacing missing country data by the country region, FALSE otherwise                                                       #
#                                                                                                                                                              #
################################################################################################################################################################

# This function uses previously defined function exposure() to compute PFASs' exposure for one individual with known food consumption in mass/day.
# It uses exposure() for handling food PFASs' concentrations and individual food consumptions on form of a data frame

exposure.mass<-function(food.mass.individual,pfass,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),n=10000,data=data_pfas,enable.regions=TRUE) {
  
  warn<-character(0)
  results<-data.frame()
  
  fmnl<-nrow(food.mass.individual)
  food.mass.individual<-aggregate(list("mass"=food.mass.individual$mass),by=list("expoHierarchyCode.x"=food.mass.individual$expoHierarchyCode.x,"expoHierarchyCode.2"=food.mass.individual$expoHierarchyCode.2),FUN="sum")
  if(fmnl!=nrow(food.mass.individual)) {
    warn<-c(warn,"In food.mass.individual, multiple rows for same food items were aggregated.")
  }
  
  for(pfas in pfass) {
    
    contrib<-data.frame("gmean"=numeric(0),"gsd"=numeric(0),weight=numeric(0))
    for(i in 1:nrow(food.mass.individual)) {
      contrib<-rbind(contrib,c(exposure(food.mass.individual$expoHierarchyCode.x[i],food.mass.individual$expoHierarchyCode.2[i],pfas,region=region,country=country,year=year,start=start,end=end,fit=fit,method="analytical",data=data,enable.regions=enable.regions),"weight"=food.mass.individual$mass[i]))
    }
    rownames(contrib)<-food.mass.individual$expoHierarchyCode.x
    colnames(contrib)<-c("gmean","gsd","weight")
    
    if(method=="numeric") {
      ## Monte-Carlo
      cs<-numeric(0)
      for(l in 1:n) {
        c<-0
        for(row in 1:nrow(contrib)) {
          if(!is.na(rowSums(contrib)[row])) {
            
            c<-c+contrib$weight[row]*rlnorm(1,log(contrib$gmean[row]),log(contrib$gsd[row]))
            
          } else  {
            warn<-c(warn,paste0("No concentration of ",pfas," was found for the food item ",rownames(contrib)[row],"."))
          }
        }
        cs<-c(cs,c)
      }
      addline<-as.data.frame(t(quantile(cs,q)))
      colnames(addline)<-q
      results<-rbind(results,addline)
    }
    
    if(method=="analytical") {
      concentrationmis<-log(contrib$gmean*contrib$weight)[which(!is.na(rowSums(contrib)))]
      concentrationsds<-log(contrib$gsd  *             1)[which(!is.na(rowSums(contrib)))]
      
      # Renaming to keep consistency with combined() function.
      mis<-concentrationmis
      sds<-concentrationsds
      
      # Central moments of the new distribution based on https://arxiv.org/pdf/1502.03619, formulas (6) and (7) if there is no mutual covariance
      totalm<-sum(exp(mis)*exp((sds^2)/2),na.rm=TRUE)
      totals<-sqrt(sum(exp(2*mis+(sds^2))*(exp(sds^2)-1),na.rm=TRUE))
      
      # Deriving characteristics of the resulting log-normal distribution based on central moments
      mi<-log(totalm/sqrt((totals/totalm)^2+1))
      sd<-sqrt(log(((totals/totalm)^2+1)))
      
      # Storing as results
      addline<-data.frame("gmean"=exp(mi),"gsd"=exp(sd))
      results<-rbind(results,addline)
    }
  }
  rownames(results)<-pfass
  
  ############
  # Warnings #
  ############ 
  
  for(i in warn) {
    warning(i,call.=FALSE)
  }
  return(results)
}



#########################################################################################################################################################################
#                                                                                                                                                                       #
#  Function for the whole cohort: exposure.cohort()                                                                                                                     #
#                                                                                                                                                                       #
#  Arguments:                                                                                                                                                           #
#  food.freq        data.frame with weekly food frequencies (columns "expoHierarchyCode.x", "expoHierarchyCode.2", and one column with participant id for each person)  #
#  portions         data.frame with portion sizes in grams in form of log-normal distributions (columns "expo.level.x", "expo.level.2", "gmean", "gsd")                 #
#  pfass            vector of pfas compound codes as in the data data.frame                                                                                             #
#  region           one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"                                                                   #
#  country          ISO 3166 alpha-2 code                                                                                                                               #
#  year             if given, uses trend (+ CI) for estimate exposure in that year                                                                                      #
#  start            lower limit of pfas data used for the estimation                                                                                                    #
#  end              upper limit of pfas data used for the estimation                                                                                                    #
#  fit              either "moments" for GM & GSD or "MLE" for most likelihood                                                                                          #
#  method           either "numeric" for summing randomly generated distributions or "analytical" for using Fenton-Wilkinson approximation                              #
#  n                number of generated values in case of method="numeric"                                                                                              #
#  data             pfas food exposure data (data.frame as provided by Daria Sapunova)                                                                                  #
#  enable.regions   TRUE for replacing missing country data by the country region, FALSE otherwise                                                                      #
#                                                                                                                                                                       #
#########################################################################################################################################################################

exposure.cohort<-function(food.freq,portions,pfass,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),n=10000,data=data_pfas,enable.regions=TRUE) {
  
  if(method=="numeric") {
    results<-as.data.frame(matrix(NA,0,length(q)+2))
    colnames(results)<-c("id","pfas",paste0("q_",as.character(q)))
  }
  
  if(method=="analytical") {
    results<-as.data.frame(matrix(NA,0,4))
    colnames(results)<-c("id","pfas","gmean","gsd")
  }
  
  for(id in colnames(food.freq)[-c(1,2)]) {
    food.freq.individual<-food.freq[,c("expoHierarchyCode.x","expoHierarchyCode.2",as.character(id))]
    colnames(food.freq.individual)[ncol(food.freq.individual)]<-"freq"
    
    e<-exposure.freq(food.freq.individual=food.freq.individual,portions=portions,pfass=pfass,country=country,year=year,fit=fit,method=method,q=q,n=n,data=data,enable.regions=enable.regions)
    e$id<-id
    e$pfas=rownames(e)
    
    results<-rbind(results,e)
  }
  
  results<-cbind(results[,c("id","pfas")],results[,-c(which(colnames(results)%in%c("id","pfas")))])
  rownames(results)<-NULL
  
  return(results)
}

