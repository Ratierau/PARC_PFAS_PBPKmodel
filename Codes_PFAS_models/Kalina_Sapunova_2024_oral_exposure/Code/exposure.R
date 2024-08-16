#########################################################################################################################
#                                                                                                                       #
#  Function for computing content of PFASs in food: exposure()                                                          #
#                                                                                                                       #
#  Arguments:                                                                                                           #
#  expoHierarchyCode.x   EFSA expo hierarchy code level > 2, if not provided, expoHierarchyCode.2 is used instead       #
#  expoHierarchyCode.2   EFSA expo hierarchy code level 2, used only if expoHierarchyCode.x is NA                       #
#  pfas                  compound code as in the data data.frame                                                        #
#  region                one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"              #
#  country               ISO 3166 alpha-2 code                                                                          #
#  year                  if given, uses trend (+ CI) for estimate exposure in that year                                 #
#  start                 lower limit of pfas data used for the estimation                                               #
#  end                   upper limit of pfas data used for the estimation                                               #
#  fit                   either "moments" for GM & GSD or "MLE" for most likelihood                                     #
#  method                either "numeric" for randomly generated distribution or "analytic" for returning GM and GSD    #
#  q                     quantiles to return if method is set to "numeric"                                              # 
#  data                  pfas food exposure data (data.frame as provided by Daria Sapunova)                             #
#  enable.regions        TRUE for replacing missing country data by the country region, FALSE otherwise                 #
#                                                                                                                       #
#########################################################################################################################

# This function computes exposure for individual food items in given country and given time (based on time trend) or over the whole period of selected data (if year=NA).
# The function uses MLE or simple 1st and 2nd moments estimation to define a lognormal distribution based on dataset of PFASs' food concentrations collected by Daria Sapunova (food_concentrations.rds).
# If there is not enough data for selected country and enable.regions=TRUE, it uses broader regions (4 parts of Europe according to HBM4EU) or the whole dataset instead of missing data.

exposure<-function(expoHierarchyCode.x=NA,expoHierarchyCode.2=substr(expoHierarchyCode.x,1,10),pfas,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),data=food_concentrations,enable.regions=TRUE) {
  
  warn<-character(0)
  
  #############################
  # Test of input consistency #
  #############################
  if(!is.na(expoHierarchyCode.x)) {
    if(expoHierarchyCode.2!=substr(expoHierarchyCode.x,1,10)) {
      stop(paste0("There is a mismatch between ",expoHierarchyCode.x," and ",expoHierarchyCode.2,"."))
    }
  }
  
  #####################
  # Selection of data #
  #####################  
  
  # First the function tries to find appropriate expoHierarchyCode.x records within given region,
  # if there are < 3 tries expoHierarchyCode.2 in given region and
  # if still there are < 3 tries expoHierarchyCode.2 in the whole dataset.
  if(is.na(country)) {
    selection<-data[which(data$expoHierarchyCode.x==expoHierarchyCode.x&data$chemical==pfas&data$Europe_region==region&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
    if(nrow(selection)<3) {
      if(!is.na(expoHierarchyCode.x)) {
        if(nchar(expoHierarchyCode.x)>nchar(expoHierarchyCode.2)) {
          warn<-c(warn,paste0("There is less than 3 records for",pfas," and ",expoHierarchyCode.x," in given period in ",region,". I am using ",expoHierarchyCode.2," instead."))
        }
      }
      selection<-data[which(data$expoHierarchyCode.2==expoHierarchyCode.2&data$chemical==pfas&data$Europe_region==region&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
      if(nrow(selection)<3) {
        selection<-data[which(data$expoHierarchyCode.2==expoHierarchyCode.2&data$chemical==pfas&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
        if(nrow(selection)<3) {
          warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in given period. No values to return."))
        } else {
          warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",region," and given period. I am using the whole dataset instead."))
        }
      }
    }
  } 
  
  # First the function tries to find appropriate expoHierarchyCode.x records within given country,
  # if there are < 3 tries expoHierarchyCode.2 in given country  and
  # if still there are < 3 tries expoHierarchyCode.2 in given region and
  # if still there are < 3 tries expoHierarchyCode.2 in the whole dataset.
  if(!is.na(country)) {
    selection<-data[which(data$expoHierarchyCode.x==expoHierarchyCode.x&data$chemical==pfas&grepl(country,data$country)&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
    if(nrow(selection)<3) {
      if(!is.na(expoHierarchyCode.x)) {
        if(nchar(expoHierarchyCode.x)>nchar(expoHierarchyCode.2)) {
          warn<-c(warn,paste0("There is less than 3 records for ",pfas," and ",expoHierarchyCode.x," in given period in ",country,". I am using ",expoHierarchyCode.2," instead."))
        }
      }
      selection<-data[which(data$expoHierarchyCode.2==expoHierarchyCode.2&data$chemical==pfas&grepl(country,data$country)&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
      if(nrow(selection)<3&enable.regions==TRUE) {
        region.of.country<-unique(data$Europe_region[which(data$country==country)])[1]
        selection<-data[which(data$expoHierarchyCode.2==expoHierarchyCode.2&data$chemical==pfas&data$Europe_region==region.of.country&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
        if(nrow(selection)<3) {
          selection<-data[which(data$expoHierarchyCode.2==expoHierarchyCode.2&data$chemical==pfas&!is.na(data$median)&data$date_midpoint>=start&data$date_midpoint<=end),]
          if(nrow(selection)<3) {
            warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in given period. No values to return."))
          } else {
            warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",country," neither the region of ",region.of.country," in given period. I am using the whole dataset instead."))
          }
        } else {
          warn<-c(warn,paste0("There is less then 3 records for ",pfas," and ",expoHierarchyCode.2," in ",country," in given period. I am using the region of ",region.of.country," instead."))
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
      
      if(method=="numeric") {
        results<-q*NA
        names(results)<-q
      }
      
      if(method=="analytic") {
        results<-c(NA,NA)
        names(results)<-c("gmean","gsd")
      }
    } else {
      y<-log10(subselection$median)
      x<-subselection$date_midpoint
      model<-lm(y~x)
      
      if(method=="numeric") {
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
      
      if(method=="analytic") {
        results<-c("gmean"=10^(model$coefficients[1]+model$coefficients[2]*year),"gsd"=10^((predict(model,newdata=data.frame("x"=year),level=0.682,interval="confidence")[1,"upr"]-predict(model,newdata=data.frame("x"=year),level=0.682,interval="confidence")[1,"lwr"])/2))
        names(results)<-c("gmean","gsd")
      }
    }
  } else {
    
    ##########
    # Static #
    ########## 
    
    subselection<-selection[which(!is.na(selection$median)),]
    if(nrow(subselection)<3) {
      warn<-c(warn,paste0("There is less than 3 points. It is not possible to derive any distribution."))
      
      if(method=="numeric") {
        results<-q*NA
        names(results)<-q
      }
      
      if(method=="analytic") {
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
        if(length(unique(subselection$median))>1) {
          fitted.distribution<-fitdist(subselection$median,"lnorm")
          mean<-fitted.distribution$estimate["meanlog"]/log(10)
          sd  <-fitted.distribution$estimate["sdlog"]  /log(10)
        } else {
          warn<-c(warn,paste0("There is not enough values for MLE, using moments instead."))
          mean<-mean(log10(subselection$median))
          sd  <-  sd(log10(subselection$median))
        }
      }
      
      
      ## Generating results
      
      if(method=="numeric") {
        results<-numeric(0)*NA
        for(quantile in q) {
          results<-c(results,10^qnorm(quantile,mean,sd))
        }
        names(results)<-q
      }
      
      if(method=="analytic") {
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