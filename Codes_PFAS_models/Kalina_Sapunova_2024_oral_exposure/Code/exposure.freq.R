###################################################################################################################################################################################################
#                                                                                                                                                                                                 #
#  Function for individual exposure including individual food frequencies: exposure.freq()                                                                                                        #
#                                                                                                                                                                                                 #
#  Arguments:                                                                                                                                                                                     #
#  food.freq.individual   data.frame with individual weekly food frequencies (columns "expoHierarchyCode.x", "expoHierarchyCode.2", "freq")                                                       #
#  portions               data.frame with portion sizes in grams in form of log-normal distributions (columns "expoHierarchyCode.x", "expoHierarchyCode.2", "gmean", "gsd", "region", "country")  #
#  pfass                  vector of pfas compound codes as in the data data.frame                                                                                                                 #
#  region                 one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"                                                                                       #
#  country                ISO 3166 alpha-2 code                                                                                                                                                   #
#  year                   if given, uses trend (+ CI) for estimate exposure in that year                                                                                                          #
#  start                  lower limit of pfas data used for the estimation                                                                                                                        #
#  end                    upper limit of pfas data used for the estimation                                                                                                                        #
#  fit                    either "moments" for GM & GSD or "MLE" for most likelihood                                                                                                              #
#  method                 either "numeric" for summing randomly generated distributions or "analytic" for using Fenton-Wilkinson approximation                                                    #
#  n                      number of generated values in case of method="numeric"                                                                                                                  #
#  data                   pfas food exposure data (data.frame as provided by Daria Sapunova)                                                                                                      #
#  enable.regions         TRUE for replacing missing country data by the country region, FALSE otherwise                                                                                          #
#                                                                                                                                                                                                 #
###################################################################################################################################################################################################

# This function uses previously defined function exposure() to compute PFASs' exposure for one individual with known food frequencies.
# It uses exposure() for handling food PFASs' concentrations, portions.rds for a distribution of a typical portion size and individual food frequencies on form of a data frame

exposure.freq<-function(food.freq.individual,portions,pfass,region=NA,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),n=10000,data=data_pfas,enable.regions=TRUE) {
  
  # Selects only portions for specified region and country.
  if(is.na(country)) {
    if(is.na(region)) {
      # aggregating portions over the whole dataset
      portions<-data.frame(      aggregate(list("gmean"=portions$gmean),by=list("expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE),
                           "gsd"=aggregate(list("gsd"  =portions$gsd  ),by=list("expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE)$gsd)
      rownames(portions)<-portions$expoHierarchyCode.x
    }
    else {
      # aggregating portions according to region
      portions<-portions[which(portions$region==region),]     
      portions<-data.frame(      aggregate(list("gmean"=portions$gmean),by=list("region"=portions$region,"expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE),
                           "gsd"=aggregate(list("gsd"  =portions$gsd  ),by=list("region"=portions$region,"expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE)$gsd)
      rownames(portions)<-portions$expoHierarchyCode.x
    }
  } else {
    # aggregating portions according to country
    portions<-portions[which(portions$country==country),]
    portions<-data.frame(      aggregate(list("gmean"=portions$gmean),by=list("country"=portions$country,"expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE),
                         "gsd"=aggregate(list("gsd"  =portions$gsd  ),by=list("country"=portions$country,"expoHierarchyCode.x"=portions$expoHierarchyCode.x,"expoHierarchyCode.2"=portions$expoHierarchyCode.2),FUN="median",na.rm=TRUE)$gsd)
    rownames(portions)<-portions$expoHierarchyCode.x
  }
  
  # prepare empty output variables
  warn<-character(0)
  results<-data.frame()
  
  # food frequency number of lines
  ffnl<-nrow(food.freq.individual)
  food.freq.individual<-aggregate(list("freq"=food.freq.individual$freq),by=list("expoHierarchyCode.x"=food.freq.individual$expoHierarchyCode.x,"expoHierarchyCode.2"=food.freq.individual$expoHierarchyCode.2),FUN="sum")
  if(ffnl!=nrow(food.freq.individual)) {
    warn<-c(warn,"In food.freq.individual, multiple rows for same food items were aggregated.")
  }
  
  for(pfas in pfass) {
    
    contrib<-data.frame("gmean"=numeric(0),"gsd"=numeric(0),weight=numeric(0))
    for(i in 1:nrow(food.freq.individual)) {
      contrib<-rbind(contrib,c(exposure(food.freq.individual$expoHierarchyCode.x[i],food.freq.individual$expoHierarchyCode.2[i],pfas,region=region,country=country,year=year,fit=fit,method="analytic",data=data,enable.regions=enable.regions),"weight"=food.freq.individual$freq[i]))
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
    
    if(method=="analytic") {
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