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