#############################################################################################################################################################
#                                                                                                                                                           #
#  Function for individual exposure based on food consumption in mass/day: exposure.mass()                                                                  #
#                                                                                                                                                           #
#  Arguments:                                                                                                                                               #
#  food_mass_individual   data.frame with individual daily food consumption in grams (columns "expoHierarchyCode.x", "expoHierarchyCode.2", "gmean","gsd)   #
#  pfass                  vector of pfas compound codes as in the data data.frame                                                                           #
#  region                 one of "Northern Europe", "Western Europe", "Eastern Europe" or "Southern Europe"                                                 #
#  country                ISO 3166 alpha-2 code                                                                                                             #
#  year                   if given, uses trend (+ CI) for estimate exposure in that year                                                                    #
#  start                  lower limit of pfas data used for the estimation                                                                                  #
#  end                    upper limit of pfas data used for the estimation                                                                                  #
#  fit                    either "moments" for GM & GSD or "MLE" for most likelihood                                                                        #
#  method                 either "numeric" for summing randomly generated distributions or "analytic" for using Fenton-Wilkinson approximation              #
#  n                      number of generated values in case of method="numeric"                                                                            #
#  data                   pfas food exposure data (data.frame as provided by Daria Sapunova)                                                                #
#  enable.regions         TRUE for replacing missing country data by the country region, FALSE otherwise                                                    #
#                                                                                                                                                           #
#############################################################################################################################################################

# This function uses previously defined function exposure() to compute PFASs' exposure for one individual with known food consumption in mass/day.
# It uses exposure() for handling food PFASs' concentrations and individual food consumptions on form of a data frame

exposure.mass<-function(food_mass_individual,pfass,region,country=NA,year=NA,start=-Inf,end=Inf,fit="moments",method="numeric",q=c(0.05,0.25,0.50,0.75,0.95),n=10000,data=food_concentrations,enable.regions=TRUE) {
  
  warn<-character(0)
  results<-data.frame()
  
  # Using Fenton-Wilkinson approximation for summing duplicate records
  fmnrows<-nrow(food_mass_individual)
  food_mass_individual_aggr<-food_mass_individual[0,]
  unique<-unique(food_mass_individual[,c("expoHierarchyCode.x","expoHierarchyCode.2")])
  for(i in 1:nrow(unique)) {
    mus<-log(food_mass_individual[which(food_mass_individual$expoHierarchyCode.x==unique$expoHierarchyCode.x[i]&food_mass_individual$expoHierarchyCode.2==unique$expoHierarchyCode.2[i]),"gmean"])
    sds<-log(food_mass_individual[which(food_mass_individual$expoHierarchyCode.x==unique$expoHierarchyCode.x[i]&food_mass_individual$expoHierarchyCode.2==unique$expoHierarchyCode.2[i]),"gsd"])

    food_mass_individual_aggr<-rbind(food_mass_individual_aggr,data.frame("expoHierarchyCode.x"=unique$expoHierarchyCode.x[i],
                                                                          "expoHierarchyCode.2"=unique$expoHierarchyCode.2[i],
                                                                          "gmean"=exp(fenton.wilkinson(mus,sds)$mu),
                                                                          "gsd"=exp(fenton.wilkinson(mus,sds)$sd)))
  }
  food_mass_individual<-food_mass_individual_aggr
  
  if(fmnrows!=nrow(food_mass_individual)) {
    warn<-c(warn,"In food_mass_individual, multiple rows for same food items were aggregated using Fenton-Wilkinson approximation.")
  }
  
  # Calculation itself
  for(pfas in pfass) {
    
    contrib<-data.frame("gmean.food"=numeric(0),"gsd.food"=numeric(0),"gmean.mass"=numeric(0),"gsd.mass"=numeric(0))
    for(i in 1:nrow(food_mass_individual)) {
      contrib<-rbind(contrib,c(exposure(food_mass_individual$expoHierarchyCode.x[i],food_mass_individual$expoHierarchyCode.2[i],pfas,region=region,country=country,year=year,start=start,end=end,fit=fit,method="analytic",data=data,enable.regions=enable.regions),"gmean.food"=food_mass_individual$gmean[i],"gsd.food"=food_mass_individual$gsd[i]))
    }
    rownames(contrib)<-food_mass_individual$expoHierarchyCode.x
    colnames(contrib)<-c("gmean.food","gsd.food","gmean.mass","gsd.mass")
    
    if(method=="numeric") {
      ## Monte-Carlo
      cs<-numeric(0)
      for(l in 1:n) {
        c<-0
        for(row in 1:nrow(contrib)) {
          if(!is.na(rowSums(contrib)[row])) {
            
            c<-c+rlnorm(1,log(contrib$gmean.food[row]),log(contrib$gsd.food[row]))*rlnorm(1,log(contrib$gmean.mass[row]),log(contrib$gsd.mass[row]))
            
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
      mus<-log(contrib$gmean.food*contrib$gmean.mass)[which(!is.na(rowSums(contrib)))]
      sds<-sqrt(log(contrib$gsd.food)^2+log(contrib$gsd.mass)^2)[which(!is.na(rowSums(contrib)))]
      
      # Deriving characteristics of the resulting log-normal distribution based on central moments
      mu<-fenton.wilkinson(mus,sds)$mu
      sd<-fenton.wilkinson(mus,sds)$sd
      
      # Storing as results
      addline<-data.frame("gmean"=exp(mu),"gsd"=exp(sd))
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