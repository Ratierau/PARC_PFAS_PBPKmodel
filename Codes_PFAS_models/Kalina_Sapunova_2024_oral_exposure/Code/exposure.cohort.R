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