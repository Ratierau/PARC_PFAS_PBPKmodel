library("deSolve")
library("tidyverse")
library("xlsx")
library("truncnorm")
library("snow")
library("doSNOW")
library("parallel")

################################################################
#                                                              #
#  Function trunclnorm() for truncated lognormal distribution  #
#                                                              #
################################################################

rtrunclnorm<-function(x,a=0,b=Inf,meanlog=0,sdlog=1) {
  result<-exp(rtruncnorm(x,log(a),log(b),meanlog,sdlog))
  return(result)
}

###########################################
#                                         #
#  Reverse modelling - scenario for PFOS  #
#                                         #
###########################################

popul<-read.xlsx("population.xlsx",1)          # reads the file with individuals' characteristics (BW, height, sex, age) and blood concentration
param<-read.xlsx("scenario_iv_pfos.xlsx",2)    # parametrization file derived from the MCRA format
results<-data.frame("individual"=character(0),
                    "conc_blood"=numeric(0),
                    "edi"=numeric(0))
loops<-100

# PFOS 
compound<-"PFOS"
sel<-param[which((param$Compound==compound|is.na(param$Compound))&!is.na(param$Param1)),]

#####################
#                   #
#  Parallelization  #
#                   #
#####################

workers<-makeCluster(4,type="SOCK")
registerDoSNOW(workers)

results<-foreach(person=1:4,.verbose=TRUE,.packages=c("deSolve","tidyverse","dplyr","truncnorm")) %dopar% {
  source("pfas_pbk_lifetime_model.R")
  source("faster_convergence_full_param.R")
  
  edi.person<-numeric(0)
  for(loop in 1:loops) {
    
    factors<-character(0)
    parames<-numeric(0)
    
    for(i in 1:nrow(sel)) {
      factors<-c(factors,sel$Parameter[i])
      if(sel$DistributionType[i]=="const") {
        parames<-c(parames,sel$Param1[i])
      } 
      if(sel$DistributionType[i]=="unif") {
        parames<-c(parames,runif(1,sel$Param1[i],sel$Param2[i]))
      }
      if(sel$DistributionType[i]=="norm") {
        parames<-c(parames,rnorm(1,sel$Param1[i],sel$Param2[i]))
      }
      if(sel$DistributionType[i]=="lnorm") {
        parames<-c(parames,rlnorm(1,log(sel$Param1[i]),sel$Param2[i])) # First argument is exp(meanlog), second is sdlog! 
      }
      if(sel$DistributionType[i]=="truncnorm") {
        parames<-c(parames,rtruncnorm(1,sel$Param3[i],sel$Param4[i],sel$Param1[i],sel$Param2[i]))
      }
      if(sel$DistributionType[i]=="trunclnorm") {
        parames<-c(parames,rtrunclnorm(1,sel$Param3[i],sel$Param4[i],log(sel$Param1[i]),sel$Param2[i])) # Non-systematic (according to Klára): first argument is exp(meanlog), second is sdlog! 
      }
    }
    
    # removing redundant factors
    redundant<-c("VBloodC","VSkC","hcell")
    removethese<-which(is.element(factors,redundant))
    if(length(removethese)>0) {
      factors<-factors[-removethese]
      parames<-parames[-removethese]
    }
    
    # adding concentration, age and sex
    c<-popul$conc_blood[person]
    age<-popul$age[person]
    sex<-popul$sex[person]
    weight<-popul$weight[person]
    height<-popul$height[person]/100 # height in m (recalculated from cm)
    factors<-c("c","age","Gender","BW","height",factors)
    parames<-c(c,age,paste0("'",sex,"'"),weight,height,parames)
    names(parames)<-factors
    
    # modelling
    to.evaluate<-paste0("faster.convergence.full.param(",paste0("'",names(parames),"'=",parames,collapse=","),",FUN='joost.model',max.loops=100,max.diff=0.01)")
    edi.person<-c(edi.person,eval(parse(text=to.evaluate)))
    
  }
  rbind(data.frame("individual"=popul$individual[person],
                   "conc_blood"=popul$conc_blood[person],
                   "edi"=edi.person))
}
stopCluster(workers)

results.a<-do.call("rbind", results)

edi.persons<-aggregate(list("edi"=results.a$edi),by=list("individual"=results.a$individual),FUN="median")$edi

##################
#                #
#  Save results  #
#                #
##################

popul$edi<-NA
popul$edi<-edi.persons
write.xlsx(popul,file="results_pfos.xlsx")
write.xlsx(results.a,file="results_pfos_individuals.xlsx",row.names=FALSE)