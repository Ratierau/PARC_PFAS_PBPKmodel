##############################################################################################
#                                                                                            #
#  Function for summing log-normal distributions based on Fenton-Wilkinson approximation     #
#  All computations are based on https://arxiv.org/pdf/1502.03619.                           #
#                                                                                            #
#  Arguments:                                                                                #
#  mus   vector of distribution central values (in logarithmic form, ie. mus = log(gmeans))  #
#  sds   vector of distribution variances (in logarithmic form, ie. sds = log(gsds))         #
#                                                                                            #
##############################################################################################

fenton.wilkinson<-function(mus,sds) {
  
  # Central moments of the new distribution based on https://arxiv.org/pdf/1502.03619, formulas (6) and (7) if there is no mutual covariance
  #totalm<-sum(exp(mus)*exp((sds^2)/2))
  #totals<-sqrt(sum(exp(2*mus+(sds^2))*(exp(sds^2)-1)))

  # Note Jiri Kalina 3034-06-30: this approach is very prone to overflow if there is a distribution with extreme low gm^2/gsd ratio.
  # To handle all variety of distributions, I made this adjustment:
  
  prelim.m <-exp(mus)*exp((sds^2)/2)            # summands for totalm
  prelim.s2<-exp(2*mus+(sds^2))*(exp(sds^2)-1)  # summands for totals2
  
  ratio<-prelim.m/prelim.s2                     # ratio of m summands to s2 summands as a measure of extremness
  extremes<-which(ratio<1)                      # identifies extremely high variances
  
  sds2<--mus[extremes]/1.5
  #sds2[which(sds2<0)]<-0
  sds[extremes]<-sqrt(abs(sds2))                # fixing the variance at a level providing maximal estimate (tricky, more research needed))
  
  totalm<-sum(exp(mus)*exp((sds^2)/2))
  totals<-sqrt(sum(exp(2*mus+(sds^2))*(exp(sds^2)-1)))
  
  # Deriving characteristics of the resulting log-normal distribution based on central moments
  mu<-log(totalm/sqrt((totals/totalm)^2+1))
  sd<-sqrt(log(((totals/totalm)^2+1)))
  
  return(list("mu"=mu,"sd"=sd,"prelim.m"=prelim.m,"prelim.s2"=prelim.s2,"totalm"=totalm,"totals"=totals))
}




















