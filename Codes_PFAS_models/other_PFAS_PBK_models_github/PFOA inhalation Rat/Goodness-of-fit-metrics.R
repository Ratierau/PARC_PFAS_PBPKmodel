
# SODI function the returns the SODI index described in Tsiros et al.2024

# predictions: list of vectors containing the predicted data

# names of the compartments

SODI <- function(observations, predictions, comp.names =NULL){
  
  # Check if the user provided the correct input format
  
  if (!is.list(observations) || !is.list(predictions)){
    
    stop(" The observations and predictions must be lists")
    
  }
  
  # Check if the user provided equal length lists
  
  if (length(observations) != length(predictions)){
    
    stop(" The observations and predictions must have the same compartments")
    
  }
  
  Ncomp <- length(observations) # Number of compartments
  
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  
  #loop over the compartments
  
  for (i in 1:Ncomp){
    
    Et <- 0 #relative error with observations
    
    St <- 0  #relative error with simulations
    
    N <- length(observations[[i]]) # number of observations for compartment i
    
    # Check if observations and predictions have equal length
    
    if(N != length(predictions[[i]])){
      
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
      
    }
    
    N_obs[i] <- N # populate the N_obs vector
    
    for (j in 1:N){
      
      # sum of relative squared errors (error = observations - predictions)
      
      Et <- Et + ( abs(observations[[i]][j] - predictions[[i]][j])  / observations[[i]][j] )  ^2
      
      St <- St + ( abs(observations[[i]][j] - predictions[[i]][j])  / predictions[[i]][j] )  ^2
      
    }
    
    
    
    # root mean of the square of observed values
    
    RMEt <- sqrt(Et/N)
    
    # root mean of the square of simulated values
    
    RMSt <- sqrt( St/N)
    
    
    
    I[i] <- (RMEt + RMSt)/2   
    
  }
  
  # Total number of observations
  
  Ntot <- sum(N_obs)
  
  # Initialise the consolidated discrepancy index
  
  Ic <-0
  
  for (i in 1:Ncomp){
    
    # Give weight to compartments with more observations (more information)
    
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
    
  }
  
  # Name the list of compartment discrepancy indices
  
  if ( !is.null(comp.names)){
    
    names(I) <- comp.names
    
  }else if (!is.null(names(observations))){
    
    names(I) <- names(observations)
    
  } else if (!is.null(names(predictions)) && is.null(comp.names) ){
    
    names(I) <- names(predictions)
    
  }
  
  return(Ic)
  
  #return(list(Total_index = Ic, Compartment_index= I))
  
}



#  absolute average fold error
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

#  average fold error
AFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- (log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}
