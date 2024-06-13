library(deSolve) # Imports ODE Solvers
library(ggplot2) # Creating plots
library(nloptr)  # Optimization algorithms

# Change the working directory
wd = '/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/R_examples' 
setwd(wd)

#-----------------------------------------------#
# Custom functions for the solution of the ODEs #
#-----------------------------------------------#

# create.params() returns a list with all user defined parameters

create.params <- function(user_input){
  with(as.list(user_input),{
    k12 = k12_value 
    k21 = k21_value
    ke = 0.3
    
    return(list('k12' = k12, 'k21' = k21, 'ke'=ke,
                'Mx_initial'=Mx_initial, 'My_initial'=My_initial, 
                'IV_intakes'=IV_intakes, 'IV_times'=IV_times))
  })
}

# ode.func() is a function with ODEs of the model
ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    dMx = k21 * My - k12 * Mx - ke * Mx
    
    dMy = - k21 * My + k12 * Mx 
    
    dM_eliminated = ke * Mx
    
    # Important: The derivatives must be returned with the same row
    # that they are reported 
    return(list(c('dMx'=dMx, 'dMy'=dMy, 'dM_eliminated'=dM_eliminated)))
  })
}


# create.inits() returns a vector with the initial conditions
# provided to the ODEs
create.inits <- function(parameters){
  with(as.list(parameters),{
    
    Mx <- Mx_initial
    My <- My_initial
    M_eliminated = 0
    
    return(c('Mx'= Mx, 'My' = My, 'M_eliminated' = M_eliminated))
  })
}

# create.events() returns a dataframe with the dosing plan 
create.events <- function(parameters){
  with(as.list(parameters),{
    
    IV_intakes <- IV_intakes
    IV_times <- IV_times
    IV_ltimes <- length(IV_times)
    
    events <- data.frame(var = c(rep('Mx', IV_ltimes)), 
                         time = c(IV_times),
                         value = c(IV_intakes),
                         method = rep('add',IV_ltimes))
    #events <- events[order(events$time),] 
    return(list(data=events))
  })
}

#------------------------#
# Load Experimental Data #
#------------------------#
exp_data <- openxlsx::read.xlsx(xlsxFile = 'Dataset/PK_data.xlsx')

#--------------------#
# Mean-Squared Error #
#--------------------#

mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

#-------------------------------#
# Define the Objective function # 
#-------------------------------#

# obj.func(): A function that takes as input a vector of values
# that will be assigned to the unknown parameters. It returns 
# a the value of a a metric function, that quantifies how well are 
# the experimental data described by the model, using the 
# current values of x.

obj.func <- function(x){
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  
  # Assign the x values to the corresponding parameters
  user_input <- list('k12_value'=x[1], 
                     'k21_value'=x[2], 
                     'Mx_initial'= 1,
                     'My_initial'= 0.1,
                     'IV_intakes'=NULL,
                     'IV_times'=NULL)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  # events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.01)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                      #events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for x and y for the time points 
  # at which we have available data. 

  preds <- solution[solution$time %in% exp_data$time, ]
  
  # Initiate a variable to save the score of goodness of fit for each variable
  score_x <- mse_custom(exp_data$Mx, preds$Mx)
  score_y <- mse_custom(exp_data$My, preds$My)
    
  mean_score <- mean(c(score_x, score_y))
  
  return(mean_score)
}

#---------------------------------#
# Set up the Optimization process #
#---------------------------------#

# Initial values of the unknown parameters k12 and k21
x0 <- c(100,100)

# Maximum number of iterations of the optimization algorithm
N_iter <- 1000

# Extra options for the optimization algorithm
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",  #"NLOPT_LN_SBPLX" ,
              "xtol_rel" = 1e-05,
              "ftol_rel" = 1e-05,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = N_iter,
              "print_level" = 1 )

# Call the optimization algorithm and provide him with the input data
optimization <- nloptr::nloptr(x0 = x0,
                               eval_f = obj.func,
                               lb	= rep(1e-03, length(x0)),
                               ub = rep(1e+03, length(x0)),
                               opts = opts)

# The "optimization" object contains the results of the optimization process
# as well as the final values of the optimized parameters

# The minimized value of the objective function
optimization$objective

# The values of the optimized params 
x_opt <- optimization$solution

#---------------------------------------------#
# Plot predictions over the experimental data #
#---------------------------------------------#

# Step 1: Solve the ODEs using the optimized values of parameters k12 and k21
user_input <- list('k12_value'=x_opt[1], 
                   'k21_value'=x_opt[2], 
                   'Mx_initial'= 1,
                   'My_initial'= 0.1,
                   'IV_intakes'=NULL,
                   'IV_times'=NULL)
params <- create.params(user_input)
inits <- create.inits(params)
# events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.01)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    #events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

compartments <- c('Mx', 'My', 'M_eliminated')
color_codes <- scales::hue_pal()(length(compartments))

plot <- ggplot()+
  geom_line(data = solution, aes(x = time, y = Mx, color='Mx'), size=1.3)+
  geom_line(data = solution, aes(x = time, y = My, color='My'), size=1.3)+
  geom_line(data = solution, aes(x = time, y = M_eliminated, color='M_eliminated'), size=1.3)+
  
  geom_point(data = exp_data, aes(x = time, y = Mx, color='Mx'), size=5)+
  geom_point(data = exp_data, aes(x = time, y = My, color='My'), size=5)+
  geom_point(data = exp_data, aes(x = time, y = M_eliminated, color='M_eliminated'), size=5)+
  
  labs(title = 'Predicted vs Observed Values',
       y = 'Mass (ug)' , x = "Time (hours)")+
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25),
        legend.text=element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
  
  scale_color_manual("Compartments", values=color_codes)+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))

print(plot)
