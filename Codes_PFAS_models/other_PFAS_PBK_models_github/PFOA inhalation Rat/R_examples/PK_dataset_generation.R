library(deSolve) # Imports ODE Solvers
library(ggplot2) # Creating plots

# check the working directory 
getwd()

# Change the working directory
# add whatever working directory you want
wd = '/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/R_examples' 
setwd(wd)
getwd()

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

#------------#
# User Input #
#------------#

user_input <- list('k12_value'=0.5, 
                   'k21_value'=0.2, 
                   'Mx_initial'= 1,
                   'My_initial'= 0.1,
                   'IV_intakes'=c(1, 1),
                   'IV_times'=c(5, 7.5))

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.01)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           #events = events,
                           method="lsodes",rtol = 1e-05, atol = 1e-05))

# Check the mass balance: the sum of each row of the solution dataframe should be
# equal to the sum of initial amounts.
#rowSums(solution[ ,-1])

# Basic Plot of Mx vs time
# plot(x = solution$time, y = solution$Mx)

# Plot with ggplot2

compartments <- c('Mx', 'My', 'M_eliminated')
color_codes <- scales::hue_pal()(length(compartments))

plot <- ggplot(data = solution)+
  geom_line( aes(x = time, y = Mx, color='Mx'), size=1.3)+
  geom_line( aes(x = time, y = My, color='My'), size=1.3)+
  geom_line( aes(x = time, y = M_eliminated, color='M_eliminated'), size=1.3)+
  
  labs(title = 'Predicted values',
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

# Save the generated data

# Define the directory at which you want to save the data
wd_data <- '/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/R_examples'
setwd(wd_data)

# Select which time points you want to save
t_selected <- seq(1,10,1)

# Save the 'solution' dataframe into an xlsx (excel) file
openxlsx::write.xlsx(x = solution[solution$time %in% t_selected, ], file = 'Dataset/PK_data.xlsx')

# Or save the 'solution' dataframe as a csv file
write.csv(x = solution, file = 'Dataset/PK_data.csv', row.names = FALSE)
