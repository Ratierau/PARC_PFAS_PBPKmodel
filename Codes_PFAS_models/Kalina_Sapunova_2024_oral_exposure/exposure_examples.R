setwd("C:/Jirka/Environment/PARC/WP7/Task7.3/Task 7.3.2/External_exposure_model/R package/Code")

#########################
#                       #
#  Attaching libraries  #
#                       #
#########################

library("fitdistrplus")  # This library is necessary for the MLE estimates.



##################
#                #
#  Loading data  #
#                #
##################

load("../Data/food_concentrations.rds")         # This is an rds file with PFASs concentrations in food items. Prepared by Daria.
#load("../Data/food_concentrations_example.rds")  # Fake data for example use.
load("../Data/food_freq.rds")                   # This is an rds file with food item frequencies. Prepared by food.consumption.generate.R file from data provided by Pavel Piler.
#load("../Data/food_freq_example.rds")            # Fake data for example use.
load("../Data/portions_masses.rds")               # This is an rds file with portion sizes (in quantiles). Prepared by Daria and adjusted by prepare.portions.R.


#######################
#                     #
#  Loading functions  #
#                     #
#######################

source("exposure.R")                        # The main function returning an exposure for one compound, one food item and one person based on estimated lognormal distribution.
source("exposure.mass.R")                   # The function returning an exposure for multiple compounds and multiple food items based on sum of exposure() results and mass/time consumption.
source("exposure.freq.R")                   # The function returning an exposure for multiple compounds and multiple food items based on sum of exposure() results and a frequency of food consumption + portion sizes.
source("exposure.cohort.R")                 # The function for modelling whole cohorts (frequency-based).
source("fenton.wilkinson.R")                # The function for modelling whole cohorts (frequency-based).


#####################################
#                                   #
#  Testing the function exposure()  #
#                                   #
#####################################

# Results for a given country (CZ)
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ",year=2004,data=food_concentrations) # Fully specified, time dependent result.
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ",year=2004)                          # If there is a variable called food_concentrations in the computing environment, it is not necessary to specify the "data" argument.
exposure("Z0008.0001.0015","Z0008.0001","ΣPFOA",country="CZ",year=2004)                     # Trying to get results for "Z0008.0001.0015" in Czechia. It is not available, so level 2 EFSA expo ("Z0008.0001") is used as a surrogate.
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ",year=2005)                          # Another year (compare!).
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ",year=2024)                          # Another year.
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ")                                    # Time independent result (all data are used, no trend computed).
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="CZ",start=2012)                         # Only data from 2012 are used to estimate the exposure.

# Results for a given region (Easter Europe)
exposure("Z0008.0001","Z0008.0001","ΣPFOA",region="Eastern Europe")                         # As in the previous examples, there is not enought data in Eastern Europe and the whole Europe is used for the estimated exposure.

# Let's try some data rich country (the Netherlands).
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="NL",year=2024)                          # Time-trend estimated result.
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="NL",start=2014,year=2024)               # Result for last 10 years are not available since all measurements are from 2021.
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="NL",start=2014)                         # But we can get static result based on data from the same period (last 10 years).
exposure("Z0008.0001","Z0008.0001","ΣPFOA",country="NL")                                    # Slightly different result over the whole dataset.

hist(food_concentrations$median[which(food_concentrations$expoHierarchyCode.x=="Z0008.0001"&food_concentrations$chemical=="ΣPFOA"&food_concentrations$country=="NL")]) # Just to see what is in the DB.

exposure("Z0008.0001","Z0008.0001","ΣPFOA",region="Western Europe")                         # The whole region of Western Europe.
hist(food_concentrations$median[which(food_concentrations$expoHierarchyCode.x=="Z0008.0001"&food_concentrations$chemical=="ΣPFOA"&food_concentrations$Europe_region=="Western Europe")],freq=FALSE) # Let's visualize it.
hist(rlnorm(1000,log(exposure("Z0008.0001","Z0008.0001","ΣPFOA",region="Western Europe",method="analytic")["gmean"]),                                                                               # Adding a simulated data based on analytic result.
                 log(exposure("Z0008.0001","Z0008.0001","ΣPFOA",region="Western Europe",method="analytic")["gsd"])),add=TRUE,col=rgb(0.6,0.6,0.2,0.5),freq=FALSE)


exposure(NA,"Z0007.0001","ΣPFOS",region="Northern Europe",fit="moments")                     
exposure(NA,"Z0007.0001","ΣPFOS",region="Northern Europe",fit="mle")

exposure(NA,"Z0007.0001","ΣPFOS",region="Northern Europe",fit="moments",method="analytic")  # Using analytic output (characteristics of the log-normal distribution) with basic fit by 1st and 2nd moment.
exposure(NA,"Z0007.0001","ΣPFOS",region="Northern Europe",fit="mle",method="analytic")      # Using analytic output (characteristics of the log-normal distribution) with MLE fit.



##########################################
#                                        #
#  Testing the function exposure.mass()  #
#                                        #
##########################################

# Defining food.mass.individual in g/day of each food
food_mass_individual<-data.frame("expoHierarchyCode.x"=c("Z0001.0001.0002.0012","Z0001.0001.0002.0012","Z0001.0002","Z0002.0008","Z0005.0001.0002.0002.0001","Z0013.0001"),
                                 "expoHierarchyCode.2"=c("Z0001.0001","Z0001.0001","Z0001.0002","Z0002.0008","Z0005.0001","Z0013.0001"),
                                 "gmean"              =c(100,100,100,100,100,100),
                                 "gsd"                =c(2,1.2,1.5,1.5,0.5,1.5))  

exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),region="Northern Europe",method="numeric",n=10000)
exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),region="Northern Europe",method="analytic")

# Defining food.mass.individual in g/day of each food
food_mass_individual<-data.frame("expoHierarchyCode.x"=c("Z0008.0001"),
                                 "expoHierarchyCode.2"=c("Z0008.0001"),
                                 "gmean"              =c(1000),
                                 "gsd"                =c(2))  

exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),region="Northern Europe",method="numeric",n=10000)
exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),region="Northern Europe",method="analytic")


# Be careful in case there are very high values of GSD in some food items, the Fenton-Wilkinson strives.
food_mass_individual<-data.frame("expoHierarchyCode.x"=c("Z0007.0001","Z0019.0001.0001.0001"),
                                 "expoHierarchyCode.2"=c("Z0007.0001","Z0019.0001"),
                                 "gmean"              =c(300,160),
                                 "gsd"                =c(1.5,1.4)) 

exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),country="NO",method="numeric",n=10000)
exposure.mass(food_mass_individual=food_mass_individual,c("ΣPFOA","ΣPFOS"),country="NO",method="analytic")



##########################################
#                                        #
#  Testing the function exposure.freq()  #
#                                        #
##########################################

# Example of individual food frequencies food.freq.individual created from food.freq cohort table
food_freq_individual<-food_freq[,c("expoHierarchyCode.x","expoHierarchyCode.2","individual.1")]
colnames(food_freq_individual)[3]<-"freq"

exposure.freq(food_freq_individual,portions_masses,c("ΣPFOA","ΣPFOS"),country="FR")
exposure.freq(food_freq_individual,portions_masses,c("ΣPFOA","ΣPFOS"),country="FR",method="analytic")
