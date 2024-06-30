########################################################################
#                                                                      #
#  This is an elementary introduction to the external exposure model.  #
#             Made by Jiří Kalina and Daria Sapunova 2024.             #
#                                                                      #
########################################################################

##################
#                #
#  Introduction  #
#                #
##################

# The model is based on functions defined in separate files and different types of input data:
#
# In version 0, full real datasets are not provided. Examples for testing purposes are available until RCX publishes the PFAS data.
#
#
# 1. Dataset of PFASs concentration in food items. Dataset made by Lisa Melymuk and Daria Sapunova (Europe) will be available in the file food_concentrations.rds.
#    An example dataset for testing purposes is available in the file food_concentrations_example.rds 
#    The strucuture of the dataset is described in the file data_pfas_descriptiopn.R.
#
# 2. Dataset of food consumption in g/day. In this version we don't have a file with real data yet, so it is necessary to create your own dataset.
#    For each individual it is a data frame with four columns:
#      "expoHierarchyCode.x" EFSA expo hierarchy code level > 2
#      "expoHierarchyCode.2" EFSA expo hierarchy code level 2
#      "gmean"               geometric mean of the consumption in g/day           
#      "gsd"                 geometric standard deviation of the consumption in g/day
#
# 3. Alternatively, the frequency of consumption may be used, together with portion sizes.
#    The frequency of consumption is stored in portions/year. Will be available in the file food_freq.rds.
#    In version 0, fake data is avaiable as the file food_freq_example.rds, containing 3 persons only.
#
# 4. Portion sizes. Contains same columns as food consumption but in g/portion.
#    Portion masses are based on EFSA data. Lognormal distribution is expected. Available in the file portions_masses.rds.


#########################
#                       #
#  Available functions  #
#                       #
#########################

# Depending on the level of detail of the exposure (and the input variables), four functions are available:
#
#
# 1. Elementary function exposure() which returns an exposure for 1 g of a food item (specified either by expoHierarchyCode.x or expoHierarchyCode.2).
#    Estimates the exposure based on data for selected country or region.
#    If "allow.regions" switch is on (implicitly), jumps from country to region and from region to the whole of Europe in case of insufficient data.
#    Allows either to use all data (static) or define concrete year of exposure (trend-based estimation) or limit the period for data (start to end years).
#    Acommodates two methods for estimating log-normal distribution of PFAss concentration in food item and country.: most-likelihood estimate (MLE) and
#    simple estimation of geometric mean and geometric standard deviation based on sample values (from table food_concentrations.rds).
#
# 2. Function exposure.mass() as a wrapper for exposure() and multiple PFASs and multiple food items. If the mass of each item consumed per time unit is available,
#    estimates the cumulative exposure. Either by Monte-Carlo approcah or by using adjusted Fenton-Wilkinson approximation for sum of independent log-normal distributions.
#
# 3. Alternative, food frequency and portion sizes may be used as arguments of exposure.freq function. For an individual, frequencies of different food items and
#    log-normal distribution of portion sizes may be used as inputs of this function. Similarly to exposure.mass() the function either uses numeric Monte-Carlo
#    simulation or uses Fenton-Wilkinson approximation for summing multiples of log-normal distributions (portion × concentration.)
#
# 4. Wrapper of exposure.freq() called exposure.cohort() allows to repeated the frequentist approach for the whole cohort. Not debugged yet, in version 0 just for
#    testing. Will be adapted for cohort data from the PARC T7.3.2 project.


##########################
#                        #
#  Functions' arguments  #
#                        #
################x#########

# Detailed description of functions' arguments and usegae is available in file exposure_examples.R.