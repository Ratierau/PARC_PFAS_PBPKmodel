
# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)
library(pksensi)

# SENSITIVITY ANALYSIS ####

## Define the distribution of the parameters that you will analyse in the sensitivity test 
## q: vector that contains all the parameters to be tested and defines their distribution 
## For now the distribution is assumed to be unified for all of them, but I propably need to change this
q <- c("qunif", # VInL
       "qunif", # VBi
       "qunif", # VIn
       "qunif", # VLi
       "qunif", # VKi
       "qunif", # VAd
       "qunif", # VSl
       "qunif", # VRa
       "qunif", # VBl
       "qunif", # PInb
       "qunif", # PLib
       "qunif", # PKib
       "qunif", # PAdb
       "qunif", # PSlb
       "qunif", # PRab
       "qunif", # fub
       "qunif", # PInb.glu
       "qunif", # PLib.glu
       "qunif", # PKib.glu
       "qunif", # PAdb.glu
       "qunif", # PSlb.glu
       "qunif", # PRab.glu
       "qunif", # fub.glu 
       "qunif", # QIn
       "qunif", # QPo
       "qunif", # QLi
       "qunif", # QHa
       "qunif", # QKi
       "qunif", # QAd
       "qunif", # QSl
       "qunif", # QRa
       "qunif", # ka
       "qunif", # kTIn
       "qunif", # kTge
       "qunif", # CL_Re
       "qunif", # CL_Bi
       "qunif", # VmaxMet.In
       "qunif", # KmGluc.In
       "qunif", # VmaxMet.Li
       "qunif", # KmGluc.Li
       "qunif"  # kHydro
)


## Set parameter distribution ##
# we use 10% change in all parameters

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit

## q.arg: lists all the values of each parameter to be tested
q.arg <- list(list(min = PARAMETERS["VInL"]*LL, max= PARAMETERS["VInL"]*UL),
              list(min = PARAMETERS["VBi"]*LL, max= PARAMETERS["VBi"]*UL),
              list(min = PARAMETERS["VIn"]*LL, max =  PARAMETERS["VIn"]*UL),
              list(min = PARAMETERS["VLi"]*LL, max =  PARAMETERS["VLi"]*UL),
              list(min = PARAMETERS["VKi"]*LL, max =  PARAMETERS["VKi"]*UL),
              list(min = PARAMETERS["VAd"]*LL, max = PARAMETERS["VAd"]*UL),
              list(min = PARAMETERS["VSl"]*LL, max = PARAMETERS["VSl"]*UL),
              list(min = PARAMETERS["VRa"]*LL, max = PARAMETERS["VRa"]*UL),
              list(min = PARAMETERS["VBl"]*LL, max = PARAMETERS["VBl"]*UL),
              list(min = PARAMETERS["PInb"]*LL, max = PARAMETERS["PInb"]*UL),
              list(min = PARAMETERS["PLib"]*LL, max = PARAMETERS["PLib"]*UL),
              list(min = PARAMETERS["PKib"]*LL, max = PARAMETERS["PKib"]*UL),
              list(min = PARAMETERS["PAdb"]*LL, max = PARAMETERS["PAdb"]*UL),
              list(min = PARAMETERS["PSlb"]*LL, max = PARAMETERS["PSlb"]*UL),
              list(min = PARAMETERS["PRab"]*LL, max = PARAMETERS["PRab"]*UL),
              list(min = PARAMETERS["fub"]*LL, max = PARAMETERS["fub"]*UL),
              list(min = PARAMETERS["PInb.glu"]*LL, max = PARAMETERS["PInb.glu"]*UL),
              list(min = PARAMETERS["PLib.glu"]*LL, max = PARAMETERS["PLib.glu"]*UL),
              list(min = PARAMETERS["PKib.glu"]*LL, max = PARAMETERS["PKib.glu"]*UL),
              list(min = PARAMETERS["PAdb.glu"]*LL, max = PARAMETERS["PAdb.glu"]*UL),
              list(min = PARAMETERS["PSlb.glu"]*LL, max = PARAMETERS["PSlb.glu"]*UL),
              list(min = PARAMETERS["PRab.glu"]*LL, max = PARAMETERS["PRab.glu"]*UL),
              list(min = PARAMETERS["fub.glu"]*LL, max = PARAMETERS["fub.glu"]*UL), 
              list(min = PARAMETERS["QIn"]*LL, max = PARAMETERS["QIn"]*UL),
              list(min = PARAMETERS["QPo"]*LL, max = PARAMETERS["QPo"]*UL),
              list(min = PARAMETERS["QLi"]*LL, max = PARAMETERS["QLi"]*UL),
              list(min = PARAMETERS["QHa"]*LL, max = PARAMETERS["QHa"]*UL),
              list(min = PARAMETERS["QKi"]*LL, max = PARAMETERS["QKi"]*UL),
              list(min = PARAMETERS["QAd"]*LL, max = PARAMETERS["QAd"]*UL),
              list(min = PARAMETERS["QSl"]*LL, max = PARAMETERS["QSl"]*UL),
              list(min = PARAMETERS["QRa"]*LL, max = PARAMETERS["QRa"]*UL),
              list(min = PARAMETERS["ka"]*LL, max = PARAMETERS["ka"]*UL),
              list(min = PARAMETERS["kTIn"]*LL, max = PARAMETERS["kTIn"]*UL),
              list(min = PARAMETERS["kTge"]*LL, max = PARAMETERS["kTge"]*UL),
              list(min = PARAMETERS["CL_Re"]*LL, max = PARAMETERS["CL_Re"]*UL),
              list(min = PARAMETERS["CL_Bi"]*LL, max = PARAMETERS["CL_Bi"]*UL),
              list(min = PARAMETERS["VmaxMet.In"]*LL, max = PARAMETERS["VmaxMet.In"]*UL),
              list(min = PARAMETERS["KmGluc.In"]*LL, max = PARAMETERS["KmGluc.In"]*UL),
              list(min = PARAMETERS["VmaxMet.Li"]*LL, max = PARAMETERS["VmaxMet.Li"]*UL),
              list(min = PARAMETERS["KmGluc.Li"]*LL, max = PARAMETERS["KmGluc.Li"]*UL),
              list(min = PARAMETERS["kHydro"]*LL, max = PARAMETERS["kHydro"]*UL)
)



## Create parameter matrix ##  
params <- c("VInL",
            "VBi",
            "VIn",
            "VLi",
            "VKi",
            "VAd",
            "VSl",
            "VRa",
            "VBl",
            "PInb",
            "PLib",
            "PKib",
            "PAdb",
            "PSlb",
            "PRab",
            "fub",
            "PInb.glu",
            "PLib.glu",
            "PKib.glu",
            "PAdb.glu",
            "PSlb.glu",
            "PRab.glu",
            "fub.glu", 
            "QIn",
            "QPo",
            "QLi",
            "QHa",
            "QKi",
            "QAd",
            "QSl",
            "QRa",
            "ka",
            "kTIn",
            "kTge",
            "CL_Re",
            "CL_Bi",
            "VmaxMet.In",
            "KmGluc.In",
            "VmaxMet.Li",
            "KmGluc.Li",
            "kHydro")

length(params)==length(q)  # Check if lengths of params and q match

set.seed(1234) # For reproducibility ?! not sure what this is
x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 10, conf = 0.95) # Generate the parameter matrix
dim(x$a) # The array of c(model evaluation, replication, parameters)


## Conduct simulation ##
# out <- solve_fun(x, time=times, func = PBPKmodPFOA, initState = yini, outnames = outputs)

out <- solve_fun(x,  
                 times = times,  
                 func = PBK_ZENandGLUC,
                 initState = INITIALstate, 
                 outnames = outputs
)
# head(PBK)


saveRDS(object=out, file="out_scaled.rds") # Save the output
#out <- readRDS("out_scaled.rds")

## Output of the Uncertainty analysis ##
pdf("out_scaled.pdf") # Open a PDF device
pksim(out) # Output the results of the uncertainty analysis
dev.off() # Close the PDF device

## Output from the sensitivity analysis ##
pdf("out_scaled.pdf")
plot(out)
dev.off()

ResultsSI <- as.data.frame(print(out["tSI"])) # Extract and convert sensitivity indices to a data frame
ResultsSI$Times <- rownames(ResultsSI)  # Get times as row names
ResultsSI$Times <- as.numeric(as.character(ResultsSI$Times)) # Convert times to numeric

write.xlsx(ResultsSI,
           file =file.path(newday,"ResultsSI.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(ResultsSI,
           file = "Results/2023-09-24/ResultsSI.xlsx",
           colNames = TRUE
)

check(out) #Perform checks on the output?!
pdf("heat_check_CI_scaled.pdf") #Open a PDF device
heat_check(out, index = "CI") #Plot all heat checks
dev.off() #Close the PDF device

pdf("heat_check_all.pdf") #Open a PDF device
heat_check(out, show.all = TRUE) #Plot all heat checks
dev.off() #Close the PDF device


