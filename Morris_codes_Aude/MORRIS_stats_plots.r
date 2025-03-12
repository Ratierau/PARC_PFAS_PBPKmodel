
#========== Library ============================================================
require(sensitivity)
require(gplots)

#========== Model information ==================================================
Model     <- "LifeTimePBPK_PFAS_Aug2023_flux_Foetus_volumes_LAST_FRACINTAKE_bdw_Fmilk_tempPFNA.model"
FichInput <- "Morris_mother_PFHxS.in"

# storage vectors 
fichierSim <- "Plan_Exp_SA.txt"
fichierOut <- "SA_Morris_mother_PFHxS3.out"

#===============================================================================
# (1) Parameter data collection
#===============================================================================
Nom.Parametre.Entier <- c()          

Nom.Parametre.Reel   <- c(
  "Free_no_pgcy",
  "Free_fet",
  "Hct",
  "Hct_foetus",
  "BirthYear",
  "HalfLife",
  "BDW_Var",
  "WeightGainPregnancy_Var",
  "BDW_foetus_Var",
  "Q_inj_rate",
  "Kd_uter2pla_frac_pop",
  "SD_Kd_uter2pla_frac_pop",
  "DecreaseIntake_1",
  "DecreaseIntake_2",
  "Abs",
  "PC_0",
  "PC_1",
  "PC_2",
  "PC_3",
  "PC_4",
  "PC_5",
  "PC_6",
  "PC_7",
  "PC_8",
  "PC_9",
  "PC_10",
  "PC_11",
  "PC_12",
  "PC_13",
  "PC_14",
  "PC_15",
  "PC_16",
  "PC_17",
  "PC_18",
  "PC_27",
  "PC_29",
  "PC_foetus_0",
  "PC_foetus_1",
  "PC_foetus_2",
  "PC_foetus_3",
  "PC_foetus_4",
  "PC_foetus_5",
  "PC_foetus_6",
  "PC_foetus_7",
  "PC_foetus_8",
  "PC_foetus_9",
  "PC_foetus_10",
  "PC_foetus_11",
  "PC_foetus_12",
  "PC_foetus_13",
  "PC_foetus_14",
  "PC_foetus_15",
  "PC_foetus_16",
  "PC_foetus_17",
  "PC_foetus_18"
  )

Nom.Parametre = c(Nom.Parametre.Entier ,  Nom.Parametre.Reel) 

  # (1.2) read initial value of the parameters =================================

Table.Par <- read.csv(paste("TablePar.csv",sep=""), header=T, sep=";")

ParList   <- Table.Par[,2]
Borne.inf <- Table.Par[,3]
Borne.sup <- Table.Par[,4]
N.par     <- nrow(Table.Par)

#===============================================================================
# (2) Sampling design of the Morris method
#===============================================================================
r      <- 1000 # Number of trajectory (in Morris method)
Nb_sim <- r * (N.par+1)              # Nb_Simulation Morris (without replicate)
levels <- 6                          # Number of discretisation

print(Nb_sim) #check

morris.design <- list(type = "oat", levels = levels, grid.jump = 3)

# Names of the outputs
output.names = c("C_ven", "C_foetus_ven")
# output.names = c("C_ven", "C_foetus_ven")
Time = c(16482680, 17482680.6968323, 17728925.830276, 17901727.830276, 18401727) #Output time in minutes

output.names = paste(rep(output.names, each=length(Time)), Time, sep="_") #generate output tag name with the associate time

n.output.names = length(output.names)

  # Parameter matrix collection ================================================

#morris function
Experience <- morris(model = NULL,
                     factors = Nom.Parametre,
                     r = r,
                     design = morris.design,
                     binf = Borne.inf,
                     bsup = Borne.sup,
                     scale=TRUE)

PlanExp = Experience$X
write.table(PlanExp ,file=fichierSim,sep="\t",quote=F)

## postscript("Distr_ParFull.ps")
pdf("Distr_ParFull.pdf")
par(mfrow=c(4,3))
for( i in 1:ncol(Experience$X)){ hist( Experience$X[,i], breaks=100, col="pink", main=colnames(Experience$X)[i] ) }
dev.off()

save(Experience, file = "ExperienceFull.RData")

#===============================================================================
# (3) parallel simulation                 
#===============================================================================
#system("makemcsim LifeTimePBPK_PFAS_women_fetus_PFHxS_PFNA.model") # if run on windows
system(paste0("./mcsim.LifeTimePBPK_PFAS_Aug2023_flux_Foetus_volumes_LAST_FRACINTAKE_bdw_Fmilk_tempPFNA Morris_mother_PFHxS.in SA_Morris_mother_PFHxS3.out")) #for job ccrt

#===============================================================================
#  (4) Results: outputs compilation
#===============================================================================

Data <- read.csv(fichierOut,dec=".",header=TRUE,sep="\t",na.strings = "na")
save(Data,         file = "DataFull.RData")

#N.par <- 60
SorScal <- Data[1:nrow(Data),  (N.par+2): ncol(Data)]
Sortie <- as.matrix(SorScal) # matrix, array

#select the number of output to be analysis
	#n.output.names: select all the output 
	#1:n -> n first outputs columns 
		Sortie.selected <- 1:n.output.names 

Sortie <- Sortie[, Sortie.selected ] #numeric

save(Sortie,         file = "SortieFull.RData")


Param <- Data[1:nrow(Data),2:(N.par+1)]
Parametre <- as.matrix(Param)
save(Parametre,         file = "SortieFull.RData")

Nom.Parametre=colnames(Param)
colnames(Sortie) = output.names

  # (4.1) Calcul: Package sensitivity ==========================================

sim.results.morris = t(Sortie)

mu = mu.star = sigma = NULL

for (i in (1:length(sim.results.morris[,1])))
{
  tell(Experience, sim.results.morris[i,])
  mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
  mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
  sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))

}

rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = output.names


write.table(mu, "muFull.csv" )  
write.table(mu.star,"mu_starFull.csv" )
write.table(sigma, "sigmaFull.csv" )

save(mu, file = "muFUll.RData" )  
save(mu.star, file = "mu_starFull.RData" )
save(sigma, file = "sigmaFull.RData" )

#===============================================================================
#  (5)             Graphical representation               
#===============================================================================

load("mu_starFull.RData")
load("sigmaFull.RData")

#mu_plot = mu#[-c(9,10,15,16),]
mu.star_plot = mu.star[-c(6,9,10),] #remove C_foetus_ven line output = 0
sigma_plot = sigma[-c(6,9,10),] #remove C_foetus_ven line output = 0

  # (5.1) For each output ======================================================
                                   
png(filename = "ResFigBrut%d.png", # compression = "lzw",
    width = 200, height = 200, units = "mm", # taille de l'image en mm
    res = 400, # resolution in dpi
    pointsize = 10,   bg = "white", family = "",  type = "cairo" )

NomParConv <- read.csv("TableParName.csv",dec=".",header=TRUE,sep=";",na.strings = "na")

names_brut = colnames(mu.star)
names=NULL
for(i in 1:length(names_brut)){ names= c( names, as.character( NomParConv[names_brut[i]==NomParConv[ ,2],1]) )}
names_base = names

par(mfrow=c(1,1))
for (k in 1:nrow(sigma_plot)) {
  xmin = min(mu.star_plot[k,]) - 0.2*min(mu.star_plot[k,])
  xmax = max(mu.star_plot[k,]) + 0.2*max(mu.star_plot[k,])
  
  ymin = min(sigma_plot[k,]) - 0.2*min(sigma_plot[k,])
  ymax = max(sigma_plot[k,]) + 0.2*max(sigma_plot[k,])
  
  min = min(xmin,ymin)
  max = max(xmax,ymax)
  
  plot((mu.star_plot[k,]),(sigma_plot[k,]) , pch = 16, col = 2, xlim = c(min, max), ylim = c(min, max), 
       main = rownames(sigma_plot)[k], xlab = "mu.star", ylab = "sigma")
  
  
  pos_alea = rep( c(3,4), length(colnames(sigma_plot)) )[1: length(colnames(sigma_plot))]
  text(mu.star_plot[k,], sigma_plot[k,], colnames(sigma_plot), cex = 0.6, pos = pos_alea , srt = 50, offset=0.5)
  
  abline(b=1,a=0,col="grey", lty=4)
  
}


dev.off()

  # (5.2) Calculation of the global index ======================================
mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)

IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )

xmin = min(IGlobal_mu.star) - 0.05*min(IGlobal_mu.star)
xmax = max(IGlobal_mu.star) + 0.05*max(IGlobal_mu.star)

ymin = min(IGlobal_sigma) - 0.05*min(IGlobal_sigma)
ymax = max(IGlobal_sigma) + 0.05*max(IGlobal_sigma)

min = min(xmin,ymin)
max = max(xmax,ymax)


png(filename = "ResFigGlobal_zoom3.png",
    width = 200, height = 200, units = "mm",
    res = 400, # resolution in dpi
    pointsize = 10,   bg = "white", family = "",  type = "cairo" )


  
plot(IGlobal_sigma  ~ IGlobal_mu.star, pch = 16, col = 2, xlim=c(min,20), 
    ylim=c(min,15), xlab = "mu.star", ylab = "sigma", main = "Global Index")


#plot(IGlobal_sigma  ~ IGlobal_mu.star, pch = 16, col = 2, xlim=c(min,max), ylim=c(min,max), xlab = "mu.star", ylab = "sigma", main = "Global Index")
text(IGlobal_mu.star, IGlobal_sigma , names(IGlobal_mu.star), cex = 0.6, pos = 4)
x=seq(-5,20)
y=seq(-5,20)
fit= lm(y~x)
abline(fit,col="grey", lty=4)

dev.off()


# zoom it
png(filename = "ResFigGlobal_zoom8.png",
    width = 200, height = 200, units = "mm",
    res = 400, # resolution in dpi
    pointsize = 10,   bg = "white", family = "",  type = "cairo" )

plot(IGlobal_sigma  ~ IGlobal_mu.star, pch = 16, col = 2, xlim=c(min,6), 
     ylim=c(min,5), xlab = "mu.star", ylab = "sigma", main = "Global Index")

#plot(IGlobal_sigma  ~ IGlobal_mu.star, pch = 16, col = 2, xlim=c(min,max), ylim=c(min,max), xlab = "mu.star", ylab = "sigma", main = "Global Index")
text(IGlobal_mu.star, IGlobal_sigma , names(IGlobal_mu.star), cex = 0.6, pos = 4)
x=seq(-5,20)
y=seq(-5,20)
fit= lm(y~x)
abline(fit,col="grey", lty=4)

dev.off()
# 
# # plotly
# df_morris <- data.frame(IGlobal_sigma, IGlobal_mu.star)
# 
# ggplotly(ggplot(df_morris,aes( x= IGlobal_mu.star, y=IGlobal_sigma, label = rownames(df_morris)))+
#            geom_point()+
#            theme_minimal())
# 
# ggplot(df_morris, aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = rownames(df_morris))) +
#   geom_point() +
#   geom_text(size = 3, vjust = -0.5, hjust = 0.5) +  # Taille réduite du texte+
#   ylim(0.85,20)+
#   xlim(0.85,20)+
#   theme_minimal()


library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)

#==============================================
# (0) Select parameters 4 Sobol
#==============================================

# Sigma
selected_cols <- sigma_scale >= 1  
filtered_mat_sigma <- sigma_scale
filtered_mat_sigma[!selected_cols] <- NA
View(filtered_mat_sigma)
# Mu
selected_cols <- mu.star_scale >= 1  
filtered_mat_mu <- mu.star_scale
filtered_mat_mu[!selected_cols] <- NA
View(filtered_mat_mu)
filtered_mat_mu <- as.data.frame(filtered_mat_mu)

 
# Plot barplot sigma
df_long <- as.data.frame(sigma_scale) %>%
  rownames_to_column(var = "Sample") %>%   # Convertir les noms de ligne en colonne
  pivot_longer(cols = -Sample,             # Toutes les colonnes sauf "Sample"
               names_to = "Variable", 
               values_to = "Value") %>%
  filter(!is.na(Value))  # Supprimer les valeurs NA

for (sample_id in unique(df_long$Sample)) {
  
  df_sample <- df_long %>% filter(Sample == sample_id) %>%
    mutate(Color = ifelse(Value < 1, "Blanc", "Gris")) 
  
  p <- ggplot(df_sample, aes(x = reorder(Variable, Value), y = Value, fill = Color)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 1,  color = "red", linewidth = 1) +  # Ligne verticale rouge
    scale_fill_manual(values = c("Blanc" = "gray", "Gris" = "black")) +  # Attribution des couleurs
    coord_flip() +  # Rotation pour lisibilité
    theme_minimal() +
    labs(title = paste(sample_id),
         x = "Parameters",
         y = "sigma_scale") +
    theme(legend.position = "none")  
  
  print(p)
  ggsave(filename = paste0("barplot_sigma_", sample_id, ".png"), plot = p, width = 8, height = 5)
}

# Plot barplot mu.star
df_long <- as.data.frame(mu.star_scale) %>%
  rownames_to_column(var = "Sample") %>%   # Convertir les noms de ligne en colonne
  pivot_longer(cols = -Sample,             # Toutes les colonnes sauf "Sample"
               names_to = "Variable", 
               values_to = "Value") %>%
  filter(!is.na(Value))  # Supprimer les valeurs NA

for (sample_id in unique(df_long$Sample)) {
  
  df_sample <- df_long %>% filter(Sample == sample_id) %>%
    mutate(Color = ifelse(Value < 1, "Blanc", "Gris")) 
  
  p_s <- ggplot(df_sample, aes(x = reorder(Variable, Value), y = Value, fill = Color)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 1,  color = "red", linewidth = 1) +  # Ligne verticale rouge
    scale_fill_manual(values = c("Blanc" = "gray", "Gris" = "black")) +  # Attribution des couleurs
    coord_flip() +  # Rotation pour lisibilité
    theme_minimal() +
    labs(title = paste(sample_id),
         x = "Parameters",
         y = "mu.star_scale") +
    theme(legend.position = "none")  
  
  print(p_s)
  ggsave(filename = paste0("barplot_mustar_", sample_id, ".png"), plot = p, width = 8, height = 5)
}


# (5.2.1) Calculation of the global index for each time ========================
# chemin = "../"

#source(paste0("global_index_time.R") )

# (5.2.2) Calculation of the global index for each output ======================

#source(paste0("global_index_output.R") )

# TO VISUALIZE MANUALLY THE POINTS ==============================================

# library(ggplot2)
# library(plotly)
# IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star)
# P <- ggplot(IGlobal) +
#   aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = rownames(IGlobal)) +
#   xlim (min,5) + #to zoom in/out, determine max value
#   ylim (min,5) + #to zoom in/out, determine max value
#   geom_point() 
# 
# D <- ggplotly(P)
# D

