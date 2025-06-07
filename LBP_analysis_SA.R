# 5/23/2025
# LBP analysis 

# -- set directories 
rm(list = ls())

################  
library(rhdf5) 
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

#setwd("/Users/hyung/Dropbox/LBP/github_mOFC_signature")
dat0 <-  read.csv("phenotype.csv")
# Redefine the variable by replacing "Widespread Pain" with "Widespread Hyperalgesia"
dat0$Pain_Phenotype <- ifelse(dat0$Pain_Phenotype == "Widespread Pain",
                              "Widespread Hyperalgesia",
                              dat0$Pain_Phenotype)

summary(dat0)
Subject <- dat0$Subject
Pain_Phenotype <- dat0$Pain_Phenotype 
 
 
dat2save <- readRDS("LBP_EEG_05062025.RData")  
eeg_dat_Y <- dat2save$eeg_dat_Y
eeg_dat_J <- dat2save$eeg_dat_J
 
## figures 

# -- All conditions and subset to condition_select
conditions <- c("Back_32", "Back_128", "Back_256", "Hand_32", "Back_128", "Hand_256")
condition_ids <- c(1, 3, 4, 6)
condition_select <- conditions[condition_ids]  # Selected subset
condition_select <- unique(condition_select)   # Just to be safe, remove any duplicates

# -- ROIs
ROI_labels = roi_s <- c("rACC-lh", "dACC-lh", "S1-lh",
                        "Ins-lh", "dlPFC-lh", "mOFC-lh",
                        "rACC-rh", "dACC-rh", "S1-rh",
                        "Ins-rh", "dlPFC-rh", "mOFC-rh")

pain_levels <- c("Widespread Hyperalgesia", "Localized Pain", "No Pain")  
dat0$Pain_Phenotype <- factor(dat0$Pain_Phenotype, levels = pain_levels)

tt <- 17:43  # corresponds to [-0.3, 0.7] sec 
ff <- 1:45
TT <- length(tt)
FF <- length(ff)

# time sequence from -0.3 to 0.7 sec
tt_sec <- seq(-0.3, 0.7, length.out = TT) 
desired_times <- c(-0.2, 0, 0.2, 0.4, 0.6) 
breaks_labels <- as.character(desired_times)


#### Group_Average plots 

# mOFC_lh
region = 6
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    
    if(conditions[j] %in% condition_select && eeg_dat_J[i,j]){   
      
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Pain_Phenotype
      subject_group <- dat0$Pain_Phenotype[dat0$Subject == Subject[i]]
      
      if(length(subject_group) == 1){
        grid$group <- rep(subject_group, length(tt) * length(ff))
      }else{
        grid$group <- rep(NA, length(tt) * length(ff))
      }
      
      dat2plot <- rbind(dat2plot, grid)
    }
  }
}

dat2plot$group <- factor(dat2plot$group, levels = pain_levels) 
dat2plot$condition <- factor(dat2plot$condition, levels = condition_select)
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Back_32" = "Back 32", 
                                    "Back_256"= "Back 256", 
                                    "Hand_32" = "Hand 32",
                                    "Hand_256" = "Hand 256") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))


df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = mean(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

# Plot the results 
pp_mOFC_lh <- ggplot(df_summary, aes(x = time_sec, y = ff, fill = avg_value)) +
  geom_tile() + 
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0  
  ) + 
  facet_grid(group ~ condition) + 
  labs(x = "Time (s)",
       y = "Frequency (Hz)",
       #fill = "dB change") +
       fill = "Value") + 
  scale_x_continuous(breaks = desired_times, labels = breaks_labels) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major.x = element_blank(),   
        panel.grid.minor.x = element_blank()) + 
  ggtitle(paste(ROI_labels[region]))

plot(pp_mOFC_lh) 





# mOFC_rh
region = 12
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    
    if(conditions[j] %in% condition_select && eeg_dat_J[i,j]){   
      
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Pain_Phenotype
      subject_group <- dat0$Pain_Phenotype[dat0$Subject == Subject[i]]
      
      if(length(subject_group) == 1){
        grid$group <- rep(subject_group, length(tt) * length(ff))
      }else{
        grid$group <- rep(NA, length(tt) * length(ff))
      }
      
      dat2plot <- rbind(dat2plot, grid)
    }
  }
}

dat2plot$group <- factor(dat2plot$group, levels = pain_levels) 
dat2plot$condition <- factor(dat2plot$condition, levels = condition_select)
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Back_32" = "Back 32", 
                                    "Back_256"= "Back 256", 
                                    "Hand_32" = "Hand 32",
                                    "Hand_256" = "Hand 256") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))


df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = mean(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

# Plot the results 
pp_mOFC_rh <- ggplot(df_summary, aes(x = time_sec, y = ff, fill = avg_value)) +
  geom_tile() + 
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0  
  ) + 
  facet_grid(group ~ condition) + 
  labs(x = "Time (s)",
       y = "Frequency (Hz)",
       fill = "Value") + 
  scale_x_continuous(breaks = desired_times, labels = breaks_labels) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major.x = element_blank(),   
        panel.grid.minor.x = element_blank()) + 
  ggtitle(paste(ROI_labels[region]))

plot(pp_mOFC_rh) 

grid.arrange(pp_mOFC_lh,  
             pp_mOFC_rh,  
             nrow = 2, ncol= 1, heights = c(1,1), widths = 1) 




# dlPFC_lh
region = 5
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    
    if(conditions[j] %in% condition_select && eeg_dat_J[i,j]){   
      
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Pain_Phenotype
      subject_group <- dat0$Pain_Phenotype[dat0$Subject == Subject[i]]
      
      if(length(subject_group) == 1){
        grid$group <- rep(subject_group, length(tt) * length(ff))
      }else{
        grid$group <- rep(NA, length(tt) * length(ff))
      }
      
      dat2plot <- rbind(dat2plot, grid)
    }
  }
}

dat2plot$group <- factor(dat2plot$group, levels = pain_levels) 
dat2plot$condition <- factor(dat2plot$condition, levels = condition_select)
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Back_32" = "Back 32", 
                                    "Back_256"= "Back 256", 
                                    "Hand_32" = "Hand 32",
                                    "Hand_256" = "Hand 256") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))


df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = mean(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

# Plot the results 
pp_dlPFC_lh <- ggplot(df_summary, aes(x = time_sec, y = ff, fill = avg_value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0  
  ) + 
  facet_grid(group ~ condition) + 
  labs(x = "Time (s)",
       y = "Frequency (Hz)", 
       fill = "Value") + 
  scale_x_continuous(breaks = desired_times, labels = breaks_labels) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major.x = element_blank(),   
        panel.grid.minor.x = element_blank()) + 
  ggtitle(paste(ROI_labels[region]))

#plot(pp_dlPFC_lh) 




###############
# dlPFC_rh
region = 11
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    
    if(conditions[j] %in% condition_select && eeg_dat_J[i,j]){   
      
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Pain_Phenotype
      subject_group <- dat0$Pain_Phenotype[dat0$Subject == Subject[i]]
      
      if(length(subject_group) == 1){
        grid$group <- rep(subject_group, length(tt) * length(ff))
      }else{
        grid$group <- rep(NA, length(tt) * length(ff))
      }
      
      dat2plot <- rbind(dat2plot, grid)
    }
  }
}

dat2plot$group <- factor(dat2plot$group, levels = pain_levels) 
dat2plot$condition <- factor(dat2plot$condition, levels = condition_select)
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Back_32" = "Back 32", 
                                    "Back_256"= "Back 256", 
                                    "Hand_32" = "Hand 32",
                                    "Hand_256" = "Hand 256") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))


df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = mean(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

# Plot the results 
pp_dlPFC_rh <- ggplot(df_summary, aes(x = time_sec, y = ff, fill = avg_value)) +
  geom_tile() +
  #geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +  # <-- added line 
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0  
  ) + 
  facet_grid(group ~ condition) + 
  labs(x = "Time (s)",
       y = "Frequency (Hz)",
       fill = "Value") + 
  scale_x_continuous(breaks = desired_times, labels = breaks_labels) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major.x = element_blank(),   
        panel.grid.minor.x = element_blank()) + 
  ggtitle(paste(ROI_labels[region]))

plot(pp_dlPFC_rh)  
grid.arrange(pp_dlPFC_lh,  
             pp_dlPFC_rh,  
             nrow = 2, ncol= 1, heights = c(1,1), widths = 1) 





######################################
######################################
######################################
## modeling 
######################################
######################################
######################################

### analyze EEG data using multivariate mixed effect models 
rm(list = ls())
library(splines)
library(pracma)
library(lattice)
library(gridExtra) 
library(mcompanion)
library(MASS)
library(mvtnorm)
library(tidyr)
library(ggplot2)
library(gtools) 
library(simdd)
library(VGAM)


setwd("/Users/hyung/Dropbox/LBP/0320_EEG")
source("fun_gen.R")
source("helper.R")
source("tensor_ssl.R")

dat_all <- readRDS("LBP_EEG_05062025.RData")
dat_phenotype <- read.csv("phenotype.csv")
n <- nrow(dat_phenotype)
colnames(dat_phenotype)
dat_phenotype$Pain_Phenotype <- ifelse(dat_phenotype$Pain_Phenotype == "Widespread Pain",
                                       "Widespread Hyperalgesia",
                                       dat_phenotype$Pain_Phenotype)


# Summary of pain intensity by pain phenotype and gender  
summary_pain <- dat_phenotype %>%
  group_by(Pain_Phenotype, Gender) %>%
  summarise(
    N = n(),
    Mean = mean(Pain_Intensity, na.rm = TRUE),
    SD = sd(Pain_Intensity, na.rm = TRUE),
    Median = median(Pain_Intensity, na.rm = TRUE),
    IQR = IQR(Pain_Intensity, na.rm = TRUE)
  )

print(summary_pain)

# convert treatment indicators to factor (just in case)
dat_phenotype$central_analgesic <- as.factor(dat_phenotype$central_analgesic)
dat_phenotype$peripheral_treatment <- as.factor(dat_phenotype$peripheral_treatment)
dat_phenotype$behavioral_modulation <- as.factor(dat_phenotype$behavioral_modulation)

# summarize as percentages 
treatment_summary <- dat_phenotype %>%
  group_by(Pain_Phenotype)  %>%
  summarise(
    n = n(),
    pct_central_analgesic = 100 * mean(as.numeric(central_analgesic)-1),
    pct_peripheral_treatment = 100 * mean(as.numeric(peripheral_treatment)-1),
    pct_behavioral_modulation = 100 * mean(as.numeric(behavioral_modulation)-1)) 

print(treatment_summary)


library(tableone) 
# Ensure variable types are correct
dat_phenotype$Gender <- as.factor(dat_phenotype$Gender)
dat_phenotype$central_analgesic <- as.factor(dat_phenotype$central_analgesic)
dat_phenotype$peripheral_treatment <- as.factor(dat_phenotype$peripheral_treatment)
dat_phenotype$behavioral_modulation <- as.factor(dat_phenotype$behavioral_modulation)
dat_phenotype$Pain_Phenotype <- factor(dat_phenotype$Pain_Phenotype,
                                       levels = c("No Pain", "Localized Pain", "Widespread Hyperalgesia"))

# Specify variable types
factorVars <- c("Gender") #, "central_analgesic", "peripheral_treatment", "behavioral_modulation")
contVars <- c("Pain_Intensity", "Age")
allVars <- c(factorVars, contVars)

# Create Table 1
table1 <- CreateTableOne(vars = allVars,
                         strata = "Pain_Phenotype",
                         data = dat_phenotype,
                         factorVars = factorVars)

# Print Table 1 with means and percentages
print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)

# Convert table to data frame
table1_df <- print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

# Write to CSV file
#write.csv(table1_df, file = "Table1_PainPhenotype.csv", row.names = TRUE)

#############
#############
#############
#############




# construct a design matrix  
XX <- matrix(0, n, 6) # fit model with pain category + gender + age + pain score
XX[,1] <- 1
XX[,2] <- dat_phenotype$Pain_Phenotype == "Localized Pain"
XX[,3] <- dat_phenotype$Pain_Phenotype == "Widespread Hyperalgesia" 
XX[,4] <- as.numeric(as.factor(dat_phenotype$Gender))-1  # 1 for male; 0 for female
XX[,5] <- dat_phenotype$Age/10   # age in decade 
XX[,6] <- dat_phenotype$Pain_Intensity
#XX[which(is.na(XX[,6])),6] <- median(XX[,6], na.rm = T)  

p <- ncol(XX)
p

conditions <- c("Back_32", "Back_128","Back_256", "Hand_32", "Back_128", "Hand_256")
condition_ids <- c(1,3,4,6) 
condition_select <- conditions[condition_ids]
JJ <- dat_all$eeg_dat_J[,condition_ids]
J <- ncol(JJ) # maximum number of conditions

roi_s <- c("rACC-lh", "dACC-lh", "S1-lh",
           "Ins-lh", "dlPFC-lh", "mOFC-lh",
           "rACC-rh", "dACC-rh", "S1-rh",
           "Ins-rh", "dlPFC-rh", "mOFC-rh")


########################################## 
# we will focus only on mOFC and dlPFC
# repeat the process for each: 
roi_idx <- which(roi_s == "mOFC-lh")    #1
#roi_idx <- which(roi_s == "mOFC-rh")   #2
#roi_idx <- which(roi_s == "dlPFC-lh")  #3
#roi_idx <- which(roi_s == "dlPFC-rh")  #4 
##########################################

# time and frequency points of interest 
tt <- 25:43  # [0,0.7] sec
ff <- 1:45
TT <- length(tt)
FF <- length(ff) 
time_seq <- seq(0, 0.7, length.out = length(tt))
desired_times <- c(0, 0.2, 0.4, 0.6)
breaks_idx <- sapply(desired_times, function(t) {
  tt[which.min(abs(time_seq - t))]
})


## note, YY is on 10*log10 scale (i.e., decibel scale)
# construct responses
YY <- array(NA, dim = c(n, length(condition_ids), TT *FF))
for(i in 1:n){
  for(j in 1:length(condition_ids))
    if(JJ[i,j]){
      YY[i,j,] <- c(dat_all$eeg_dat_Y[[i]][[condition_ids[j]]][tt,ff,roi_idx])
    }
}


# specify number of basis functions
K_T <- 5
K_F <- 10  
n_burn   <- 1500   
n_sample <- 500    
params   <- list(Sigma2_delta = rep(5, p),  aa_gamma = 5, bb_gamma = 1, aa_omega = 50, bb_omega = 1)

# used later to reconstruct the fitted 2D functions 
tmp_prep <- prepocess(TT, FF, K_T, K_F)
O_tilde  <- tmp_prep$O_tilde
bsMat_tt <- tmp_prep$bsMat_tt
bsMat_ff <- tmp_prep$bsMat_ff
C_t <-  tmp_prep$C_t
C_f <-  tmp_prep$C_f

seed <- 123 # for reproducibility 

# -- R is the maximum rank (we select rank based on a threshold)
# -- save_all = TRUE is used for saving burn-in and sampling iterations (otherwise it only saves sampling iterations)
# -- type_simple = TRUE runs the homogeneous variance case and is a little faster
# -- larger threshold will only select 1 rank

mOFC_lh_tensor_model = tensor_model <- tensor.decomp.ssl(seed, YY, XX, JJ, TT, FF, R = 2, K_T, K_F, n_sample, n_burn, params = params, inits  = NULL, save_all = TRUE, type_simple = TRUE, threshold =  0.0001)
#mOFC_rh_tensor_model = tensor_model <- tensor.decomp.ssl(seed, YY, XX, JJ, TT, FF, R = 2, K_T, K_F, n_sample, n_burn, params = params, inits  = NULL, save_all = TRUE, type_simple = TRUE, threshold =  0.0001)
#dlPFC_lh_tensor_model = tensor_model <- tensor.decomp.ssl(seed, YY, XX, JJ, TT, FF, R = 2, K_T, K_F, n_sample, n_burn, params = params, inits  = NULL, save_all = TRUE, type_simple = TRUE, threshold =  0.0001)
#dlPFC_rh_tensor_model = tensor_model <- tensor.decomp.ssl(seed, YY, XX, JJ, TT, FF, R = 2, K_T, K_F, n_sample, n_burn, params = params, inits  = NULL, save_all = TRUE, type_simple = TRUE, threshold =  0.0001)

dat2save <- list(mOFC_lh_tensor_model = mOFC_lh_tensor_model,
                 mOFC_rh_tensor_model = mOFC_rh_tensor_model, 
                 dlPFC_lh_tensor_model = dlPFC_lh_tensor_model,
                 dlPFC_rh_tensor_model = dlPFC_rh_tensor_model)
#saveRDS(file = "tensor_model_05062025.RData", dat2save)   
#dat2save <- readRDS("tensor_model_05062025.RData")   ## load the saved dataset 

tensor_model <- dat2save$mOFC_lh_tensor_model
#tensor_model <- dat2save$mOFC_rh_tensor_model
#tensor_model <- dat2save$dlPFC_lh_tensor_model
#tensor_model <- dat2save$dlPFC_rh_tensor_model

roi_idx <- which(roi_s == "mOFC-lh")    #1
#roi_idx <- which(roi_s == "mOFC-rh")   #2
#roi_idx <- which(roi_s == "dlPFC-lh")  #3
#roi_idx <- which(roi_s == "dlPFC-rh")  #4 


# index of the samples 
sample_idx <- (n_burn+1):(n_burn+n_sample)  

# selected rank index 
R_indices <- tensor_model$R_indices 
R <- length(tensor_model$R_indices)

#############################################################
## plots 
## estimated fixed effects of phenotype group 
subject_idx <- c( which(dat_phenotype$Pain_Phenotype == "Widespread Hyperalgesia")[1], 
                  which(dat_phenotype$Pain_Phenotype == "Localized Pain")[1], 
                  which(dat_phenotype$Pain_Phenotype == "No Pain")[1])

dat2plot <- NULL
for(j in 1:J){ # plot for each subject 
  grid <- expand.grid(tt = tt, ff = ff)
  grid$mat_w <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[1]], 2, mean)
  grid$mat_l <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[2]], 2, mean)
  grid$mat_c <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[3]], 2, mean)
  
  grid$condition <- rep(condition_select[j], TT*FF)
  
  df_long <- pivot_longer(grid, 
                          cols = c("mat_w", "mat_l", "mat_c"),  # Columns to gather
                          names_to = "mat_type")
  df_long$mat_type <- dplyr::recode(df_long$mat_type, "mat_w" = "Widespread Hyperalgesia", "mat_l" = "Localized Pain", "mat_c" = "No Pain")
  df_long$condition <- dplyr::recode(df_long$condition, 
                                     "Back_32" = "Back 32", 
                                     "Back_256" = "Back 256", 
                                     "Hand_32" = "Hand 32",
                                     "Hand_256" = "Hand 256")
  
  dat2plot  <- rbind(dat2plot,  df_long)
}


dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))
dat2plot$mat_type <- factor(dat2plot$mat_type, levels = c("Widespread Hyperalgesia","Localized Pain","No Pain"))

pp <- ggplot(dat2plot, aes(x = tt, y = ff, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Heatmap color scale
  facet_grid(mat_type~condition) +  # One panel per group
  scale_x_continuous(
    name = "Time (s)",
    breaks = breaks_idx,
    labels = desired_times
  ) + 
  labs(x = "Time",
       y = "Frequency (Hz)",
       fill = "dB change") +    theme(panel.background = element_blank(),
                                      plot.background = element_blank(),
                                      panel.grid.major.x = element_blank(),   
                                      panel.grid.minor.x = element_blank(),
                                      strip.text = element_text(size = 12, face = "bold")
       ) 

pp 


## plot the Sigmas (variance components of the random effects)
df <- data.frame(
  value = sqrt(c(tensor_model$Sigma2_s$Sigma2_gamma_s[sample_idx],
                 tensor_model$Sigma2_s$Sigma2_omega_s[sample_idx,1],
                 tensor_model$Sigma2_s$Sigma2_epsilon_s[sample_idx])),
  group = factor(rep(c("1", "2", "3"), each = length(sample_idx)))
)

ggplot(df, aes(x = group, y = value)) +
  geom_boxplot() +
  scale_x_discrete(
    labels = c( 
      "1" = expression(sigma[gamma]),
      "2" = expression(sigma[omega]),
      "3" = expression(sigma[epsilon])
    )) +  labs(x = NULL, y = "Standard Deviation") +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 16)) +
  theme(strip.text = element_text(size = 14))  # Adjust facet title size 




## plot the tensor basis patches
dat2plot <- NULL
for(r in R_indices){
  
  grid <- expand.grid(tt = tt, ff = ff)
  val_tt_est <-  (bsMat_tt%*%C_t)%*% apply(tensor_model$u_s[sample_idx,r, ], 2, mean)
  #val_ff_est <-  (bsMat_ff%*%C_f)%*% apply(tensor_model$v_s[sample_idx,r, ], 2, mean)
  val_ff_est <- - (bsMat_ff%*%C_f)%*% apply(tensor_model$v_s[sample_idx,r, ], 2, mean)
  
  grid$CP_est <- kronecker(val_ff_est, val_tt_est)
  grid$rank <- rep(paste("Component ", which(R_indices == r), sep = ""), TT*FF)
  
  df_long <- pivot_longer(grid, 
                          cols = c("CP_est"),  # Columns to gather
                          names_to = "z") 
  
  dat2plot  <- rbind(dat2plot,  df_long)
}



dat2plot$rank[dat2plot$rank == "Component 1"] <- paste(roi_s[roi_idx], "Component") #"Component"
## can add more depending on how many ended up being selected 
p1 <- ggplot(#dat2plot[dat2plot$rank == "Component",], 
  dat2plot[dat2plot$rank == paste(roi_s[roi_idx], "Component"),], 
  aes(x = tt, y = ff, fill = value)) +
  geom_tile() + 
  scale_fill_viridis_c() +  
  facet_wrap(~rank,nrow = 1) + 
  scale_x_continuous(
    name = "Time (s)",
    breaks = breaks_idx,
    labels = desired_times
  ) + 
  labs(x = "Time",
       y = "Frequency (Hz)",
       fill = "dB change") +    theme(panel.background = element_blank(),
                                      plot.background = element_blank(),
                                      panel.grid.major.x = element_blank(),   
                                      panel.grid.minor.x = element_blank()
       ) +theme(legend.position  = "left") +  theme(strip.text = element_text(size = 12, face = "bold"))
p1


## principal functions and credible intervals 
pt_list <- list()
pf_list <- list() 

for(r in R_indices){
  
  # there's a sign ambiguity, we can simply change the sign of the marginal components appropriately to improve interpretation 
  data_t <- data.frame(
    x = tt,
    mean = apply(bsMat_tt%*%C_t%*%t(tensor_model$u_s[sample_idx,r, ]), 1, mean),
    lower_95 = apply(bsMat_tt%*%C_t%*%t(tensor_model$u_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.025)}),
    upper_95 = apply(bsMat_tt%*%C_t%*%t(tensor_model$u_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.975)})
  )
  # data_f <- data.frame(
  #   x = ff,
  #   mean = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, mean),
  #   lower_95 = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.025)}),
  #   upper_95 = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.975)})
  # )
  data_f <- data.frame(
    x = ff,
    mean = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, mean),
    lower_95 = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.025)}),
    upper_95 = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.975)})
  )
  
  rr <- which(R_indices == r)
  
  if(rr==1){
    pt_list[[rr]] = ggplot(data_t, aes(x = x)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "gray60", alpha = 0.5) +  # 95% CI
      geom_line(aes(y = mean), color = "black", linewidth= 0.6) +  # Posterior mean
      labs(y = expression(phi[1]^"*"~"(t)"),
           x = "Time") +scale_x_continuous(
             name = "Time (s)",
             breaks = breaks_idx,
             labels = desired_times
           ) 
    
    pf_list[[rr]] = ggplot(data_f, aes(x = x)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "gray60", alpha = 0.5) +  # 95% CI
      geom_line(aes(y = mean), color = "black", linewidth = 0.6) +  # Posterior mean
      labs(y = expression(psi[1]^"*"~"(f)"),
           x = "Frequency (Hz)")
  }
  
  if(rr==2){
    pt_list[[rr]] = ggplot(data_t, aes(x = x)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "gray60", alpha = 0.5) +  # 95% CI
      geom_line(aes(y = mean), color = "black", linewidth= 0.6) +  # Posterior mean
      labs(y = expression(phi[2]^"*"~"(t)"),
           x = "Time")+scale_x_continuous(
             name = "Time (s)",
             breaks = breaks_idx,
             labels = desired_times
           ) 
    pf_list[[rr]] = ggplot(data_f, aes(x = x)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "gray60", alpha = 0.5) +  # 95% CI
      geom_line(aes(y = mean), color = "black", linewidth = 0.6) +  # Posterior mean
      labs(y = expression(psi[2]^"*"~"(f)"),
           x = "Frequency (Hz)")
  }
  
}

# adjust the plotting format 
for(r in R_indices){
  rr <- which(R_indices == r)
  pf_list[[rr]] =  pf_list[[rr]] + theme_minimal() + theme(
    panel.grid = element_blank(),  # Remove all grid lines
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text = element_text(size = 12, face = "bold"),  # Adjust text size
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12) 
  )   + theme(plot.margin = margin(40, 10, 10, 10))
  pt_list[[rr]] =  pt_list[[rr]] + theme_minimal() + theme(
    panel.grid = element_blank(),  # Remove all grid lines
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text = element_text(size = 12),  # Adjust text size
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12) 
  )  +  theme(plot.margin = margin(40, 10, 10, 10))
}

plots = grid.arrange(p1,  pt_list[[1]], pf_list[[1]],
                     #p2, pt_list[[2]], pf_list[[2]], 
                     #nrow = 2, ncol= 3, heights = c(1,1), widths = c(1.3,1,1))
                     nrow = 1, ncol= 3, heights = 1, widths = c(1.3,1,1))
plots


## plot the 95% credible intervals for delta coefficients 
pp_list <- list()

for(r in R_indices){
  dat2plot <- NULL
  for(j in 1:J){
    # there's a sign ambiguity, we can simply change the sign of the marginal components appropriately to improve interpretation 
    #tmp <-  data.frame(cbind(tensor_model$delta_s[sample_idx,j,r, 1], tensor_model$delta_s[sample_idx,j,r, 1]+ tensor_model$delta_s[sample_idx,j,r, 2], tensor_model$delta_s[sample_idx,j,r, 1]+ tensor_model$delta_s[sample_idx,j,r, 3] ))
    tmp <-  - data.frame(cbind(tensor_model$delta_s[sample_idx,j,r, 1], tensor_model$delta_s[sample_idx,j,r, 1]+ tensor_model$delta_s[sample_idx,j,r, 2], tensor_model$delta_s[sample_idx,j,r, 1]+ tensor_model$delta_s[sample_idx,j,r, 3] ))
    
    tmp2 = tmp[,3] - tmp[,2] ## difference in the effect between widespread vs. localized pain
    lower2 <- quantile(tmp2, 0.05)
    upper2 <- quantile(tmp2, 0.95)
    mid2 <- mean(tmp2)
    print(c(mid2, lower2, upper2))
    
    lower <- apply(tmp, 2, function(u){quantile(u, 0.05)})
    upper <- apply(tmp, 2, function(u){quantile(u, 0.95)})
    mid <- apply(tmp, 2, mean) #apply(tmp, 2, median)
    tmp <- data.frame(cbind(lower, upper,mid,  c("No Pain", "Localized Pain", "Widespread Hyperalgesia"), rep(condition_select[j], 3)))
    dat2plot <- rbind(dat2plot, tmp)
  }
  colnames(dat2plot) <- c("lower","upper","mid", "group", "condition")
  
  dat2plot$group <- factor(dat2plot$group, c("Widespread Hyperalgesia", "Localized Pain", "No Pain"))
  dat2plot$lower  <- as.numeric(dat2plot$lower) 
  dat2plot$upper  <- as.numeric(dat2plot$upper) 
  dat2plot$mid <- as.numeric(dat2plot$mid) 
  
  dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                      "Back_32" = "Back 32", 
                                      "Back_256" = "Back 256", 
                                      "Hand_32" = "Hand 32",
                                      "Hand_256" = "Hand 256")  
  dat2plot$condition <- factor(dat2plot$condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))
  
  rr <- which(R_indices == r)
  
  pp_list[[rr]] <- ggplot(dat2plot, aes(y = condition, x = mid, color = group)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.5)) +
    theme_minimal() +
    labs(
      x = "Group-specific effect",
      y = "Condition",
      color = ""  # Legend title for the group color
    ) +
    theme(
      legend.position = "top",  # or "top", or inside like c(0.8, 0.8)
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 14, face = "bold")#,
      #plot.margin = margin(40, 10, 10, 10)
    )
}


grid.arrange(p1 , #p2, 
             pp_list[[1]],
             #pp_list[[2]],
             #nrow = 2, ncol= 2, heights = c(1,1), widths = c(1,1)) 
             nrow = 1, ncol= 2, heights = 1, widths = c(1,1)) 
## we can repeat this for the other ROIs 




contrast_vec <- c(0,1,0,1)
names(contrast_vec) <- condition_select  # Ensure correct order

pp_contrast_ci_list <- list()

for (r in R_indices) {
  dat2plot <- NULL
  
  for (g in 1:3) {
    contrast_samples <- sapply(1:length(sample_idx), function(s) {
      sum(contrast_vec * sapply(1:J, function(j) {
        if (g == 1) {
          return(tensor_model$delta_s[sample_idx[s], j, r, 1])
        } else {
          return(tensor_model$delta_s[sample_idx[s], j, r, 1] +
                   tensor_model$delta_s[sample_idx[s], j, r, g])
        }
      }))
    })
    
    lower <- quantile(contrast_samples, 0.05)
    upper <- quantile(contrast_samples, 0.95)
    mid <- median(contrast_samples)
    
    tmp <- data.frame(lower = lower,
                      upper = upper,
                      mid = mid,
                      group = c("No Pain", "Localized Pain", "Widespread Hyperalgesia")[g])
    
    dat2plot <- rbind(dat2plot, tmp)
  }
  
  dat2plot$group <- factor(dat2plot$group, c("Widespread Hyperalgesia", "Localized Pain", "No Pain"))
  rr <- which(R_indices == r)
  
  pp_contrast_ci_list[[rr]] <- ggplot(dat2plot, aes(y = group, x = mid, color = group)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,
                   position = position_dodge(width = 0.5)) +
    theme_minimal() +
    labs(title = #paste("Contrast-weighted effect (Component =", rr, ")"),
           paste("Contrast-weighted effect"),
         x = "score",#"Posterior median (95% CI)",
         y = NULL) +
    theme(
      legend.position = c(1.05, 1.2), 
      legend.title = element_blank(),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.text = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(40, 10, 10, 10)
    )
}

# Example: plot first component horizontally with consistent legend
plot(pp_contrast_ci_list[[1]])







#############
#############
# projection-based analysis 
#############
############# 

base_patterns <- NULL

for(r in R_indices){
  rr <- which(R_indices == r) 
  val_tt_est <-  (bsMat_tt%*%C_t)%*% apply(tensor_model$u_s[sample_idx,r, ], 2, mean)
  val_ff_est <-  -(bsMat_ff%*%C_f)%*% apply(tensor_model$v_s[sample_idx,r, ], 2, mean)
  base_patterns <- rbind(base_patterns, c(kronecker(val_ff_est, val_tt_est)))
}


#############################
#group membership stored separately. 
group_labels <- factor(dat_phenotype$Pain_Phenotype, levels = c("Widespread Hyperalgesia", "Localized Pain", "No Pain"))  # Factor with levels: "Widespread Pain", "Localized Pain", "No Pain"

# Initialize
res_all <- array(NA, dim = c(n, J, R))  # n = number of subjects, J = number of conditions, R = number of components

# Compute inner products
for(i in 1:n){
  for(j in 1:J){
    if(sum(is.na(YY[i,j, ])) == 0){
      for(r in 1:R){
        res_all[i,j,r] <- sum(YY[i,j, ] * base_patterns[r, ])
      }
    }
  }
}

# Create plots
pp_ip_list <- list()
dat2plot <- list()

for(r in R_indices){
  rr <- which(R_indices == r)
  
  # Flatten the res_all array to a data frame
  temp_df <- data.frame(
    Value = as.vector(res_all[,,rr]),
    Condition = rep(condition_select, each = n),
    Subject = rep(1:n, times = J),
    Group = rep(group_labels, times = J)
  )
  
  # Clean and prepare for plotting
  long_data_clean <- temp_df %>% filter(!is.na(Value))
  #long_data_clean$Condition <- factor(long_data_clean$Condition, levels = condition_select)
  
  long_data_clean$Condition <- dplyr::recode(long_data_clean$Condition, 
                                             "Back_32" = "Back 32", 
                                             "Back_256" = "Back 256", 
                                             "Hand_32" = "Hand 32",
                                             "Hand_256" = "Hand 256")  
  long_data_clean$Condition <- factor(long_data_clean$Condition, levels = c("Back 32", "Back 256", "Hand 32", "Hand 256"))
  
  
  # Create boxplot
  pp_ip_list[[rr]] <- ggplot(long_data_clean, aes(x = Condition, y = Value, fill = Group)) +
    geom_boxplot() +
    theme_minimal() +
    labs(#title = paste("Component = ", rr),
      x = "Condition",
      y = "Inner Product") +
    coord_cartesian(ylim = c(-30, 40)) + 
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 14))
}

# Display the first plot (only one if R = 1)
plot(pp_ip_list[[1]])



############ 4/14/2025 
# custom_projection  
# Time and frequency grid info
TT <- length(unique(tt))  # Number of time points in trimmed data
#19
FF <- length(unique(ff))  # Number of frequency points
#45

ROI_labels <- c("rACC-lh", "dACC-lh", "S1-lh",
                "Ins-lh", "dlPFC-lh", "mOFC-lh",
                "rACC-rh", "dACC-rh", "S1-rh",
                "Ins-rh", "dlPFC-rh", "mOFC-rh")

# Construct custom contrast pattern
custom_pattern <- matrix(0, nrow = TT, ncol = FF)

# Activate tt = 37:43 and ff = 1:12
active_tt <- 35:43   # time-window: (0.4,0.7)
active_ff_pos <- 2:7  
active_ff_neg <- 8:12
# Map to 1-based indices relative to tt and ff
tt_map <- match(active_tt, tt)
ff_pos_map <- match(active_ff_pos, ff)
ff_neg_map <- match(active_ff_neg, ff)

# Fill in the weights
custom_pattern[tt_map, ff_pos_map] <- 1
custom_pattern[tt_map, ff_neg_map] <- -1

custom_pattern

# Vectorize the pattern
custom_projection_vec <- as.vector(custom_pattern)

# Initialize res_all with only one projection (e.g., r = 1)
R <- 1  # Since we're focusing on a single custom projection
res_all <- array(NA, dim = c(n, J, R))  # [subject, condition, projection index]

# Compute projection scores as inner products
for (i in 1:n) {
  for (j in 1:J) {
    if (sum(is.na(YY[i, j, ])) == 0) {
      #res_all[i, j, 1] <- sum(YY[i, j, ] * custom_projection_vec)
      tmp <- dat_all$eeg_dat_Y[[i]][[condition_ids[j]]]
      res_all[i, j, 1]  <- 
        sum(tmp[tt, ff, which(ROI_labels == "mOFC-lh")] * custom_pattern) + 
        sum(tmp[tt, ff, which(ROI_labels == "mOFC-rh")] * custom_pattern)
    }
  }
}



########## contrast-based projection 

contrast_vec <- c(0, 1, 1, 1) 
names(contrast_vec) <- c("Back_32", "Back_256", "Hand_32", "Hand_256")

# Ensure levels of Condition match the contrast vector
condition_levels <- c("Back_32", "Back_256", "Hand_32", "Hand_256")
contrast_vec <- contrast_vec[condition_levels]  # Ensure correct order

# Initialize storage
contrast_scores <- matrix(NA, nrow = n, ncol = R)  # [subject, component]

# Apply contrast for each subject and component
for (r in 1:R) {
  for (i in 1:n) {
    subject_vals <- res_all[i, , r]
    if (!any(is.na(subject_vals))) {
      contrast_scores[i, r] <- sum(subject_vals * contrast_vec)
    }
  }
}

R_indices = 1
# Create plots
pp_contrast_list <- list()
for (r in R_indices) {
  rr <- which(R_indices == r)
  
  df_contrast <- data.frame(
    Subject = 1:n,
    ContrastValue = contrast_scores[, rr],
    Group = group_labels
  )
  
  df_contrast <- df_contrast %>% filter(!is.na(ContrastValue))
  
  pp_contrast_list[[rr]] <- ggplot(df_contrast, aes(x = Group, y = ContrastValue, fill = Group)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Inner product (Component =", rr, ")"),
         x = "Group",
         y = "Contrast Score")+ 
    coord_cartesian(ylim = c(-400, 600)) + 
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 14))
}

# Example: Display first component
plot(pp_contrast_list[[1]])


# end of the code