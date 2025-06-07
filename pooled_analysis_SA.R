# 6/3/2025
# -- set directories 
#setwd("/Users/hyung/Dropbox/LBP/github_mOFC_signature")

rm(list = ls()) 
source("fun_gen.R")
source("helper.R")
source("tensor_ssl.R")


LBP_dat2save <- readRDS("LBP_EEG_05062025.RData")   ## load the saved dataset  
CP_dat2save <- readRDS("CP_EEG_06032025.RData") 
LBP_eeg_dat_Y <- LBP_dat2save$eeg_dat_Y
length(LBP_eeg_dat_Y)
is(LBP_eeg_dat_Y)
LBP_eeg_dat_J <- LBP_dat2save$eeg_dat_J
dim(LBP_eeg_dat_J)
CP_eeg_dat_Y <- CP_dat2save$eeg_dat_Y
length(CP_eeg_dat_Y)
CP_eeg_dat_J <- CP_dat2save$eeg_dat_J
dim(CP_eeg_dat_J)


# -- ROIs 
ROI_labels <- c("rACC-lh", "dACC-lh", "S1-lh",
                "Ins-lh", "dlPFC-lh", "mOFC-lh",
                "rACC-rh", "dACC-rh", "S1-rh",
                "Ins-rh", "dlPFC-rh", "mOFC-rh")
# -- CP conditions 
CP_conditions <- c("Abdomen_32", "Abdomen_512","Forearm_32", "Forearm_512")

# -- All conditions and subset to condition_select
LBP_conditions <- c("Back_32", "Back_128", "Back_256", "Hand_32", "Back_128", "Hand_256")
LBP_condition_ids <- c(1, 3, 4, 6)
LBP_condition_select <- LBP_conditions[LBP_condition_ids]  # Selected subset
LBP_condition_select <- unique(LBP_condition_select)   # Just to be safe, remove any duplicates

LBP_dat_phenotype <- read.csv("phenotype.csv")
n <- nrow(LBP_dat_phenotype)
colnames(LBP_dat_phenotype)
LBP_dat_phenotype$Phenotype <- ifelse(LBP_dat_phenotype$Pain_Phenotype == "Widespread Pain",
                                      "Widespread Hyperalgesia",
                                      LBP_dat_phenotype$Pain_Phenotype)

LBP_dat_phenotype$Phenotype <- factor(LBP_dat_phenotype$Phenotype, 
                                      levels = c("No Pain", "Localized Pain", "Widespread Hyperalgesia"))


CP_dat_phenotype <- read_excel("25_05_08Phenotype_CP.xlsx", sheet = 1) 
colnames(CP_dat_phenotype) 

# Ensure variable types are correct
CP_dat_phenotype$Sex <- as.factor(CP_dat_phenotype$Sex)
CP_dat_phenotype$Phenotype <- factor(CP_dat_phenotype$Phenotype, 
                                     levels = c("No Hyperalgesia", "Segmental Hyperalgesia", "Widespread Hyperalgesia"))

colnames(CP_dat_phenotype)
colnames(LBP_dat_phenotype)
CP_dat_phenotype$Phenotype
LBP_dat_phenotype$Phenotype



# Define the harmonized condition names
Condition <- c("Affected Low", "Affected High", "Unaffected Low", "Unaffected High")
names(Condition) <- c("Back_32", "Back_256", "Hand_32", "Hand_256")  # For LBP
names(Condition) <- c("Abdomen_32", "Abdomen_512", "Forearm_32", "Forearm_512")  # For CP

# === Merge EEG data ===
# LBP: subset to selected conditions
LBP_eeg_dat_Y_subset <- lapply(LBP_eeg_dat_Y, function(x) x[LBP_condition_ids])
LBP_eeg_dat_Y_renamed <- lapply(LBP_eeg_dat_Y_subset, function(lst) {
  names(lst) <- Condition
  return(lst)
})

# CP: already in correct order
names(CP_eeg_dat_Y) <- NULL  # reset any prior naming
CP_eeg_dat_Y_renamed <- lapply(CP_eeg_dat_Y, function(lst) {
  names(lst) <- Condition
  return(lst)
})

# Combine EEG data
combined_eeg_dat_Y <- c(LBP_eeg_dat_Y_renamed, CP_eeg_dat_Y_renamed)
dim(combined_eeg_dat_Y[[85]][[4]]) 


# === Merge EEG availability (J) matrices ===
# Subset LBP J to selected conditions
LBP_eeg_dat_J_subset <- LBP_eeg_dat_J[, LBP_condition_ids]
colnames(LBP_eeg_dat_J_subset) <- Condition
colnames(CP_eeg_dat_J) <- Condition

# Combine J matrices
combined_eeg_dat_J <- rbind(LBP_eeg_dat_J_subset, CP_eeg_dat_J)

# === Merge phenotype data ===
# Standardize columns
LBP_phenotype <- data.frame(
  Subject = LBP_dat_phenotype$Subject,
  Sex = LBP_dat_phenotype$Gender,
  Age = LBP_dat_phenotype$Age,
  Phenotype = LBP_dat_phenotype$Phenotype,
  Pain_Intensity = LBP_dat_phenotype$Pain_Intensity, 
  Cohort = "LBP"
)

CP_phenotype <- data.frame(
  Subject = CP_dat_phenotype$record_id,
  Sex = CP_dat_phenotype$Sex,
  Age = CP_dat_phenotype$Age,
  Phenotype = CP_dat_phenotype$Phenotype,
  Pain_Intensity = CP_dat_phenotype$BPI_Baseline, 
  Cohort = "CP"
)

CP_phenotype$Subject <- paste0("CP", CP_phenotype$Subject)

LBP_phenotype$Subject
CP_phenotype$Subject

# Combine phenotype
combined_phenotype <- rbind(LBP_phenotype, CP_phenotype)
dim(combined_phenotype )
head(combined_phenotype)
tail(combined_phenotype)
# Re-factor Phenotype levels if needed
combined_phenotype$Phenotype <- factor(combined_phenotype$Phenotype, 
                                       levels = c("No Pain", "Localized Pain", "Widespread Hyperalgesia", 
                                                  "No Hyperalgesia", "Segmental Hyperalgesia"))

table(combined_phenotype$Phenotype)

# Recode the phenotype levels
combined_phenotype$Phenotype <- as.character(combined_phenotype$Phenotype)

combined_phenotype$Phenotype[combined_phenotype$Phenotype %in% 
                               c("Localized Pain", "No Hyperalgesia", "Segmental Hyperalgesia")] <-"Localized Pain"
 
n <- nrow(combined_phenotype)

# construct a design matrix  
XX <- matrix(0, n, 7) # fit model with pain category + gender + age + pain score + cohort
XX[,1] <- 1
XX[,2] <- combined_phenotype$Phenotype == "Localized Pain"
XX[,3] <- combined_phenotype$Phenotype == "Widespread Hyperalgesia" 
XX[,4] <- as.numeric(as.factor(combined_phenotype$Sex))-1  # 1 for male; 0 for female
XX[,5] <- combined_phenotype$Age/10   # age in decade 
XX[,6] <- combined_phenotype$Pain_Intensity 
XX[,7] <- as.numeric(as.factor(combined_phenotype$Cohort)) -1 

p <- ncol(XX)
p 
JJ <- combined_eeg_dat_J
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
YY <- array(NA, dim = c(n, J, TT *FF))
for(i in 1:n){
  for(j in 1:J)
    if(JJ[i,j]){
      YY[i,j,] <- c(combined_eeg_dat_Y[[i]][[j]][tt,ff,roi_idx])
    }
}


# specify number of basis functions
K_T <- 5
K_F <- 10  
n_burn   <-  1500   
n_sample <-  500    
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
#saveRDS(file = "pooled_model_06032025.RData", dat2save)   
#dat2save <- readRDS("pooled_model_06032025.RData")   ## load the saved dataset 
#dat2save <- readRDS("combined_model_05232025.RData")   ## load the saved dataset 

tensor_model <- dat2save$mOFC_lh_tensor_model
#tensor_model <- dat2save$mOFC_rh_tensor_model
#tensor_model <- dat2save$dlPFC_lh_tensor_model
#tensor_model <- dat2save$dlPFC_rh_tensor_model

roi_idx <- which(roi_s == "mOFC-lh")    #1
#roi_idx <- which(roi_s == "mOFC-rh")   #2
#roi_idx <- which(roi_s == "dlPFC-lh")  #3
#roi_idx <- which(roi_s == "dlPFC-rh")  #4 




################################
################################
################################
# index of the samples 
sample_idx <- (n_burn+1):(n_burn+n_sample)  

# selected rank index 
R_indices <- tensor_model$R_indices 
R <- length(tensor_model$R_indices)


Condition <- c("Affected Low", "Affected High", "Unaffected Low", "Unaffected High")

#############################################################
## plots 
## estimated fixed effects of phenotype group 
subject_idx <- c( which(combined_phenotype$Phenotype == "Widespread Hyperalgesia")[1], 
                  which(combined_phenotype$Phenotype == "Localized Pain")[1], 
                  which(combined_phenotype$Phenotype == "No Pain")[1])

dat2plot <- NULL
for(j in 1:J){ # plot for each subject 
  grid <- expand.grid(tt = tt, ff = ff)
  grid$mat_w <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[1]], 2, mean)
  grid$mat_l <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[2]], 2, mean)
  grid$mat_c <-  O_tilde%*%apply(tensor_model$alpha_s[sample_idx,,j,subject_idx[3]], 2, mean)
  
  grid$condition <- rep(Condition[j], TT*FF)
  
  df_long <- pivot_longer(grid, 
                          cols = c("mat_w", "mat_l", "mat_c"),  # Columns to gather
                          names_to = "mat_type")
  df_long$mat_type <- dplyr::recode(df_long$mat_type, "mat_w" = "Widespread Hyperalgesia", "mat_l" = "Localized Pain", "mat_c" = "No Pain")
  dat2plot  <- rbind(dat2plot,  df_long)
}


dat2plot$condition <- factor(dat2plot$condition, levels = c("Affected Low", "Affected High", "Unaffected Low", "Unaffected High"))
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


## principal functions and 95% credible intervals 
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
  data_f <- data.frame(
    x = ff,
    mean = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, mean),
    lower_95 = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.025)}),
    upper_95 = apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.975)})
  )
  # data_f <- data.frame(
  #   x = ff,
  #   mean = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, mean),
  #   lower_95 = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.025)}),
  #   upper_95 = - apply(bsMat_ff%*%C_f%*%t(tensor_model$v_s[sample_idx,r, ]), 1, function(u){quantile(u, 0.975)})
  # )
  
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
    tmp <- data.frame(cbind(lower, upper,mid,  c("No Pain", "Localized Pain", "Widespread Hyperalgesia"), rep(Condition[j], 3)))
    dat2plot <- rbind(dat2plot, tmp)
  }
  colnames(dat2plot) <- c("lower","upper","mid", "group", "condition")
  
  dat2plot$group <- factor(dat2plot$group, c("Widespread Hyperalgesia", "Localized Pain", "No Pain"))
  dat2plot$lower  <- as.numeric(dat2plot$lower) 
  dat2plot$upper  <- as.numeric(dat2plot$upper) 
  dat2plot$mid <- as.numeric(dat2plot$mid) 
   
  dat2plot$condition <- factor(dat2plot$condition, levels = c("Affected Low", "Affected High", "Unaffected Low", "Unaffected High")
  )
  
  rr <- which(R_indices == r)
  
  pp_list[[rr]] <- ggplot(dat2plot, aes(y = condition, x = mid, color = group)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.5)) +
    theme_minimal() +
    scale_x_continuous(breaks = c(-5,0, 5,10), 
                       #limits = c(-5, NA)
                       ) +  # Add this line
    labs(
      x = "Group-specific effect",
      y = "Condition",
      color = ""
    ) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 14, face = "bold")
    )
  
}

p1_lh <- p1
pp_list_lh <- pp_list 
#p1_rh <- p1
#pp_list_rh <- pp_list 

grid.arrange(p1_lh ,  
             pp_list_lh[[1]], 
             p1_rh,  
             pp_list_rh[[1]], 
             #pp_list[[2]],
             nrow = 2, ncol= 2, heights = c(1,1), widths = c(1,1)) 

## we can repeat this for the other ROIs 







contrast_vec <- c(0,1,1,1) 
names(contrast_vec) <- Condition #condition_select  # Ensure correct order

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


 

#########################
###### projection   #####
###### 6/3/2025    #####
#########################

Condition <- c("Affected Low", "Affected High", "Unaffected Low", "Unaffected High")
J <- length(Condition)
colnames(combined_phenotype) 
n <- nrow(combined_phenotype)  
Subject <- combined_phenotype$Subject

# Define time and frequency indices 
tt <- 25:43  # [0,0.7] sec
ff <- 1:45   # full frequency index

active_tt <- 35:43   # time-window: (0.4,0.7)
active_ff_pos <- 2:7  
active_ff_neg <- 8:12

# Map to index relative to tt and ff
tt_map <- match(active_tt, tt)
ff_pos_map <- match(active_ff_pos, ff)
ff_neg_map <- match(active_ff_neg, ff)

TT <- length(tt)
FF <- length(ff)

# Construct contrast weight matrix (TT x FF)
custom_pattern <- matrix(0, nrow = TT, ncol = FF)
custom_pattern[tt_map, ff_pos_map] <-  1
custom_pattern[tt_map, ff_neg_map] <- -1

# Vectorize projection vector
custom_projection_vec <- as.vector(custom_pattern)

# Initialize result storage


res_proj = res_proj_lh = res_proj_rh <- matrix(NA, nrow = n, ncol = J) 


# Loop through subjects and conditions
for (i in 1:n) {
  for (j in 1:J) {
    if (combined_eeg_dat_J[i,j]) {
      res_proj_lh[i,j] <- sum(combined_eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-lh")] * custom_pattern)
      res_proj_rh[i,j] <- sum(combined_eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-rh")] * custom_pattern)
      res_proj[i,j] <- sum(combined_eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-lh")] * custom_pattern) + 
        sum(combined_eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-rh")] * custom_pattern)
    }
  }
}


# -- Initialize output data
proj_plot_data <- NULL

for (i in 1:n) {
  subj_id <- combined_phenotype$Subject[i]
  subj_group <- combined_phenotype$Phenotype[combined_phenotype$Subject == subj_id]
  
  for (j in 1:J) {
    if (combined_eeg_dat_J[i,j]) {
      score <- res_proj[i,j]
      score_lh <- res_proj_lh[i,j]
      score_rh <- res_proj_rh[i,j]
      proj_plot_data <- rbind(proj_plot_data, data.frame(
        Subject = subj_id,
        Group = subj_group,
        Condition = Condition[j],
        ProjectionScore = score,
        ProjectionScore_lh = score_lh,
        ProjectionScore_rh = score_rh
      ))
    }
  }
}



region = 6
# Clean levels
proj_plot_data$Group <- factor(proj_plot_data$Group,
                               levels = c("Widespread Hyperalgesia", "Localized Pain", "No Pain"))
proj_plot_data$Condition <- factor(proj_plot_data$Condition, levels = Condition)
colnames(proj_plot_data)


# Plot: 3 groups side-by-side within each of 4 conditions
library(ggplot2)

gg_proj_plot <- ggplot(proj_plot_data, 
                       aes(x = Condition, 
                           y = ProjectionScore, #y = ProjectionScore_lh
                           fill = Group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  labs(title = paste("Projection Score Comparison:", ROI_labels[region]),
       x = "Condition",
       y = "Projection Score") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

print(gg_proj_plot)



# ---- Define contrast weights over conditions 
contrast_vec <- c(0,1,1, 1)
condition_levels = names(contrast_vec) <- Condition 
# Ensure conditions match
contrast_vec <- contrast_vec[condition_levels]

#res_proj_lh[i,j]

# ---- Compute contrast-based projection score per subject
projection_scores = projection_scores_lh = projection_scores_rh <- rep(NA, n)

for (i in 1:n) {
  projection_scores[i]    <- sum(res_proj[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_lh[i] <- sum(res_proj_lh[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_rh[i] <- sum(res_proj_rh[i, ] * contrast_vec, na.rm = TRUE)
}


colnames(combined_phenotype)

df_projection <- data.frame(
  Subject = combined_phenotype$Subject,
  Group = combined_phenotype$Phenotype,
  Pain_Intensity = combined_phenotype$Pain_Intensity, 
  Age = combined_phenotype$Age, 
  Sex = combined_phenotype$Sex, 
  ProjectionScore_lh = projection_scores_lh,
  ProjectionScore_rh = projection_scores_rh,
  ProjectionScore = projection_scores
)


df_projection <- df_projection %>%
  filter(!is.na(ProjectionScore)) %>%
  mutate(Group = factor(Group,
                        levels = c("Widespread Hyperalgesia", "Localized Pain", "No Pain")))


#############
p_projection <- ggplot(df_projection, 
                       aes(x = Group, 
                           y = ProjectionScore, 
                           fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste("EEG Marker Projection Score"),
       x = "Pain Phenotype",
       y = "Projection Score") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )

print(p_projection)


#############
p_projection_lh <- ggplot(df_projection, 
                          aes(x = Group, 
                              y = ProjectionScore_lh, 
                              fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(#title = paste("EEG Marker Projection Score:", "mOFC-lh"),
    title = paste("Projection for pooled data:", "mOFC-lh"),
    x = "Pain Phenotype",
    y = "Projection Score") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
print(p_projection_lh)

p_projection_rh <- ggplot(df_projection, 
                          aes(x = Group, 
                              y = ProjectionScore_rh, 
                              fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste("Projection for pooled data:", "mOFC-rh"),
       x = "Pain Phenotype",
       y = "Projection Score") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )

# Show plot
print(p_projection_rh)


grid.arrange(p_projection_lh,  
             p_projection_rh, 
             nrow = 1, ncol= 2, heights = 1, widths = c(1,1)) 


############### Mean and SEM
# ---- Summarize: compute mean and SE for LH and RH
df_summary_lh <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore_lh, na.rm = TRUE),
    SE = sd(ProjectionScore_lh, na.rm = TRUE) / sqrt(n())
  )

df_summary_rh <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore_rh, na.rm = TRUE),
    SE = sd(ProjectionScore_rh, na.rm = TRUE) / sqrt(n())
  )

# ---- Plot LH
p_projection_lh <- ggplot(df_summary_lh, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(title = "Projection for pooled data: mOFC-lh",
       x = "Pain Phenotype",
       y = "Mean ± SE") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15),
    legend.position = "none"
  )

# ---- Plot RH
p_projection_rh <- ggplot(df_summary_rh, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(#title = "Projection of CP EEG data: mOFC-rh",
    title = "Projection for pooled data: mOFC-rh",
    #title = "CP EEG data projection: mOFC-rh",
    x = "Pain Phenotype",
    y = "Mean ± SE") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15),
    legend.position = "none"
  )

# ---- Arrange side by side
grid.arrange(p_projection_lh, p_projection_rh, nrow = 1, widths = c(1, 1))

pooled_p_projection_lh <- p_projection_lh 
pooled_p_projection_rh <- p_projection_rh 



grid.arrange(cp_p_projection_lh, cp_p_projection_rh, 
             pooled_p_projection_lh, pooled_p_projection_rh, 
             nrow = 2, ncol= 2, heights = c(1,1), widths = c(1,1)) 
#nrow = 1, ncol= 2, heights = 1, widths = c(1,1)) 



###########
# Run one-way ANOVA for left hemisphere
anova_lh <- aov(ProjectionScore_lh ~ Group, data = df_projection)
summary(anova_lh)

# Post-hoc: Tukey HSD for left hemisphere
tukey_lh <- TukeyHSD(anova_lh)
print(tukey_lh)

# Run one-way ANOVA for right hemisphere
anova_rh <- aov(ProjectionScore_rh ~ Group, data = df_projection)
summary(anova_rh)

# Post-hoc: Tukey HSD for right hemisphere
tukey_rh <- TukeyHSD(anova_rh)
print(tukey_rh)


# Run one-way ANOVA  (for both hemespheres)
anova_results <- aov(ProjectionScore  ~ Group, data = df_projection)
summary(anova_results)

# Post-hoc: Tukey HSD   (for both hemespheres)
tukey_results <- TukeyHSD(anova_results)
print(tukey_results)


##################
# Run ANCOVA for left hemisphere
ancova_lh <- aov(ProjectionScore_lh ~ Group + Pain_Intensity + Age + Sex, data = df_projection)
summary(ancova_lh)

# Run ANCOVA for right hemisphere
ancova_rh <- aov(ProjectionScore_rh ~ Group + Pain_Intensity + Age + Sex, data = df_projection)
summary(ancova_rh)

# Run ANCOVA for right hemisphere
ancova_results <- aov(ProjectionScore ~ Group + Pain_Intensity + Age + Sex, data = df_projection)
summary(ancova_results)






#####
###
library(dplyr)
library(ggplot2)

# ---- Summarize data: compute mean and SE per group
df_summary <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore, na.rm = TRUE),
    SE = sd(ProjectionScore, na.rm = TRUE) / sqrt(n())
  )

# ---- Plot mean ± SE as point and error bar
p_projection <- ggplot(df_summary, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(
    title = "EEG Marker Projection Score",
    x = "Pain Phenotype",
    y = "Mean ± SE"
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12)
  )

print(p_projection)





## lh
# ---- Summarize data: compute mean and SE per group
df_summary <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore_lh, na.rm = TRUE),
    SE = sd(ProjectionScore_lh, na.rm = TRUE) / sqrt(n())
  )

# ---- Plot mean ± SE as point and error bar
p_projection_lh <- ggplot(df_summary, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(
    title = paste("EEG Marker Projection Score:", "mOFC-lh"), #"EEG Marker Projection Score",
    x = "Pain Phenotype",
    y = "Mean ± SE"
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12)
  )

print(p_projection_lh)




## rh
# ---- Summarize data: compute mean and SE per group
df_summary <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore_rh, na.rm = TRUE),
    SE = sd(ProjectionScore_lh, na.rm = TRUE) / sqrt(n())
  )

# ---- Plot mean ± SE as point and error bar
p_projection_rh <- ggplot(df_summary, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(
    title = paste("EEG Marker Projection Score:", "mOFC-rh"), #"EEG Marker Projection Score",
    x = "Pain Phenotype",
    y = "Mean ± SE"
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12)
  )

print(p_projection_rh)


grid.arrange(p_projection_lh,  
             p_projection_rh, 
             nrow = 1, ncol= 2, heights = 1, widths = c(1,1)) 

## end of the code