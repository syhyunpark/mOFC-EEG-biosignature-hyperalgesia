##  6/3/2025
## Chronic pancreatitis (CP) cohort analysis 

# -- set directories 
rm(list = ls())

#setwd("/Users/hyung/Dropbox/LBP/github_mOFC_signature") # set a working directory where the files are saved 
library(readxl) 
library(tableone) 
library(rhdf5) 
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr) 
library(gridExtra)
library(grid)
library(pROC)
library(caret)
 
file_path <- "25_05_08Phenotype_CP.xlsx"
df <- read_excel(file_path, sheet = 1) 
colnames(df) 

# Ensure variable types are correct
df$Sex <- as.factor(df$Sex)
df$Phenotype <- factor(df$Phenotype, 
                       levels = c("No Hyperalgesia", "Segmental Hyperalgesia", "Widespread Hyperalgesia"))

# Specify variable types
factorVars <- c("Sex") #, "central_analgesic", "peripheral_treatment", "behavioral_modulation")
contVars <- c("BPI_Baseline", "Age")
allVars <- c(factorVars, contVars)

# Create Table 1
table1 <- CreateTableOne(vars = allVars,
                         strata = "Phenotype",
                         data = df,
                         factorVars = factorVars)

# Print Table 1 with means and percentages
print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)

# Convert table to data frame
table1_df <- print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

# Write to CSV file
#write.csv(table1_df, file = "Table1_CP_PainPhenotype.csv", row.names = TRUE)


#dat2save <- readRDS("CP_EEG_06032025.RData")  
eeg_dat_Y <- dat2save$eeg_dat_Y
eeg_dat_J <- dat2save$eeg_dat_J

# -- all conditions 
conditions <- c("Abdomen_32", "Abdomen_512","Forearm_32", "Forearm_512")

# -- ROIs 
ROI_labels <- c("rACC-lh", "dACC-lh", "S1-lh",
                "Ins-lh", "dlPFC-lh", "mOFC-lh",
                "rACC-rh", "dACC-rh", "S1-rh",
                "Ins-rh", "dlPFC-rh", "mOFC-rh")

pain_levels <- c("Widespread Hyperalgesia", "Segmental Hyperalgesia", "No Hyperalgesia")  
Phenotype = df$Phenotype <- factor(df$Phenotype, levels = pain_levels)
Subject <- df$record_id 



tt <- 17:43  # corresponds to [-0.3, 0.7] sec
ff <- 1:45
TT <- length(tt)
FF <- length(ff) 
# time sequence from -0.3 to 0.7 sec 
tt_sec <- seq(-0.3, 0.7, length.out = TT)
desired_times <- c(-0.2, 0, 0.2, 0.4, 0.6)
breaks_labels <- as.character(desired_times)

# -- Loop through brain regions

###########
# mOFC_lh
region = 6
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    print(j) 
    if(eeg_dat_J[i,j]){   
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Phenotype
      subject_group <- df$Phenotype[df$record_id == Subject[i]]
      
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
dat2plot$condition <- factor(dat2plot$condition, conditions) 
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Abdomen_32" = "Abdomen 32", 
                                    "Abdomen_512"= "Abdomen 512", 
                                    "Forearm_32" = "Forearm 32",
                                    "Forearm_512" = "Forearm 512") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Abdomen 32", "Abdomen 512", "Forearm 32", "Forearm 512"))

df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = median(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

pp_mOFC_lh <- ggplot(df_summary, aes(x = time_sec, y = frequency, fill = avg_value)) +
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

print(pp_mOFC_lh) 

###########
# mOFC_rh
region = 12
dat2plot <- NULL

for(i in 1:length(Subject)){
  print(i) 
  for(j in 1:length(conditions)){
    print(j) 
    if(eeg_dat_J[i,j]){   
      grid <- expand.grid(tt = tt, ff = ff) 
      tmp <- eeg_dat_Y[[i]][[j]][tt, ff, region]  
      grid$mat <- c(tmp)  
      grid$condition <- rep(conditions[j], length(tt) * length(ff))
      
      # assign group based on Phenotype
      subject_group <- df$Phenotype[df$record_id == Subject[i]]
      
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
dat2plot$condition <- factor(dat2plot$condition, conditions) 
dat2plot$condition <- dplyr::recode(dat2plot$condition, 
                                    "Abdomen_32" = "Abdomen 32", 
                                    "Abdomen_512"= "Abdomen 512", 
                                    "Forearm_32" = "Forearm 32",
                                    "Forearm_512" = "Forearm 512") 
dat2plot$condition <- factor(dat2plot$condition, levels = c("Abdomen 32", "Abdomen 512", "Forearm 32", "Forearm 512"))

df_summary <- dat2plot %>%
  group_by(tt, ff, condition, group) %>%
  summarise(avg_value = median(mat, na.rm = TRUE), .groups = "drop")

# match tt to seconds
df_summary <- df_summary %>%
  mutate(
    time_sec = tt_sec[tt - min(tt) + 1],
    frequency = ff
  )

pp_mOFC_rh <- ggplot(df_summary, aes(x = time_sec, y = frequency, fill = avg_value)) +
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

print(pp_mOFC_rh) 

grid.arrange(pp_mOFC_lh,  
             pp_mOFC_rh,  
             nrow = 2, ncol= 1, heights = c(1,1), widths = 1)  





##########################
###### projection   ######
##########################

# Define time and frequency indices
tt <- 17:43  # full trimmed time index  
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
conditions <- c("Abdomen_32", "Abdomen_512","Forearm_32", "Forearm_512")
J <- length(conditions)

n <- length(Subject) 
res_proj = res_proj_lh = res_proj_rh <- matrix(NA, nrow = n, ncol = J) 


# Loop through subjects and conditions
for (i in 1:n) {
  for (j in 1:J) {
    if (eeg_dat_J[i,j]) {
      res_proj_lh[i,j] <- sum(eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-lh")] * custom_pattern)
      res_proj_rh[i,j] <- sum(eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-rh")] * custom_pattern)
      res_proj[i,j] <- sum(eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-lh")] * custom_pattern) + 
        sum(eeg_dat_Y[[i]][[j]][tt, ff, which(ROI_labels == "mOFC-rh")] * custom_pattern)
    }
  }
}
 


# -- Initialize output data
proj_plot_data <- NULL

for (i in 1:n) {
  subj_id <- Subject[i]
  subj_group <- df$Phenotype[df$record_id == subj_id]
  
  for (j in 1:J) {
    if (eeg_dat_J[i,j]) {
      score <- res_proj[i,j]
      score_lh <- res_proj_lh[i,j]
      score_rh <- res_proj_rh[i,j]
      proj_plot_data <- rbind(proj_plot_data, data.frame(
        Subject = subj_id,
        Group = subj_group,
        Condition = conditions[j],
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
                               levels = c("Widespread Hyperalgesia", "Segmental Hyperalgesia", "No Hyperalgesia"))
proj_plot_data$Condition <- factor(proj_plot_data$Condition, levels = conditions)

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
contrast_vec <- c(0, 1, 1, 1)  # custom contrast
names(contrast_vec) <- c("Abdomen_32", "Abdomen_512", "Forearm_32", "Forearm_512")

# Ensure conditions match
condition_levels <- conditions  # should be c("Abdomen_32", ..., etc.)
contrast_vec <- contrast_vec[condition_levels] 

# ---- Compute contrast score per subject
projection_scores = projection_scores_lh = projection_scores_rh <- rep(NA, n)

for (i in 1:n) {
  projection_scores[i]    <- sum(res_proj[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_lh[i] <- sum(res_proj_lh[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_rh[i] <- sum(res_proj_rh[i, ] * contrast_vec, na.rm = TRUE) 
}

# ---- Prepare plotting data
df_sub <- df[df$record_id %in% Subject, ]  # filter to match Subject
df_sub <- df_sub[match(Subject, df_sub$record_id), ]  # ensure order matches

df_projection <- data.frame(
  Subject = Subject,
  Group = df_sub$Phenotype,
  BPI_Baseline = df_sub$BPI_Baseline, 
  Age = df_sub$Age, 
  Sex = df_sub$Sex, 
  ProjectionScore_lh = projection_scores_lh,
  ProjectionScore_rh = projection_scores_rh,
  ProjectionScore = projection_scores
)

df_projection <- df_projection %>%
  filter(!is.na(ProjectionScore)) %>%
  mutate(Group = factor(Group,
                        levels = c("Widespread Hyperalgesia", "Segmental Hyperalgesia", "No Hyperalgesia")))



#############
p_projection_lh <- ggplot(df_projection, 
                          aes(x = Group, 
                              y = ProjectionScore_lh, 
                              fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste("Projection onto CP EEG data:", "mOFC-lh"),
       x = "Pain Phenotype",
       y = "Projection Score") +
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

# Show plot
print(p_projection_lh)

p_projection_rh <- ggplot(df_projection, 
                          aes(x = Group, 
                              y = ProjectionScore_rh, 
                              fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste("Projection onto CP EEG data:", "mOFC-rh"),
       x = "Pain Phenotype",
       y = "Projection Score") +
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

# Show plot
print(p_projection_rh)


grid.arrange(p_projection_lh,  
             p_projection_rh, 
             nrow = 1, ncol= 2, heights = 1, widths = c(1,1)) 




# ---- Compute contrast score per subject
projection_scores = projection_scores_lh = projection_scores_rh <- rep(NA, n)

for (i in 1:n) {
  projection_scores[i]    <- sum(res_proj[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_lh[i] <- sum(res_proj_lh[i, ] * contrast_vec, na.rm = TRUE)
  projection_scores_rh[i] <- sum(res_proj_rh[i, ] * contrast_vec, na.rm = TRUE)
}

# ---- Prepare plotting data
df_sub <- df[df$record_id %in% Subject, ]
df_sub <- df_sub[match(Subject, df_sub$record_id), ]

df_projection <- data.frame(
  Subject = Subject,
  Group = df_sub$Phenotype,
  BPI_Baseline = df_sub$BPI_Baseline,
  Age = df_sub$Age,
  Sex = df_sub$Sex,
  ProjectionScore_lh = projection_scores_lh,
  ProjectionScore_rh = projection_scores_rh,
  ProjectionScore = projection_scores
)

df_projection <- df_projection %>%
  filter(!is.na(ProjectionScore)) %>%
  mutate(Group = factor(Group,
                        levels = c("Widespread Hyperalgesia", "Segmental Hyperalgesia", "No Hyperalgesia")))

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
  labs(title = "Projection for CP data: mOFC-lh",
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
  labs(title = "Projection for CP data: mOFC-rh", 
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

 



##########


# ---- Summarize: compute mean and SE   

df_summary <- df_projection %>%
  group_by(Group) %>%
  summarize(
    Mean = mean(ProjectionScore, na.rm = TRUE),
    SE = sd(ProjectionScore, na.rm = TRUE) / sqrt(n())
  )

# ---- Plot LH
p_projection  <- ggplot(df_summary, aes(x = Group, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme_minimal() +
  labs(#title = "EEG Marker Projection Score",
    title = "Projection for CP data",
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

p_projection





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


# Run one-way ANOVA for both
anova_both <- aov(ProjectionScore ~ Group, data = df_projection)
summary(anova_both)

# Post-hoc: Tukey HSD for both
tukey_both <- TukeyHSD(anova_both)
print(tukey_both)


df_projection$Widespread_vs_Other <- ifelse(df_projection$Group == "Widespread Hyperalgesia", 1, 0)
t.test_lh <- t.test(ProjectionScore_lh ~ Widespread_vs_Other, data = df_projection)
t.test_lh
wilcox_test_lh <- wilcox.test(ProjectionScore_lh ~ Widespread_vs_Other, data = df_projection)
wilcox_test_lh

t.test_rh <- t.test(ProjectionScore_rh ~ Widespread_vs_Other, data = df_projection)
t.test_rh
wilcox_test_rh <- wilcox.test(ProjectionScore_rh ~ Widespread_vs_Other, data = df_projection)
wilcox_test_rh

wilcox_test_both <- wilcox.test(ProjectionScore  ~ Widespread_vs_Other, data = df_projection)
wilcox_test_both


# Run ANCOVA for left hemisphere
ancova_lh <- aov(ProjectionScore_lh ~ Group + BPI_Baseline + Age + Sex, data = df_projection)
summary(ancova_lh)

# Run ANCOVA for right hemisphere
ancova_rh <- aov(ProjectionScore_rh ~ Group + BPI_Baseline + Age + Sex, data = df_projection)
summary(ancova_rh)

# Run ANCOVA for both 
ancova <- aov(ProjectionScore ~ Group + BPI_Baseline + Age + Sex, data = df_projection)
summary(ancova)



####
# Step 2: Prepare the subset to merge — must match Subject IDs
df_sub <- df[df$record_id %in% df_projection$Subject, ]
df_sub <- df_sub[match(df_projection$Subject, df_sub$record_id), ]  # ensure correct order

# Step 3: Merge with df_contrast
df_projection_merged <- cbind(df_projection,  
                              BPI_3mo = df_sub$BPI_3mo,
                              Diff_BPI = df_sub$Diff_BPI)
colnames(df_projection_merged)

df_projection_merged$Group
df_projection_merged$Any_Hyperalgesia <- ifelse(
  df_projection_merged$Group == "No Hyperalgesia", 0, 1)


## predictive models   
model_base <- lm(Diff_BPI ~ Widespread_vs_Other, data = df_projection_merged)
summary(model_base) 

model_base2 <- lm(Diff_BPI ~ Group, data = df_projection_merged)
summary(model_base2)
anova(model_base, model_base2)

model_extended2 <- lm(Diff_BPI ~ Widespread_vs_Other * ProjectionScore, data = df_projection_merged)
summary(model_extended2)
 
summary(model_base)$adj.r.squared  # Widespread_vs_Other
summary(model_base2)$adj.r.squared # group
summary(model_extended2)$adj.r.squared # interaction  

summary(model_base)$r.squared   # Widespread_vs_Other
summary(model_base2)$r.squared  # group
summary(model_extended2)$r.squared # interaction  
 


library(ggplot2)
 
# Convert indicator to factor for plotting
df_projection_merged$Widespread_vs_Other_F <- factor(
  df_projection_merged$Widespread_vs_Other,
  levels = c(0, 1),
  labels = c("Segm/No Hyperalgesia", "Widespread Hyperalgesia")
) 
# Plot with regression lines stratified by Hyperalgesia group
ggplot(df_projection_merged, aes(x = ProjectionScore, y = Diff_BPI, color = Widespread_vs_Other_F)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  labs(
    title = "Interaction Between EEG Marker and Hyperalgesia on Pain Change", 
    x = "EEG Nociplasticity Marker",
    y = "Pain Change", #"Δ BPI (3mo - Baseline)",
    color = "Hyperalgesia Group"
  ) +
  scale_color_manual(values = c("blue", "darkorange")) +
  #scale_x_continuous(limits = c(-200, 300)) +
  scale_x_continuous(limits = c(-500, 1500)) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")



 
###
library(pROC)
library(caret)
library(dplyr)

# Prepare the data
df <- df_projection_merged %>%
  filter(complete.cases(Diff_BPI, Any_Hyperalgesia, Widespread_vs_Other, ProjectionScore)) %>%
  mutate(
    pain_reduction_percent = if_else(BPI_Baseline > 0,
                                     (Diff_BPI / BPI_Baseline) * 100,
                                     0),
    Responder = as.numeric(pain_reduction_percent <= -30)
  )

evaluate_model_loocv <- function(df, formula) {
  n <- nrow(df)
  fitted_scores <- numeric(n)
  
  for (i in 1:n) {
    train_data <- df[-i, ]
    test_data <- df[i, , drop = FALSE]
    
    X_train <- model.matrix(formula, data = train_data)
    y_train <- train_data$Diff_BPI
    X_test <- model.matrix(formula, data = test_data)
    
    missing_cols <- setdiff(colnames(X_train), colnames(X_test))
    for (col in missing_cols) {
      X_test <- cbind(X_test, setNames(rep(0, nrow(X_test)), col))
    }
    X_test <- X_test[, colnames(X_train), drop = FALSE]
    
    fit <- lm.fit(X_train, y_train)
    fitted_scores[i] <- X_test %*% fit$coefficients
  }
  
  df_out <- df %>%
    mutate(
      fitted_scores = fitted_scores,
      predicted_pain_reduction_percent = if_else(BPI_Baseline > 0,
                                                 (fitted_scores / BPI_Baseline) * 100,
                                                 0),
      predicted_Responder = as.numeric(predicted_pain_reduction_percent <= -30)
    )
  
  # ROC and AUC
  roc_obj <- roc(response = df_out$Responder, predictor = -df_out$fitted_scores)
  auc_val <- auc(roc_obj)
  
  # Confusion matrix
  conf_mat <- confusionMatrix(
    factor(df_out$predicted_Responder, levels = c(0, 1)),
    factor(df_out$Responder, levels = c(0, 1)),
    positive = "1"
  )
  
  # PPV and NPV
  tp <- conf_mat$table[2, 2]
  fp <- conf_mat$table[2, 1]
  fn <- conf_mat$table[1, 2]
  tn <- conf_mat$table[1, 1]
  
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  
  return(list(
    AUC = round(auc_val, 3),
    Accuracy = round(conf_mat$overall["Accuracy"], 3),
    Sensitivity = round(conf_mat$byClass["Sensitivity"], 3),
    Specificity = round(conf_mat$byClass["Specificity"], 3),
    PPV = round(ppv, 3),
    NPV = round(npv, 3),
    ROC = roc_obj #,
    #Predictions = df_out$predicted_Responder,
    #Truth = df_out$Responder
  ))
}

# Run all models
results_base <- evaluate_model_loocv(df, Diff_BPI ~ Widespread_vs_Other)
results_group <- evaluate_model_loocv(df, Diff_BPI ~ Group)
results_interact <- evaluate_model_loocv(df, Diff_BPI ~ Widespread_vs_Other * ProjectionScore)

# Print summary
print(results_base)
print(results_group)
print(results_interact)

## end of the code 
