library(readxl)
library(dplyr)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(janitor)
library(cowplot)
library(emmeans)
library(lme4)
library(tidyr)
library(ggplot2)
library(broom.mixed)
library(robustlmm)
library(ggdist)
library(ggdist)
library(rlang)
library(patchwork)
library(ggbeeswarm) 

#load data
d = read_excel("dataset_github.xlsx", sheet = "20-40-60") %>%
  clean_names() %>%
  mutate(
    mu_id = as.factor(mu_id),
    intensity = as.factor(intensity),
    contraction = as.factor(contraction),
    group = as.factor(group),
    participant = as.factor(participant),
    mu_id_intensity = paste(mu_id, intensity),
    )

#calculate mean
d_mean = d %>% group_by(participant,intensity,group,sex) %>%
  summarise(age = mean(age, na.rm = TRUE),
            recruitment_threshold = mean(recruitment_threshold, na.rm = TRUE),
            deltaf = mean(deltaf, na.rm = TRUE),
            handgrip = mean(handgrip, na.rm = TRUE),
            tug = mean(tug, na.rm = TRUE),
            five_sts_power_bm = mean(x5sts_power_bm, na.rm = TRUE),
            handgrip = mean(handgrip, na.rm = TRUE),
            tug = mean(tug, na.rm = TRUE),
            four_m_usual = mean(x4m_usual, na.rm = TRUE),
            four_m_fast = mean(x4m_fast, na.rm = TRUE),
            four_square = mean(x4_square, na.rm = TRUE),
            body_mass = mean(body_mass, na.rm = TRUE),
            force = mean(force, na.rm = TRUE),
            torque = ((((force*226.8)/(20*5*5000))*1000)*9.81*0.1405),
            torque_bm = (torque/body_mass)
            )
d_mean


#################
###############
#############
###########
##########
#########
########
#######
######
#####
####
###
##
#Prediction models
# Remove rows with missing delta_f values
d_clean <- d[!is.na(d$deltaf), ]

#5x sts power vs deltaf

fit_rlmer_df_5sts <- rlmer(
  formula = deltaf ~ x5sts_power_bm * sex*intensity + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,  # Increase tolerance
  max.iter = 1000  # Increase maximum iterations
)
summary(fit_rlmer_df_5sts)

vcov_matrix <- vcov(fit_rlmer_df_5sts)

# Extract coefficients and variance-covariance matrix
fixed_effects <- fixef(fit_rlmer_df_5sts)
vcov_matrix <- vcov(fit_rlmer_df_5sts)

# Define formulas for each group (for x5sts_power_bm_scaled)
slopes <- list(
  beta_f_20 = "x5sts_power_bm",
  beta_f_40 = "x5sts_power_bm + x5sts_power_bm:intensityi40%",
  beta_f_60 = "x5sts_power_bm + x5sts_power_bm:intensityi60%",
  beta_m_20 = "x5sts_power_bm + x5sts_power_bm:sexmale",
  beta_m_40 = "x5sts_power_bm + x5sts_power_bm:sexmale + x5sts_power_bm:intensityi40% + x5sts_power_bm:sexmale:intensityi40%",
  beta_m_60 = "x5sts_power_bm + x5sts_power_bm:sexmale + x5sts_power_bm:intensityi60% + x5sts_power_bm:sexmale:intensityi60%"
)

# Function to calculate slopes, SEs, and confidence intervals
compute_slope <- function(formula, fixed_effects, vcov_matrix) {
  # Parse the formula
  terms <- unlist(strsplit(formula, " + ", fixed = TRUE))
  terms <- gsub(" ", "", terms)  # Clean spaces
  
  # Compute slope
  slope <- sum(fixed_effects[terms])
  
  # Compute SE
  var_terms <- vcov_matrix[terms, terms, drop = FALSE]
  se <- sqrt(sum(var_terms))
  
  # Confidence intervals
  ci_lower <- slope - 1.96 * se
  ci_upper <- slope + 1.96 * se
  
  return(c(beta = slope, se = se, ci_lower = ci_lower, ci_upper = ci_upper))
}

# Compute slopes and CIs for all groups
results <- lapply(slopes, compute_slope, fixed_effects = fixed_effects, vcov_matrix = vcov_matrix)

# Convert to a data frame for easier viewing
results_df_5sts <- do.call(rbind, results)
rownames(results_df_5sts) <- names(slopes)

# Print results
print(results_df_5sts)



# Convert results_df to a data frame and add rownames as a column
results_df_5sts <- as.data.frame(results_df_5sts)
results_df_5sts$group <- rownames(results_df_5sts)

# Reorder the 'group' column as a factor with the desired order
results_df_5sts$group <- factor(
  results_df_5sts$group,
  #levels = c("beta_f_20", "beta_f_40", "beta_f_60", "beta_m_20", "beta_m_40", "beta_m_60")
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_5sts$sex <- ifelse(
  grepl("beta_f", results_df_5sts$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_5sts, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +                                         # Plot points with confidence intervals
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Add a red dashed line at y = 0
  coord_flip() +                                              # Flip the axes for a horizontal layout
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),  # Frame around intensity labels
    strip.text = element_text(face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  ) +
  labs(
    title = "ΔF by 5STS",                                            # Set the title
    x = NULL,                        # Set x-axis label
    y = "Slope (β)"                                        # Set y-axis label
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")  # Define colors for sex
  ) +
  scale_x_discrete(                                           # Customize x-axis labels
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_5sts

plot_df_5sts



# For torque
fit_rlmer_df_torque <- rlmer(
  formula = deltaf ~ torque * sex*intensity + body_mass + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,
  max.iter = 1000
)
summary(fit_rlmer_df_torque)

# Extract coefficients and variance-covariance matrix
fixed_effects <- fixef(fit_rlmer_df_torque)
vcov_matrix <- vcov(fit_rlmer_df_torque)

# Define formulas for each group (for torque_bm_scaled)
slopes <- list(
  beta_f_20 = "torque",
  beta_f_40 = "torque + torque:intensityi40%",
  beta_f_60 = "torque + torque:intensityi60%",
  beta_m_20 = "torque + torque:sexmale",
  beta_m_40 = "torque + torque:sexmale + torque:intensityi40% + torque:sexmale:intensityi40%",
  beta_m_60 = "torque + torque:sexmale + torque:intensityi60% + torque:sexmale:intensityi60%"
)


# Compute slopes and CIs for all groups
results <- lapply(slopes, compute_slope, fixed_effects = fixed_effects, vcov_matrix = vcov_matrix)

# Convert to a data frame for easier viewing
results_df_torque <- do.call(rbind, results)
rownames(results_df_torque) <- names(slopes)

# Print results
print(results_df_torque)

# Convert results_df to a data frame and add rownames as a column
results_df_torque <- as.data.frame(results_df_torque)
results_df_torque$group <- rownames(results_df_torque)

# Reorder the 'group' column as a factor with the desired order
results_df_torque$group <- factor(
  results_df_torque$group,
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_torque$sex <- ifelse(
  grepl("beta_f", results_df_torque$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_torque, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  ) +
  labs(
    title = "ΔF by Peak Torque", 
    x = NULL,
    y = "Slope (β)"
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")  # Define colors for sex
  ) +
  scale_x_discrete(
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_torque

plot_df_torque


# For tug_scaled
fit_rlmer_df_tug <- rlmer(
  formula = deltaf ~ tug * sex*intensity + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,
  max.iter = 1000
)
summary(fit_rlmer_df_tug)

fixed_effects <- fixef(fit_rlmer_df_tug)
vcov_matrix <- vcov(fit_rlmer_df_tug)

slopes <- list(
  beta_f_20 = "tug",
  beta_f_40 = "tug + tug:intensityi40%",
  beta_f_60 = "tug + tug:intensityi60%",
  beta_m_20 = "tug + tug:sexmale",
  beta_m_40 = "tug + tug:sexmale + tug:intensityi40% + tug:sexmale:intensityi40%",
  beta_m_60 = "tug + tug:sexmale + tug:intensityi60% + tug:sexmale:intensityi60%"
)


results <- lapply(slopes, compute_slope, fixed_effects = fixed_effects, vcov_matrix = vcov_matrix)
results_df_tug <- do.call(rbind, results)
rownames(results_df_tug) <- names(slopes)

# Print results
print(results_df_tug)

results_df_tug <- as.data.frame(results_df_tug)
results_df_tug$group <- rownames(results_df_tug)
results_df_tug$group <- factor(
  results_df_tug$group,
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_tug$sex <- ifelse(
  grepl("beta_f", results_df_tug$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_tug, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  ) +
  labs(
    title = "ΔF by TUG", 
    x = NULL,
    y = "Slope (β)"
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")  # Define colors for sex
  ) +
  scale_x_discrete(
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_tug

plot_df_tug


# For x4m_fast_scaled
fit_rlmer_df_x4m_fast <- rlmer(
  formula = deltaf ~ x4m_fast * sex*intensity + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,
  max.iter = 1000
)
summary(fit_rlmer_df_x4m_fast)

# Extract coefficients and variance-covariance matrix
fixed_effects <- fixef(fit_rlmer_df_x4m_fast)
vcov_matrix <- vcov(fit_rlmer_df_x4m_fast)

# Define formulas for each group (for x4m_fast_scaled)
slopes <- list(
  beta_f_20 = "x4m_fast",
  beta_f_40 = "x4m_fast + x4m_fast:intensityi40%",
  beta_f_60 = "x4m_fast + x4m_fast:intensityi60%",
  beta_m_20 = "x4m_fast + x4m_fast:sexmale",
  beta_m_40 = "x4m_fast + x4m_fast:sexmale + x4m_fast:intensityi40% + x4m_fast:sexmale:intensityi40%",
  beta_m_60 = "x4m_fast + x4m_fast:sexmale + x4m_fast:intensityi60% + x4m_fast:sexmale:intensityi60%"
)

# Compute slopes and CIs for all groups
results <- lapply(slopes, compute_slope, fixed_effects = fixed_effects, vcov_matrix = vcov_matrix)

# Convert to a data frame for easier viewing
results_df_x4m_fast <- do.call(rbind, results)
rownames(results_df_x4m_fast) <- names(slopes)

# Print results
print(results_df_x4m_fast)

# Convert results_df to a data frame and add rownames as a column
results_df_x4m_fast <- as.data.frame(results_df_x4m_fast)
results_df_x4m_fast$group <- rownames(results_df_x4m_fast)

# Reorder the 'group' column as a factor with the desired order
results_df_x4m_fast$group <- factor(
  results_df_x4m_fast$group,
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_x4m_fast$sex <- ifelse(
  grepl("beta_f", results_df_x4m_fast$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_x4m_fast, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  ) +
  labs(
    title = "ΔF by FGS", 
    x = NULL,
    y = "Slope (β)"
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")  # Define colors for sex
  ) +
  scale_x_discrete(
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_x4m_fast

plot_df_x4m_fast

# For x4m_usual_scaled
fit_rlmer_df_x4m_usual <- rlmer(
  formula = deltaf ~ x4m_usual * sex * intensity + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,
  max.iter = 1000
)
summary(fit_rlmer_df_x4m_usual)

# Extract coefficients and variance-covariance matrix
fixed_effects_usual <- fixef(fit_rlmer_df_x4m_usual)
vcov_matrix_usual <- vcov(fit_rlmer_df_x4m_usual)

# Define formulas for each group (for x4m_usual_scaled)
slopes_usual <- list(
  beta_f_20 = "x4m_usual",
  beta_f_40 = "x4m_usual + x4m_usual:intensityi40%",
  beta_f_60 = "x4m_usual + x4m_usual:intensityi60%",
  beta_m_20 = "x4m_usual + x4m_usual:sexmale",
  beta_m_40 = "x4m_usual + x4m_usual:sexmale + x4m_usual:intensityi40% + x4m_usual:sexmale:intensityi40%",
  beta_m_60 = "x4m_usual + x4m_usual:sexmale + x4m_usual:intensityi60% + x4m_usual:sexmale:intensityi60%"
)

# Compute slopes and CIs for all groups
results_usual <- lapply(slopes_usual, compute_slope, fixed_effects = fixed_effects_usual, vcov_matrix = vcov_matrix_usual)

# Convert to a data frame for easier viewing
results_df_x4m_usual <- do.call(rbind, results_usual)
rownames(results_df_x4m_usual) <- names(slopes_usual)

# Print results
print(results_df_x4m_usual)

# Convert results_df to a data frame and add rownames as a column
results_df_x4m_usual <- as.data.frame(results_df_x4m_usual)
results_df_x4m_usual$group <- rownames(results_df_x4m_usual)

# Reorder the 'group' column as a factor with the desired order
results_df_x4m_usual$group <- factor(
  results_df_x4m_usual$group,
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_x4m_usual$sex <- ifelse(
  grepl("beta_f", results_df_x4m_usual$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_x4m_usual, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "ΔF by UGS", 
    x = NULL,
    y = "Slope (β)"
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")
  ) +
  scale_x_discrete(
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_x4m_usual

plot_df_x4m_usual


# For x4_square_scaled
fit_rlmer_df_x4_square <- rlmer(
  formula = deltaf ~ x4_square * sex*intensity + (1 | participant/mu_id_intensity),
  data = d_clean,
  rel.tol = 1e-8,
  max.iter = 1000
)
summary(fit_rlmer_df_x4_square)

fixed_effects <- fixef(fit_rlmer_df_x4_square)
vcov_matrix <- vcov(fit_rlmer_df_x4_square)

slopes <- list(
  beta_f_20 = "x4_square",
  beta_f_40 = "x4_square + x4_square:intensityi40%",
  beta_f_60 = "x4_square + x4_square:intensityi60%",
  beta_m_20 = "x4_square + x4_square:sexmale",
  beta_m_40 = "x4_square + x4_square:sexmale + x4_square:intensityi40% + x4_square:sexmale:intensityi40%",
  beta_m_60 = "x4_square + x4_square:sexmale + x4_square:intensityi60% + x4_square:sexmale:intensityi60%"
)

results <- lapply(slopes, compute_slope, fixed_effects = fixed_effects, vcov_matrix = vcov_matrix)
results_df_x4_square <- do.call(rbind, results)
rownames(results_df_x4_square) <- names(slopes)

results_df_x4_square <- as.data.frame(results_df_x4_square)
results_df_x4_square$group <- rownames(results_df_x4_square)

# Print results
print(results_df_x4m_fast)

results_df_x4_square$group <- factor(
  results_df_x4_square$group,
  levels = c("beta_m_60", "beta_m_40", "beta_m_20", "beta_f_60", "beta_f_40", "beta_f_20")
)

# Add a new column to classify male and female groups
results_df_x4_square$sex <- ifelse(
  grepl("beta_f", results_df_x4_square$group), "Female", "Male"
)

# Updated ggplot with color based on sex
ggplot(results_df_x4_square, aes(x = group, y = beta, ymin = ci_lower, ymax = ci_upper, color = sex)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  ) +
  labs(
    title = "ΔF by 4SST", 
    x = NULL,
    y = "Slope (β)"
  ) +
  scale_color_manual(
    values = c("Female" = "deeppink", "Male" = "darkblue")  # Define colors for sex
  ) +
  scale_x_discrete(
    labels = c(
      "beta_f_20" = "Female, i20%",
      "beta_f_40" = "Female, i40%",
      "beta_f_60" = "Female, i60%",
      "beta_m_20" = "Male, i20%",
      "beta_m_40" = "Male, i40%",
      "beta_m_60" = "Male, i60%"
    )
  ) -> plot_df_x4_square

plot_df_x4_square


# Combine the plots into a single figure with labels A to E
combined_corr_plot <- plot_grid(
  plot_df_x4m_usual + labs(tag = "A"),  
  plot_df_x4m_fast + labs(tag = "B"),
  plot_df_x4_square + labs(tag = "C"), 
  plot_df_tug + labs(tag = "D"),       
  plot_df_5sts + labs(tag = "E"),      
  plot_df_torque + labs(tag = "F"),    
  ncol = 2,                            # Arrange the plots into 2 columns
  rel_heights = c(.5, .5),               # Equal height for rows
  rel_widths = c(2, 2)                 # Equal width for columns
)

# Save the combined plot as an image file (optional)
ggsave("combined_corr_plot.png", combined_corr_plot, width = 10, height = 12, dpi = 300)

# Display the combined plot
print(combined_corr_plot)

#######
######
#####
####
###
##
# Prediction model scatterplots

#Delta F vs UGS
# Step 1: Create min/max per sex × intensity group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(x4m_usual, na.rm = TRUE),
    x_max = max(x4m_usual, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create restricted prediction grid
newdata <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    x4m_usual = seq(.$x_min, .$x_max, length.out = 100)
  )) %>%
  ungroup()

newdata$intensity <- factor(newdata$intensity, levels = levels(d_clean$intensity))

# Step 3: Predict fixed effects only
X <- model.matrix(~ x4m_usual * sex * intensity, newdata)
beta <- fixef(fit_rlmer_df_x4m_usual)
vcov_mat <- vcov(fit_rlmer_df_x4m_usual)
pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata$fit <- as.vector(pred)
newdata$ci_lower <- newdata$fit - 1.96 * se
newdata$ci_upper <- newdata$fit + 1.96 * se

# Step 4: Plot
p_deltaf_ugs <- ggplot() +
  geom_point(data = d_clean, aes(x = x4m_usual, y = deltaf, color = sex), alpha = 0.1, size = 1) +
  geom_ribbon(data = newdata, aes(x = x4m_usual, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata, aes(x = x4m_usual, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by Usual Gait Speed, Sex, and Intensity",
    x = expression("UGS (m·s"^-1*")"),
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )

#Delta F vs FGS

# Step 1: Create min/max per sex × intensity group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(x4m_fast, na.rm = TRUE),
    x_max = max(x4m_fast, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create restricted prediction grid
newdata <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    x4m_fast = seq(.$x_min, .$x_max, length.out = 100)
  )) %>%
  ungroup()

newdata$intensity <- factor(newdata$intensity, levels = levels(d_clean$intensity))

# Step 3: Predict fixed effects only
X <- model.matrix(~ x4m_fast * sex * intensity, newdata)
beta <- fixef(fit_rlmer_df_x4m_fast)
vcov_mat <- vcov(fit_rlmer_df_x4m_fast)
pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata$fit <- as.vector(pred)
newdata$ci_lower <- newdata$fit - 1.96 * se
newdata$ci_upper <- newdata$fit + 1.96 * se

# Step 4: Plot
p_deltaf_fgs <- ggplot() +
  geom_point(data = d_clean, aes(x = x4m_fast, y = deltaf, color = sex), alpha = 0.1, size = 1) +
  geom_ribbon(data = newdata, aes(x = x4m_fast, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata, aes(x = x4m_fast, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by Fast Gait Speed, Sex, and Intensity",
    x = expression("FGS (m·s"^-1*")"),
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )

#Delta f vs 4SST
# Step 1: Create min/max per sex × intensity group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(x4_square, na.rm = TRUE),
    x_max = max(x4_square, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create restricted prediction grid
newdata <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    x4_square = seq(.$x_min, .$x_max, length.out = 100)
  )) %>%
  ungroup()

newdata$intensity <- factor(newdata$intensity, levels = levels(d_clean$intensity))

# Step 3: Predict fixed effects only
X <- model.matrix(~ x4_square * sex * intensity, newdata)
beta <- fixef(fit_rlmer_df_x4_square)
vcov_mat <- vcov(fit_rlmer_df_x4_square)
pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata$fit <- as.vector(pred)
newdata$ci_lower <- newdata$fit - 1.96 * se
newdata$ci_upper <- newdata$fit + 1.96 * se

# Step 4: Plot
p_deltaf_4sst <- ggplot() +
  geom_point(data = d_clean, aes(x = x4_square, y = deltaf, color = sex), alpha = 0.1, size = 1) +
  geom_ribbon(data = newdata, aes(x = x4_square, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata, aes(x = x4_square, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by 4SST, Sex, and Intensity",
    x = "4SST (s)",
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )

# Step 1: Create min/max TUG per sex × intensity group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(tug, na.rm = TRUE),
    x_max = max(tug, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create restricted prediction grid
newdata <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    tug = seq(.$x_min, .$x_max, length.out = 100)
  )) %>%
  ungroup()

newdata$intensity <- factor(newdata$intensity, levels = levels(d_clean$intensity))

# Step 3: Predict fixed effects only
X <- model.matrix(~ tug * sex * intensity, newdata)
beta <- fixef(fit_rlmer_df_tug)
vcov_mat <- vcov(fit_rlmer_df_tug)
pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata$fit <- as.vector(pred)
newdata$ci_lower <- newdata$fit - 1.96 * se
newdata$ci_upper <- newdata$fit + 1.96 * se

# Step 4: Plot
p_deltaf_tug <- ggplot() +
  geom_point(data = d_clean, aes(x = tug, y = deltaf, color = sex), alpha = 0.1, size = 1) +
  geom_ribbon(data = newdata, aes(x = tug, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata, aes(x = tug, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by TUG, Sex, and Intensity",
    x = "TUG (s)",
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )


#Delta F vs 5STS
# Step 1: Define observed min/max for each group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(x5sts_power_bm, na.rm = TRUE),
    x_max = max(x5sts_power_bm, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create prediction grid per sex × intensity group
newdata_5sts <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    x5sts_power_bm = seq(.$x_min, .$x_max, length.out = 100)
  )) %>%
  ungroup()

newdata_5sts$intensity <- factor(newdata_5sts$intensity, levels = levels(d_clean$intensity))

# Step 3: Generate predictions (fixed effects only)
X <- model.matrix(~ x5sts_power_bm * sex * intensity, newdata_5sts)
beta <- fixef(fit_rlmer_df_5sts)
vcov_mat <- vcov(fit_rlmer_df_5sts)

pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata_5sts$fit <- as.vector(pred)
newdata_5sts$ci_lower <- newdata_5sts$fit - 1.96 * se
newdata_5sts$ci_upper <- newdata_5sts$fit + 1.96 * se

# Step 4: Plot without annotation labels
p_deltaf_5stst <- ggplot() +
  geom_point(data = d_clean, aes(x = x5sts_power_bm, y = deltaf, color = sex), alpha = 0.1, size = 1.2) +
  geom_ribbon(data = newdata_5sts, aes(x = x5sts_power_bm, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata_5sts, aes(x = x5sts_power_bm, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by 5STS Power, Sex, and Intensity",
    x = expression("5STS (W·kg"^-1*")"),
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )

#Delta F vs Peak torque

# Step 1: Create min/max torque per sex × intensity group
ranges <- d_clean %>%
  group_by(sex, intensity) %>%
  summarise(
    x_min = min(torque, na.rm = TRUE),
    x_max = max(torque, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create restricted prediction grid
newdata <- ranges %>%
  rowwise() %>%
  do(data.frame(
    sex = .$sex,
    intensity = .$intensity,
    torque = seq(.$x_min, .$x_max, length.out = 100),
    body_mass = mean(d_clean$body_mass[d_clean$sex == .$sex & d_clean$intensity == .$intensity], na.rm = TRUE)
  )) %>%
  ungroup()

newdata$intensity <- factor(newdata$intensity, levels = levels(d_clean$intensity))

# Step 3: Predict fixed effects only
X <- model.matrix(~ torque * sex * intensity + body_mass, newdata)
beta <- fixef(fit_rlmer_df_torque)
vcov_mat <- vcov(fit_rlmer_df_torque)
pred <- X %*% beta
se <- sqrt(diag(X %*% vcov_mat %*% t(X)))

newdata$fit <- as.vector(pred)
newdata$ci_lower <- newdata$fit - 1.96 * se
newdata$ci_upper <- newdata$fit + 1.96 * se

# Step 4: Plot
p_deltaf_torque <- ggplot() +
  geom_point(data = d_clean, aes(x = torque, y = deltaf, color = sex), alpha = 0.1, size = 1) +
  geom_ribbon(data = newdata, aes(x = torque, ymin = ci_lower, ymax = ci_upper, fill = sex), alpha = 0.3) +
  geom_line(data = newdata, aes(x = torque, y = fit, color = sex), linewidth = 1.2) +
  facet_wrap(~ intensity) +
  scale_color_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  scale_fill_manual(values = c("female" = "deeppink", "male" = "darkblue")) +
  labs(
    #title = "ΔF by peak torque, Sex, and Intensity",
    x = expression("Peak torque (N·m"^-1*")"),
    y = "ΔF (pps)"
  ) +
  theme_bw(base_size = 14) +
  ylim(-1.0, 12.0) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "gray98", color = "black"),
    legend.position = "none"
  )


# Combine plots into 2 columns, 3 rows
combined_df_function_plot <- plot_grid(
  p_deltaf_ugs    + labs(tag = "A"),
  p_deltaf_fgs    + labs(tag = "B"),
  p_deltaf_4sst   + labs(tag = "C"),
  p_deltaf_tug    + labs(tag = "D"),
  p_deltaf_5stst  + labs(tag = "E"),
  p_deltaf_torque + labs(tag = "F"),
  ncol = 2,        # 2 columns
  align = "hv",    # Horizontal and vertical alignment
  axis = "tblr"    # Align axes on all sides
)

# Optional: Save figure
ggsave("combined_df_function_plot_2col.png", combined_df_function_plot, width = 14, height = 18, dpi = 300)

# Display
print(combined_df_function_plot)

