library(readxl)
library(janitor)
library(dplyr)
library(MASS)
library(broom.mixed)
library(emmeans)
library(robustbase)

#load data
d = read_excel("dataset_github.xlsx", sheet = "20-40-60") %>%
  clean_names() %>%
  mutate(
    mu_id = as.factor(mu_id),
    intensity = as.factor(intensity),
    contraction = as.factor(contraction),
    group = as.factor(group),
    participant = as.factor(participant),
    normalized_brace_height = normalized_brace_height*100,
    mu_id_intensity = paste(mu_id, intensity),
    )


#calculate mean
df_mean = d %>% group_by(participant,group,sex) %>%
  summarise(age = mean(age, na.rm = TRUE),
            recruitment_threshold = mean(recruitment_threshold, na.rm = TRUE),
            deltaf = mean(deltaf, na.rm = TRUE),
            maximal_discharge_rate_of_the_svr = mean(maximal_discharge_rate_of_the_svr, na.rm = TRUE),
            normalized_brace_height = mean(normalized_brace_height, na.rm = TRUE),
            handgrip = mean(handgrip, na.rm = TRUE),
            tug = mean(tug, na.rm = TRUE),
            five_sts_power_bm = mean(x5sts_power_bm, na.rm = TRUE),
            five_sts = mean(x5sts, na.rm = TRUE),
            tug = mean(tug, na.rm = TRUE),
            four_m_usual = mean(x4m_usual, na.rm = TRUE),
            four_m_fast = mean(x4m_fast, na.rm = TRUE),
            four_square = mean(x4_square, na.rm = TRUE),
            body_mass = mean(body_mass, na.rm = TRUE),
            force = mean(force, na.rm = TRUE),
            torque = ((((force*226.8)/(20*5*5000))*1000)*9.81*0.1405),
            torque_bm = (torque/body_mass),
            height = mean(height, na.rm = TRUE),
            bmi = mean(bmi, na.rm = TRUE),
            fat_perc = mean(fm_percent_weight_number_1, na.rm = TRUE),
            llfdi_fdc = mean(llfdi_dc_total_frequency_dimension_score_scaled, na.rm = TRUE),
            llfdi_ldc = mean(llfdi_dc_total_limitation_dimension_scaled, na.rm = TRUE),
            llfdi_fc = mean(llfdi_fc_overall_function_scaled, na.rm = TRUE),
            katz = (mean(katz_1, na.rm = TRUE) + mean(katz_2, na.rm = TRUE) + mean(katz_3, na.rm = TRUE) +
                      mean(katz_4, na.rm = TRUE) + mean(katz_5, na.rm = TRUE) + mean(katz_6, na.rm = TRUE)),
            smm = (mean(right_leg_smm_kg_number_1, na.rm = TRUE) + 
                     mean(left_leg_smm_kg_number_1, na.rm = TRUE) +
                     mean(right_arm_smm_kg_number_1, na.rm = TRUE) +
                     mean(left_arm_smm_kg_number_1, na.rm = TRUE)),
            smm_index = (smm/height^2),
            total_smm = (mean(smm_kg_number_1, na.rm = TRUE)),
            ipaq_e_q1 = (mean(ipaq_e_q1, na.rm = TRUE)*60),
            ipaq_e_q2 = ((mean(ipaq_e_q2a, na.rm = TRUE)*(mean(ipaq_e_q2b, na.rm = TRUE))/7)),
            ipaq_e_q3 = ((mean(ipaq_e_q3a, na.rm = TRUE)*(mean(ipaq_e_q3b, na.rm = TRUE))/7)),
            ipaq_e_q4 = ((mean(ipaq_e_q4a, na.rm = TRUE)*(mean(ipaq_e_q4b, na.rm = TRUE))/7))
  )
d_mean


# Define BMI categories
df_mean <- df_mean %>%
  mutate(
    bmi_category = case_when(
      bmi < 18.5 ~ "Underweight",
      bmi >= 18.5 & bmi < 25 ~ "Normal weight",
      bmi >= 25 & bmi < 30 ~ "Overweight",
      bmi >= 30 ~ "Obese"
    )
  )

# Count participants per group and sex
bmi_counts <- df_mean %>%
  group_by(group, sex, bmi_category) %>%
  summarise(count = n(), .groups = "drop")

# Print results
print(bmi_counts)

# Fit the robust model - change the dependent as needed.
fit <- rlm(five_sts ~ group + sex, data = d_mean)


# Create a tidy summary of fixed effects
tidy_output <- tidy(fit, effects = "fixed", conf.int = TRUE)

# Extract only main effects and interactions
main_effects_interactions <- tidy_output 
# View the main effects and interactions
print(main_effects_interactions)

# Compute confidence intervals
confint(fit)
summary(fit)

# Compute estimated marginal means
emm <- emmeans(fit, ~ group)
summary(emm, infer = TRUE)

# Pairwise comparisons with confidence intervals
pairs(emm, infer = TRUE)
confint(emm)

# Descriptive statistics for group * sex
descriptive_stats <- d_mean %>%
  group_by(group) %>%
  summarize(
    Median = median(katz),
    Q25 = quantile(katz, 0.25),
    Q75 = quantile(katz, 0.75),
    .groups = "drop"
  )

# View results
print(descriptive_stats)

#load mu number data
d_mu = read_excel("Motor unit numbers study 1.xlsx", sheet = "stats2") %>%
  clean_names() %>%
  mutate(
    group = as.factor(group),
    sex = as.factor(sex),
    intensity = as.factor(intensity),
    participant = as.factor(participant)
  )

# Descriptive statistics for group * sex
descriptive_stats <- d_mu %>%
  group_by(group, sex, intensity) %>%
  summarize(
    sum = sum(mus),
    Median = median(mus),
    Q25 = quantile(mus, 0.25),
    Q75 = quantile(mus, 0.75),
    .groups = "drop"
  )

descriptive_stats <- d_mu %>%
  summarize(
    sum_mus = sum(mus),
    sum_trains = sum(trains),
    Median = median(mus),
    Q25 = quantile(mus, 0.25),
    Q75 = quantile(mus, 0.75),
    .groups = "drop"
  )

# View results
print(descriptive_stats)

table(d_mu$group, d_mu$sex, d_mu$intensity)

# Fit the robust model
fit <- rlm(mus ~ group*sex*intensity, data = d_mu)


# Create a tidy summary of fixed effects
tidy_output <- tidy(fit, effects = "fixed", conf.int = TRUE)

# Extract only main effects and interactions
main_effects_interactions <- tidy_output 
# View the main effects and interactions
print(main_effects_interactions)

# Compute confidence intervals
confint(fit)
summary(fit)

# Compute estimated marginal means
emm <- emmeans(fit, ~ group|intensity)
summary(emm, infer = TRUE)

# Pairwise comparisons with confidence intervals
pairs(emm, infer = TRUE)
confint(emm)

#Two-part regression for IPAQ-E
# Logistic regression to model the probability of reporting activity
fit_logit <- glmrob(
  I(ipaq_e_q1 > 0) ~ group,
  family = binomial(link = "logit"),
  data = df_mean,
  method = "Mqle"
)

# Summary of the logistic model
summary(fit_logit)

# Extract coefficients
coef_logit <- coef(fit_logit)

# Compute odds ratios
exp_coef <- exp(coef_logit)
exp_coef

coef_logit <- coef(fit_logit)
se_logit <- sqrt(diag(vcov(fit_logit)))  # Standard errors

z_value <- qnorm(0.975)  # 95% CI
lower_log <- coef_logit - z_value * se_logit
upper_log <- coef_logit + z_value * se_logit

lower_or <- exp(lower_log)
upper_or <- exp(upper_log)

# Combine results into a table
odds_ratios <- data.frame(
  Coefficient = coef_logit,
  Odds_Ratio = exp_coef,
  Lower_95_CI = lower_or,
  Upper_95_CI = upper_or
)

print(odds_ratios)


#Step 2: Gamma Regression for Non-Zero Values

# Filter non-zero values
df_nonzero <- subset(df_mean, ipaq_e_q4 > 0)

# Fit a gamma regression model for non-zero values
fit_gamma_nonzero <- glmrob(
  ipaq_e_q1 ~ group,
  family = Gamma(link = "log"),
  data = df_nonzero,
  method = "Mqle"
)

# Summary of the gamma model
summary(fit_gamma_nonzero)

coef_gamma <- coef(fit_gamma_nonzero)  # Coefficients (log-scale)
se_gamma <- sqrt(diag(vcov(fit_gamma_nonzero)))  # Standard errors

z_value <- qnorm(0.975)  # For 95% CI

lower_log <- coef_gamma - z_value * se_gamma  # Lower bound (log scale)
upper_log <- coef_gamma + z_value * se_gamma  # Upper bound (log scale)

lower_or <- exp(lower_log)
upper_or <- exp(upper_log)

# Combine results into a table
or_results <- data.frame(
  Coefficient = coef_gamma,
  Odds_Ratio = exp(coef_gamma),
  Lower_95_CI = lower_or,
  Upper_95_CI = upper_or
)

print(or_results)
summary(fit_gamma_nonzero)

