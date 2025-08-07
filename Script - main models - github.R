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
    normalized_brace_height = normalized_brace_height*100,
    mu_id_intensity = paste(mu_id, intensity),
    )

# Recruitment threshold binning
d$rt_bin <- cut(
  d$recruitment_threshold,
  breaks = c(-5, 20, 40, 65),
  right = TRUE,
  include.lowest = TRUE,
  labels = c("0–20", "20–40", "40–60")
)


#calculate mean
d_mean = d %>% group_by(participant,intensity,group,sex,rt_bin) %>%
  summarise(age = mean(age, na.rm = TRUE),
            recruitment_threshold = mean(recruitment_threshold, na.rm = TRUE),
            deltaf = mean(deltaf, na.rm = TRUE),
            maximal_discharge_rate_of_the_svr = mean(maximal_discharge_rate_of_the_svr, na.rm = TRUE),
            normalized_brace_height = mean(normalized_brace_height, na.rm = TRUE),
            handgrip = mean(handgrip, na.rm = TRUE),
            tug = mean(tug, na.rm = TRUE),
            five_sts_power_bm = mean(x5sts_power_bm, na.rm = TRUE),
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
            ipaq_e_q1 = (mean(ipaq_e_q1, na.rm = TRUE)*60),
            ipaq_e_q2 = ((mean(ipaq_e_q2a, na.rm = TRUE)*(mean(ipaq_e_q2b, na.rm = TRUE))/7)),
            ipaq_e_q3 = ((mean(ipaq_e_q3a, na.rm = TRUE)*(mean(ipaq_e_q3b, na.rm = TRUE))/7)),
            ipaq_e_q4 = ((mean(ipaq_e_q4a, na.rm = TRUE)*(mean(ipaq_e_q4b, na.rm = TRUE))/7))
  )
d_mean

#########################
#########################
####### DELTA F #########
#########################
#########################

# Initialize storage for deltaf
all_contrasts_deltaf <- list()
all_summaries_deltaf <- list()
all_emmeans_deltaf <- list()
all_intensity_contrasts_deltaf <- list()
all_models_deltaf <- list()

for (bin in levels(d$rt_bin)) {
  cat("Running model for deltaf in rt_bin:", bin, "\n")
  
  d_bin <- d %>% filter(rt_bin == bin)
  if (n_distinct(d_bin$group) < 2) {
    cat("  ➔ Skipping: less than 2 groups present\n")
    next
  }
  
  model_formula <- if (bin == "40–60") {
    deltaf ~ group + sex +  (1 | participant / mu_id_intensity)
  } else {
    deltaf ~ group * intensity + sex +  (1 | participant / mu_id_intensity)
  }

## Exploratory model adjusted by age    
#  model_formula <- if (bin == "40–60") {
#    deltaf ~ group + sex + scale(age) +  (1 | participant / mu_id_intensity)
#  } else {
#    deltaf ~ group * intensity + sex + scale(age) + (1 | participant / mu_id_intensity)
#  }
  
  emm_options(lmerTest.limit = 15000)
  fit <- rlmer(model_formula, data = d_bin)
  all_models_deltaf[[bin]] <- fit
  
  model_summary <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  model_summary$rt_bin <- bin
  all_summaries_deltaf[[bin]] <- model_summary
  
  # === GROUP CONTRASTS ===
  if (bin == "40–60") {
    d_bin$intensity <- factor("i60%", levels = levels(d$intensity))
    ref <- ref_grid(fit, data = d_bin, at = list(intensity = "i60%"))
    emms <- emmeans(ref, ~ group)
    pairwise_contrasts <- pairs(emms)
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    ci_boot$intensity <- "i60%"
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- "i60%"
    emms_df <- as.data.frame(emms)
    emms_df <- tibble::add_column(emms_df, intensity = "i60%", .after = "df")
    emms_df$rt_bin <- bin
    all_emmeans_deltaf[[bin]] <- emms_df
  } else {
    emms <- emmeans(fit, ~ group | intensity, pbkrtest.limit = 15000)
    pairwise_contrasts <- contrast(emms, method = "pairwise")
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- as.data.frame(pairwise_contrasts)$intensity
  }
  
  ci_cols <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot))
  ci_cols <- c(ci_cols[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot))[1])
  stopifnot(length(ci_cols) == 2)
  
  pairwise_contrasts_df <- pairwise_contrasts_df %>%
    left_join(
      ci_boot %>%
        dplyr::select(contrast, intensity, !!sym(ci_cols[1]), !!sym(ci_cols[2])) %>%
        dplyr::rename(boot.LCL = !!sym(ci_cols[1]), boot.UCL = !!sym(ci_cols[2])),
      by = c("contrast", "intensity")
    )
  pairwise_contrasts_df$contrast <- as.character(pairwise_contrasts_df$contrast)
  pairwise_contrasts_df$rt_bin <- bin
  
  sigma_value <- sigma(fit)
  edf_value <- nobs(fit) - length(fixef(fit))
  effect_sizes <- eff_size(emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
  effect_sizes_df <- as.data.frame(effect_sizes)
  effect_sizes_df$contrast <- as.character(effect_sizes_df$contrast)
  if (!"intensity" %in% names(effect_sizes_df)) {
    effect_sizes_df$intensity <- pairwise_contrasts_df$intensity[match(effect_sizes_df$contrast, pairwise_contrasts_df$contrast)]
  }
  ci_low_col <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_df))[1]
  ci_high_col <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_df))[1]
  if (!is.null(ci_low_col) && !is.null(ci_high_col)) {
    pairwise_contrasts_df <- pairwise_contrasts_df %>%
      left_join(
        effect_sizes_df %>%
          dplyr::select(contrast, intensity, effect.size, !!sym(ci_low_col), !!sym(ci_high_col)) %>%
          dplyr::rename(
            hedges_g = effect.size,
            CI_low = !!sym(ci_low_col),
            CI_high = !!sym(ci_high_col)
          ),
        by = c("contrast", "intensity")
      )
  }
  
  all_contrasts_deltaf[[bin]] <- pairwise_contrasts_df
  
  if (bin != "40–60") {
    emms_df <- as.data.frame(emms)
    emms_df$rt_bin <- bin
    all_emmeans_deltaf[[bin]] <- emms_df
  }
  
  # === INTENSITY CONTRASTS ===
  if (bin != "40–60") {
    intensity_emms <- emmeans(fit, ~ intensity | group)
    intensity_contrasts <- contrast(intensity_emms, method = "pairwise")
    ci_boot_int <- as.data.frame(confint(intensity_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    
    ci_cols_int <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot_int))
    ci_cols_int <- c(ci_cols_int[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot_int))[1])
    stopifnot(length(ci_cols_int) == 2)
    
    intensity_contrasts_df <- as.data.frame(intensity_contrasts) %>%
      left_join(
        ci_boot_int %>%
          dplyr::select(contrast, group, !!sym(ci_cols_int[1]), !!sym(ci_cols_int[2])) %>%
          dplyr::rename(boot.LCL = !!sym(ci_cols_int[1]), boot.UCL = !!sym(ci_cols_int[2])),
        by = c("contrast", "group")
      )
    
    intensity_contrasts_df$contrast <- gsub("i", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$contrast <- gsub("%", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$rt_bin <- bin
    intensity_contrasts_df$type <- "intensity_within_group"
    
    effect_sizes_int <- eff_size(intensity_emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
    effect_sizes_int_df <- as.data.frame(effect_sizes_int)
    effect_sizes_int_df$contrast <- gsub("i", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$contrast <- gsub("%", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$rt_bin <- bin
    effect_sizes_int_df$type <- "intensity_within_group"
    
    ci_low_col_int <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_int_df))[1]
    ci_high_col_int <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_int_df))[1]
    
    if (!is.null(ci_low_col_int) && !is.null(ci_high_col_int)) {
      intensity_contrasts_df <- intensity_contrasts_df %>%
        left_join(
          effect_sizes_int_df %>%
            dplyr::select(group, contrast, effect.size, !!sym(ci_low_col_int), !!sym(ci_high_col_int)) %>%
            dplyr::rename(
              hedges_g = effect.size,
              CI_low = !!sym(ci_low_col_int),
              CI_high = !!sym(ci_high_col_int)
            ),
          by = c("group", "contrast")
        )
    }
    
    all_intensity_contrasts_deltaf[[bin]] <- intensity_contrasts_df
  }
}

combined_contrasts_deltaf <- bind_rows(all_contrasts_deltaf) %>%
  filter(!(rt_bin == "40–60" & intensity != "i60%"),
         !(rt_bin == "20–40" & intensity == "i20%"))

if (length(all_intensity_contrasts_deltaf) > 0) {
  combined_intensity_contrasts_deltaf <- bind_rows(all_intensity_contrasts_deltaf)
  combined_intensity_contrasts_deltaf$significant <- with(
    combined_intensity_contrasts_deltaf,
    ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
  )
}

combined_emmeans_deltaf <- bind_rows(all_emmeans_deltaf)
combined_emmeans_deltaf$intensity <- factor(combined_emmeans_deltaf$intensity, levels = c("i20%", "i40%", "i60%"))
combined_summaries_deltaf <- bind_rows(all_summaries_deltaf)
combined_contrasts_deltaf$significant <- with(
  combined_contrasts_deltaf,
  ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
)

print(combined_summaries_deltaf, n = 30)

# Group contrasts - mean difference
plot_deltaf_md <- ggplot(combined_contrasts_deltaf, aes(y = contrast, x = estimate)) +
  geom_point(aes(color = significant), position = position_dodge(width = 1.0), size = 2.5) +
  geom_errorbarh(
    aes(xmin = boot.LCL, xmax = boot.UCL, color = significant),
    height = 0.0, linewidth = 1.0, position = position_dodge(0.0)
  ) +
  facet_grid(
    rt_bin ~ intensity,
    labeller = labeller(rt_bin = c("0–20" = "rt0–20", "20–40" = "rt20–40", "40–60" = "rt40–60"))
  ) +
  scale_color_manual(values = c("yes" = "forestgreen", "no" = "firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_y_discrete(labels = c("20 - 40" = "i20% - i40%", "20 - 60" = "i20% - i60%", "40 - 60" = "i40% - i60%")) +
  labs(x = "ΔF Estimated Mean Difference (pps)", y = "Group Contrasts") +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(color = "black", fill = "gray98"),
    legend.position = "none"
  )

# Group contrasts - effect sizes
plot_deltaf_es <- ggplot(combined_contrasts_deltaf, aes(y = contrast, x = hedges_g)) +
  geom_point(color = "black", size = 2.5) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.0, linewidth = 1, color = "black") +
  facet_grid(
    rt_bin ~ intensity,
    labeller = labeller(rt_bin = c("0–20" = "rt0–20", "20–40" = "rt20–40", "40–60" = "rt40–60"))
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hedges' g (Effect Size)", y = "Group Contrasts") +
  scale_y_discrete(labels = c("20 - 40" = "i20% - i40%", "20 - 60" = "i20% - i60%", "40 - 60" = "i40% - i60%")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(color = "black", fill = "gray98")
  )

# Intensity within-group contrasts - mean difference
plot_deltaf_md_intensity <- ggplot(combined_intensity_contrasts_deltaf, aes(y = contrast, x = estimate)) +
  geom_point(aes(color = significant), size = 2.5) +
  geom_errorbarh(aes(xmin = boot.LCL, xmax = boot.UCL, color = significant), height = 0.0, linewidth = 1) +
  facet_grid(
    rt_bin ~ group,
    labeller = labeller(rt_bin = c("0–20" = "rt0–20", "20–40" = "rt20–40", "40–60" = "rt40–60"))
  ) +
  scale_color_manual(values = c("yes" = "forestgreen", "no" = "firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "ΔF Estimated Mean Difference (pps)", y = "Intensity Contrasts") +
  scale_y_discrete(labels = c("20 - 40" = "i20% - i40%", "20 - 60" = "i20% - i60%", "40 - 60" = "i40% - i60%")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(color = "black", fill = "gray98"),
    legend.position = "none"
  )

# Intensity within-group contrasts - effect sizes
plot_deltaf_es_intensity <- ggplot(combined_intensity_contrasts_deltaf, aes(y = contrast, x = hedges_g)) +
  geom_point(color = "black", size = 2.5) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.0, linewidth = 1, color = "black") +
  facet_grid(
    rt_bin ~ group,
    labeller = labeller(rt_bin = c("0–20" = "rt0–20", "20–40" = "rt20–40", "40–60" = "rt40–60"))
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hedges' g (Effect Size)", y = "Intensity Contrasts") +
  scale_y_discrete(labels = c("20 - 40" = "i20% - i40%", "20 - 60" = "i20% - i60%", "40 - 60" = "i40% - i60%")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(color = "black", fill = "gray98")
  )




# Bin 0–20
d_bin1 <- d %>% filter(rt_bin == "0–20")
d_mean1 <- d_mean %>% filter(rt_bin == "0–20")
mar_bin1 <- combined_emmeans_deltaf %>% filter(rt_bin == "0–20")

# Bin 20–40
d_bin2 <- d %>% filter(rt_bin == "20–40", intensity %in% c("i40%", "i60%"))
d_mean2 <- d_mean %>% filter(rt_bin == "20–40", intensity %in% c("i40%", "i60%"))
mar_bin2 <- combined_emmeans_deltaf %>% filter(rt_bin == "20–40", intensity %in% c("i40%", "i60%"))

# Bin 40–60
d_bin3 <- d %>% filter(rt_bin == "40–60", intensity == "i60%")
d_mean3 <- d_mean %>% filter(rt_bin == "40–60", intensity == "i60%")
mar_bin3 <- combined_emmeans_deltaf %>% filter(rt_bin == "40–60", intensity == "i60%")

# Define plotting function
plot_deltaf_bin <- function(d_raw, d_avg, mar_emmeans, rt_label) {
  ggplot(data = mar_emmeans, aes(x = group, y = emmean)) +
    # MU-level raw points
    geom_quasirandom(
      data = d_raw,
      aes(x = group, y = deltaf),
      color = "gray90",
      width = 0.4,
      dodge.width = 0.3,
      alpha = 0.4,
      size = 1.5,
      shape = 19
    ) +
    # Participant-level averages
    geom_quasirandom(
      data = d_avg,
      aes(x = group, y = deltaf, colour = sex),
      width = 0.1,
      dodge.width = 0.3,
      alpha = 1,
      size = 3,
      shape = 19
    ) +
    # Marginal means ± CI
    geom_errorbar(
      aes(ymin = asymp.LCL, ymax = asymp.UCL),
      position = position_nudge(x = -0.3),
      linewidth = 1,
      width = 0
    ) +
    geom_point(
      position = position_nudge(x = -0.3),
      size = 3,
      colour = "black"
    ) +
    facet_grid(~ factor(intensity, levels = c('i20%', 'i40%', 'i60%'),
                        labels = c('i20%', 'i40%', 'i60%'))) +
    theme_bw(base_size = 14) +
    labs(y = "ΔF (pps)", x = "Groups", title = paste("rt", rt_label)) +
    scale_color_manual(values = c("deeppink", "darkblue")) +
    ylim(-1.0, 12.0) +
    theme(
      axis.text = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 20),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_rect(fill = "gray98", size = 1, color = "black"),
      legend.position = "none"
    )
}

# Generate bin-specific plots
plot1 <- plot_deltaf_bin(d_bin1, d_mean1, mar_bin1, "0–20")
plot2 <- plot_deltaf_bin(d_bin2, d_mean2, mar_bin2, "20–40")
plot3 <- plot_deltaf_bin(d_bin3, d_mean3, mar_bin3, "40–60")

# Layout arrangement
deltaf_top_row <- plot_grid(
  plot2, plot3,
  ncol = 2,
  rel_widths = c(2, 1),
  labels = c("B)", "C)"),
  label_size = 14
)

deltaf_layout <- plot_grid(plot1, deltaf_top_row, ncol = 1, rel_heights = c(1, 1), labels = c("A)"), label_size = 14)

es_layout <- plot_grid(
  plot_deltaf_md,
  plot_deltaf_md_intensity,
  ncol = 1,
  rel_heights = c(1.2, 1),
  labels = c("D)", "E)"),
  label_size = 14
)

final_plot_deltaf <- plot_grid(
  deltaf_layout,
  es_layout,
  ncol = 1,
  rel_heights = c(2.2, 1.2)
)

# Save
ggsave("final_plot_deltaf.png", plot = final_plot_deltaf, width = 12, height = 19, dpi = 300)
ggsave("final_plot_deltaf.eps", plot = final_plot_deltaf, width = 12, height = 19, device = "eps")
ggsave("final_plot_deltaf.eps", plot = final_plot_deltaf, width = 12, height = 19, device = cairo_ps)



################
##BRACE HEIGHT##
################

# Initialize storage for normalized_brace_height
all_contrasts_brace <- list()
all_summaries_brace <- list()
all_emmeans_brace <- list()
all_intensity_contrasts_brace <- list()
all_models_brace <- list()

for (bin in levels(d$rt_bin)) {
  cat("Running model for normalized_brace_height in rt_bin:", bin, "\n")
  
  d_bin <- d %>% filter(rt_bin == bin)
  if (n_distinct(d_bin$group) < 2) {
    cat("  ➔ Skipping: less than 2 groups present\n")
    next
  }
  
  model_formula <- if (bin == "40–60") {
    normalized_brace_height ~ group + sex + (1 | participant / mu_id_intensity)
  } else {
    normalized_brace_height ~ group * intensity + sex + (1 | participant / mu_id_intensity)
  }
  
  emm_options(lmerTest.limit = 15000)
  fit <- rlmer(model_formula, data = d_bin)
  all_models_brace[[bin]] <- fit
  
  model_summary <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  model_summary$rt_bin <- bin
  all_summaries_brace[[bin]] <- model_summary
  
  if (bin == "40–60") {
    d_bin$intensity <- factor("i60%", levels = levels(d$intensity))
    ref <- ref_grid(fit, data = d_bin, at = list(intensity = "i60%"))
    emms <- emmeans(ref, ~ group)
    pairwise_contrasts <- pairs(emms)
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    ci_boot$intensity <- "i60%"
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- "i60%"
    emms_df <- as.data.frame(emms)
    emms_df <- tibble::add_column(emms_df, intensity = "i60%", .after = "df")
    emms_df$rt_bin <- bin
    all_emmeans_brace[[bin]] <- emms_df
  } else {
    emms <- emmeans(fit, ~ group | intensity, pbkrtest.limit = 15000)
    pairwise_contrasts <- contrast(emms, method = "pairwise")
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- as.data.frame(pairwise_contrasts)$intensity
  }
  
  ci_cols <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot))
  ci_cols <- c(ci_cols[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot))[1])
  stopifnot(length(ci_cols) == 2)
  
  pairwise_contrasts_df <- pairwise_contrasts_df %>%
    left_join(
      ci_boot %>%
        dplyr::select(contrast, intensity, !!sym(ci_cols[1]), !!sym(ci_cols[2])) %>%
        dplyr::rename(boot.LCL = !!sym(ci_cols[1]), boot.UCL = !!sym(ci_cols[2])),
      by = c("contrast", "intensity")
    )
  pairwise_contrasts_df$contrast <- as.character(pairwise_contrasts_df$contrast)
  pairwise_contrasts_df$rt_bin <- bin
  
  sigma_value <- sigma(fit)
  edf_value <- nobs(fit) - length(fixef(fit))
  effect_sizes <- eff_size(emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
  effect_sizes_df <- as.data.frame(effect_sizes)
  effect_sizes_df$contrast <- as.character(effect_sizes_df$contrast)
  if (!"intensity" %in% names(effect_sizes_df)) {
    effect_sizes_df$intensity <- pairwise_contrasts_df$intensity[match(effect_sizes_df$contrast, pairwise_contrasts_df$contrast)]
  }
  
  ci_low_col <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_df))[1]
  ci_high_col <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_df))[1]
  if (!is.null(ci_low_col) && !is.null(ci_high_col)) {
    pairwise_contrasts_df <- pairwise_contrasts_df %>%
      left_join(
        effect_sizes_df %>%
          dplyr::select(contrast, intensity, effect.size, !!sym(ci_low_col), !!sym(ci_high_col)) %>%
          dplyr::rename(
            hedges_g = effect.size,
            CI_low = !!sym(ci_low_col),
            CI_high = !!sym(ci_high_col)
          ),
        by = c("contrast", "intensity")
      )
  }
  
  all_contrasts_brace[[bin]] <- pairwise_contrasts_df
  
  if (bin != "40–60") {
    emms_df <- as.data.frame(emms)
    emms_df$rt_bin <- bin
    all_emmeans_brace[[bin]] <- emms_df
  }
  
  if (bin != "40–60") {
    intensity_emms <- emmeans(fit, ~ intensity | group)
    intensity_contrasts <- contrast(intensity_emms, method = "pairwise")
    ci_boot_int <- as.data.frame(confint(intensity_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    
    ci_cols_int <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot_int))
    ci_cols_int <- c(ci_cols_int[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot_int))[1])
    stopifnot(length(ci_cols_int) == 2)
    
    intensity_contrasts_df <- as.data.frame(intensity_contrasts) %>%
      left_join(
        ci_boot_int %>%
          dplyr::select(contrast, group, !!sym(ci_cols_int[1]), !!sym(ci_cols_int[2])) %>%
          dplyr::rename(boot.LCL = !!sym(ci_cols_int[1]), boot.UCL = !!sym(ci_cols_int[2])),
        by = c("contrast", "group")
      )
    
    intensity_contrasts_df$contrast <- gsub("i", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$contrast <- gsub("%", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$rt_bin <- bin
    intensity_contrasts_df$type <- "intensity_within_group"
    
    effect_sizes_int <- eff_size(intensity_emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
    effect_sizes_int_df <- as.data.frame(effect_sizes_int)
    effect_sizes_int_df$contrast <- gsub("i", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$contrast <- gsub("%", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$rt_bin <- bin
    effect_sizes_int_df$type <- "intensity_within_group"
    
    ci_low_col_int <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_int_df))[1]
    ci_high_col_int <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_int_df))[1]
    
    if (!is.null(ci_low_col_int) && !is.null(ci_high_col_int)) {
      intensity_contrasts_df <- intensity_contrasts_df %>%
        left_join(
          effect_sizes_int_df %>%
            dplyr::select(group, contrast, effect.size, !!sym(ci_low_col_int), !!sym(ci_high_col_int)) %>%
            dplyr::rename(
              hedges_g = effect.size,
              CI_low = !!sym(ci_low_col_int),
              CI_high = !!sym(ci_high_col_int)
            ),
          by = c("group", "contrast")
        )
    }
    
    all_intensity_contrasts_brace[[bin]] <- intensity_contrasts_df
  }
}

# Combine and annotate results
combined_contrasts_brace <- bind_rows(all_contrasts_brace) %>%
  filter(!(rt_bin == "40–60" & intensity != "i60%"),
         !(rt_bin == "20–40" & intensity == "i20%"))
combined_contrasts_brace$significant <- with(
  combined_contrasts_brace,
  ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
)

if (length(all_intensity_contrasts_brace) > 0) {
  combined_intensity_contrasts_brace <- bind_rows(all_intensity_contrasts_brace)
  combined_intensity_contrasts_brace$significant <- with(
    combined_intensity_contrasts_brace,
    ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
  )
}

combined_emmeans_brace <- bind_rows(all_emmeans_brace)
combined_emmeans_brace$intensity <- factor(combined_emmeans_brace$intensity, levels = c("i20%", "i40%", "i60%"))
combined_summaries_brace <- bind_rows(all_summaries_brace)

# Optional print
print(combined_summaries_brace, n = 30)


##############
###Peak DR####
##############

# Initialize storage for maximal_discharge_rate_of_the_svr
all_contrasts_mdr <- list()
all_summaries_mdr <- list()
all_emmeans_mdr <- list()
all_intensity_contrasts_mdr <- list()
all_models_mdr <- list()

for (bin in levels(d$rt_bin)) {
  cat("Running model for maximal_discharge_rate_of_the_svr in rt_bin:", bin, "\n")
  
  d_bin <- d %>% filter(rt_bin == bin)
  if (n_distinct(d_bin$group) < 2) {
    cat("  ➔ Skipping: less than 2 groups present\n")
    next
  }
  
  model_formula <- if (bin == "40–60") {
    maximal_discharge_rate_of_the_svr ~ group + sex + (1 | participant / mu_id_intensity)
  } else {
    maximal_discharge_rate_of_the_svr ~ group * intensity + sex + (1 | participant / mu_id_intensity)
  }
  
  emm_options(lmerTest.limit = 15000)
  fit <- rlmer(model_formula, data = d_bin)
  all_models_mdr[[bin]] <- fit
  
  model_summary <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  model_summary$rt_bin <- bin
  all_summaries_mdr[[bin]] <- model_summary
  
  # === GROUP CONTRASTS ===
  if (bin == "40–60") {
    d_bin$intensity <- factor("i60%", levels = levels(d$intensity))
    ref <- ref_grid(fit, data = d_bin, at = list(intensity = "i60%"))
    emms <- emmeans(ref, ~ group)
    pairwise_contrasts <- pairs(emms)
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    ci_boot$intensity <- "i60%"
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- "i60%"
    emms_df <- as.data.frame(emms)
    emms_df <- tibble::add_column(emms_df, intensity = "i60%", .after = "df")
    emms_df$rt_bin <- bin
    all_emmeans_mdr[[bin]] <- emms_df
  } else {
    emms <- emmeans(fit, ~ group | intensity, pbkrtest.limit = 15000)
    pairwise_contrasts <- contrast(emms, method = "pairwise")
    ci_boot <- as.data.frame(confint(pairwise_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    pairwise_contrasts_df <- as.data.frame(pairwise_contrasts)
    pairwise_contrasts_df$intensity <- as.data.frame(pairwise_contrasts)$intensity
  }
  
  ci_cols <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot))
  ci_cols <- c(ci_cols[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot))[1])
  stopifnot(length(ci_cols) == 2)
  
  pairwise_contrasts_df <- pairwise_contrasts_df %>%
    left_join(
      ci_boot %>%
        dplyr::select(contrast, intensity, !!sym(ci_cols[1]), !!sym(ci_cols[2])) %>%
        dplyr::rename(boot.LCL = !!sym(ci_cols[1]), boot.UCL = !!sym(ci_cols[2])),
      by = c("contrast", "intensity")
    )
  pairwise_contrasts_df$contrast <- as.character(pairwise_contrasts_df$contrast)
  pairwise_contrasts_df$rt_bin <- bin
  
  sigma_value <- sigma(fit)
  edf_value <- nobs(fit) - length(fixef(fit))
  effect_sizes <- eff_size(emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
  effect_sizes_df <- as.data.frame(effect_sizes)
  effect_sizes_df$contrast <- as.character(effect_sizes_df$contrast)
  if (!"intensity" %in% names(effect_sizes_df)) {
    effect_sizes_df$intensity <- pairwise_contrasts_df$intensity[match(effect_sizes_df$contrast, pairwise_contrasts_df$contrast)]
  }
  
  ci_low_col <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_df))[1]
  ci_high_col <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_df))[1]
  if (!is.null(ci_low_col) && !is.null(ci_high_col)) {
    pairwise_contrasts_df <- pairwise_contrasts_df %>%
      left_join(
        effect_sizes_df %>%
          dplyr::select(contrast, intensity, effect.size, !!sym(ci_low_col), !!sym(ci_high_col)) %>%
          dplyr::rename(
            hedges_g = effect.size,
            CI_low = !!sym(ci_low_col),
            CI_high = !!sym(ci_high_col)
          ),
        by = c("contrast", "intensity")
      )
  }
  
  all_contrasts_mdr[[bin]] <- pairwise_contrasts_df
  
  if (bin != "40–60") {
    emms_df <- as.data.frame(emms)
    emms_df$rt_bin <- bin
    all_emmeans_mdr[[bin]] <- emms_df
  }
  
  # === INTENSITY CONTRASTS ===
  if (bin != "40–60") {
    intensity_emms <- emmeans(fit, ~ intensity | group)
    intensity_contrasts <- contrast(intensity_emms, method = "pairwise")
    ci_boot_int <- as.data.frame(confint(intensity_contrasts, level = 0.95, bootstrap = TRUE, B = 2000))
    
    ci_cols_int <- intersect(c("lower.CL", "asymp.LCL"), names(ci_boot_int))
    ci_cols_int <- c(ci_cols_int[1], intersect(c("upper.CL", "asymp.UCL"), names(ci_boot_int))[1])
    stopifnot(length(ci_cols_int) == 2)
    
    intensity_contrasts_df <- as.data.frame(intensity_contrasts) %>%
      left_join(
        ci_boot_int %>%
          dplyr::select(contrast, group, !!sym(ci_cols_int[1]), !!sym(ci_cols_int[2])) %>%
          dplyr::rename(boot.LCL = !!sym(ci_cols_int[1]), boot.UCL = !!sym(ci_cols_int[2])),
        by = c("contrast", "group")
      )
    
    intensity_contrasts_df$contrast <- gsub("i", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$contrast <- gsub("%", "", intensity_contrasts_df$contrast)
    intensity_contrasts_df$rt_bin <- bin
    intensity_contrasts_df$type <- "intensity_within_group"
    
    effect_sizes_int <- eff_size(intensity_emms, sigma = sigma_value, edf = edf_value, method = "pairwise", type = "hedges_g")
    effect_sizes_int_df <- as.data.frame(effect_sizes_int)
    effect_sizes_int_df$contrast <- gsub("i", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$contrast <- gsub("%", "", effect_sizes_int_df$contrast)
    effect_sizes_int_df$rt_bin <- bin
    effect_sizes_int_df$type <- "intensity_within_group"
    
    ci_low_col_int <- intersect(c("lower.CL", "asymp.LCL"), names(effect_sizes_int_df))[1]
    ci_high_col_int <- intersect(c("upper.CL", "asymp.UCL"), names(effect_sizes_int_df))[1]
    
    if (!is.null(ci_low_col_int) && !is.null(ci_high_col_int)) {
      intensity_contrasts_df <- intensity_contrasts_df %>%
        left_join(
          effect_sizes_int_df %>%
            dplyr::select(group, contrast, effect.size, !!sym(ci_low_col_int), !!sym(ci_high_col_int)) %>%
            dplyr::rename(
              hedges_g = effect.size,
              CI_low = !!sym(ci_low_col_int),
              CI_high = !!sym(ci_high_col_int)
            ),
          by = c("group", "contrast")
        )
    }
    
    all_intensity_contrasts_mdr[[bin]] <- intensity_contrasts_df
  }
}

# Combine results
combined_contrasts_mdr <- bind_rows(all_contrasts_mdr) %>%
  filter(!(rt_bin == "40–60" & intensity != "i60%"),
         !(rt_bin == "20–40" & intensity == "i20%"))

if (length(all_intensity_contrasts_mdr) > 0) {
  combined_intensity_contrasts_mdr <- bind_rows(all_intensity_contrasts_mdr)
  combined_intensity_contrasts_mdr$significant <- with(
    combined_intensity_contrasts_mdr,
    ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
  )
}

combined_emmeans_mdr <- bind_rows(all_emmeans_mdr)
combined_emmeans_mdr$intensity <- factor(combined_emmeans_mdr$intensity, levels = c("i20%", "i40%", "i60%"))
combined_summaries_mdr <- bind_rows(all_summaries_mdr)
combined_contrasts_mdr$significant <- with(
  combined_contrasts_mdr,
  ifelse(boot.LCL > 0 | boot.UCL < 0, "yes", "no")
)

print(combined_summaries_mdr, n = 30)
