# ==============================================================================
# PART 1 CRN MODEL
#   How temperature (natal and adult, quadratic), rainfall (natal and adult,
#   linear), and population density (natal and adult, linear) modulate the
#   carry-over correlation between natal chick weight and adult body mass.
#
# Requires crn_wytham_pairedcontext.stan on STAN_MODEL_PATH
# ==============================================================================

library(cmdstanr)
library(dplyr)
library(posterior)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load prepared data
# ══════════════════════════════════════════════════════════════════════════════

#Load stan_dl_part1 created and saved in DATA PREP.R
stan_dl_part1 <- readRDS("ENTERYOURDIRECTORY/stan_dl_part1.RDS")

stopifnot(stan_dl_part1$P_y == 9)
cat("Loaded stan_dl_part1 — N:", stan_dl_part1$N, " M:", stan_dl_part1$M,
    " C:", stan_dl_part1$C, " P_y:", stan_dl_part1$P_y, "\n")

# Parameter coding for B_cpc / B_v (columns of X):
#   [2] natal_temp    [3] natal_temp2
#   [4] adult_temp    [5] adult_temp2
#   [6] natal_rain    [7] adult_rain
#   [8] natal_density [9] adult_density

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Compile
# ══════════════════════════════════════════════════════════════════════════════
crn_mod_part1 <- cmdstan_model(STAN_MODEL_PATH, stanc_options = list("O1"))
cat("Model compiled\n")

dir.create(OUTPUT_DIR, showWarnings = FALSE) #create an output directory to save the files as model runs
# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Run model
# ══════════════════════════════════════════════════════════════════════════════
crn_fit_part1_full <- crn_mod_part1$sample(
  data            = stan_dl_part1,
  output_dir      = OUTPUT_DIR, 
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 3,
  parallel_chains = 3,
  adapt_delta     = 0.95,
  max_treedepth   = 15,
  init            = 0.01,
  refresh         = 50,
  seed            = 1234
)

saveRDS(crn_fit_part1_full, "crn_fit_part1_full.RDS")
cat("Full run saved\n")

crn_fit_part1_full$diagnostic_summary()

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Extract and save posteriors
# ══════════════════════════════════════════════════════════════════════════════
params_part1 <- c("B_m_chick", "B_m_adult", "B_cpc", "B_v",
                  "sd_season_chick", "sd_season_adult",
                  "sd_breedbox", "sd_brood_chick")

post_draws_part1 <- crn_fit_part1_full$draws(variables = params_part1, format = "df")
post_crn_part1   <- crn_fit_part1_full$draws(variables = c("B_cpc", "B_v"), format = "df")

saveRDS(post_draws_part1, "posterior_draws_part1.RDS")
saveRDS(post_crn_part1, "posterior_crn_part1.RDS")

full_summary_part1 <- crn_fit_part1_full$summary(variables = params_part1)
saveRDS(full_summary_part1, "parameter_summary_part1.RDS")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Posterior summaries
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Part 1 — CRN slopes (B_cpc) ──\n")
slopes_part1 <- list(
  list("B_cpc[1,1]", "Correlation intercept (baseline)"),
  list("B_cpc[2,1]", "Correlation x natal_temp"),
  list("B_cpc[3,1]", "Correlation x natal_temp2"),
  list("B_cpc[4,1]", "Correlation x adult_temp"),
  list("B_cpc[5,1]", "Correlation x adult_temp2"),
  list("B_cpc[6,1]", "Correlation x natal_rain"),
  list("B_cpc[7,1]", "Correlation x adult_rain"),
  list("B_cpc[8,1]", "Correlation x natal_density"),
  list("B_cpc[9,1]", "Correlation x adult_density")
)
for (s in slopes_part1) {
  draws <- post_crn_part1[[s[[1]]]]
  cat(sprintf("%-40s median = %6.3f  [%6.3f, %6.3f]  P(direction) = %.3f\n",
              s[[2]], median(draws), quantile(draws, 0.10), quantile(draws, 0.90),
              max(mean(draws > 0), mean(draws < 0))))
}

cat("\n── Part 1 — variance reaction norm slopes (B_v) ──\n")
slopes_v <- list(
  list("B_v[1,1]", "Chick wt variance intercept"),
  list("B_v[2,1]", "Chick wt variance x natal_temp"),
  list("B_v[3,1]", "Chick wt variance x natal_temp2"),
  list("B_v[4,1]", "Chick wt variance x adult_temp"),
  list("B_v[5,1]", "Chick wt variance x adult_temp2"),
  list("B_v[6,1]", "Chick wt variance x natal_rain"),
  list("B_v[7,1]", "Chick wt variance x adult_rain"),
  list("B_v[8,1]", "Chick wt variance x natal_density"),
  list("B_v[9,1]", "Chick wt variance x adult_density"),
  list("B_v[1,2]", "Adult wt variance intercept"),
  list("B_v[2,2]", "Adult wt variance x natal_temp"),
  list("B_v[3,2]", "Adult wt variance x natal_temp2"),
  list("B_v[4,2]", "Adult wt variance x adult_temp"),
  list("B_v[5,2]", "Adult wt variance x adult_temp2"),
  list("B_v[6,2]", "Adult wt variance x natal_rain"),
  list("B_v[7,2]", "Adult wt variance x adult_rain"),
  list("B_v[8,2]", "Adult wt variance x natal_density"),
  list("B_v[9,2]", "Adult wt variance x adult_density")
)
for (s in slopes_v) {
  draws <- post_crn_part1[[s[[1]]]]
  cat(sprintf("%-38s median = %6.3f  [%6.3f, %6.3f]\n",
              s[[2]], median(draws), quantile(draws, 0.10), quantile(draws, 0.90)))
}

cat("\n── Fixed effects on chick weight (B_m_chick) ──\n")
fixed_chick <- lapply(seq_len(stan_dl_part1$P_chick), function(j) {
  draws <- post_draws_part1[[paste0("B_m_chick[", j, "]")]]
  c(median = median(draws), lo = quantile(draws, 0.10), hi = quantile(draws, 0.90))
})
names(fixed_chick) <- colnames(stan_dl_part1$X_chick)
print(fixed_chick)

cat("\n── Fixed effects on adult weight (B_m_adult) ──\n")
fixed_adult <- lapply(seq_len(stan_dl_part1$P_adult), function(j) {
  draws <- post_draws_part1[[paste0("B_m_adult[", j, "]")]]
  c(median = median(draws), lo = quantile(draws, 0.10), hi = quantile(draws, 0.90))
})
names(fixed_adult) <- colnames(stan_dl_part1$X_adult)
print(fixed_adult)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Plots to visualise CRN slopes
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(ggdist)
library(patchwork)

#Load what is needed for making plots if running on it's own

context_key <- readRDS(file.path(YOURDIRECTORY, "context_key.RDS")) #from DATA PREP.R
scaling_part1 <- readRDS(file.path(YOURDIRECTORY, "scaling_part1.RDS")) #from DATA PREP.R
post_crn_part1 <- readRDS(file.path(YOURDIRECTORY, "posterior_crn_part1.RDS")) #from step 4

B_cpc <- post_crn_part1 %>% select(starts_with("B_cpc[")) %>% as.matrix()

# environmental variable x life stage, b_idx / quad_idx are column positions in B_cpc
#quad_idx = NA is a linear-only predictor
predictor_specs <- list(
  list(id = "natal_density", raw_col = "natal_density",  b_idx = 8, quad_idx = NA,
       mean_key = "natal_dens_mean", sd_key = "natal_dens_sd",
       xlab = "Natal density (pairs)",              colour = "#3B6D11", fill = "#97C459", family = "density"),
  list(id = "adult_density", raw_col = "adult_density",  b_idx = 9, quad_idx = NA,
       mean_key = "adult_dens_mean", sd_key = "adult_dens_sd",
       xlab = "Adult density (pairs)",               colour = "#3B6D11", fill = "#97C459", family = "density"),
  list(id = "natal_temp",    raw_col = "natal_temp_ctx", b_idx = 2, quad_idx = 3,
       mean_key = "natal_temp_mean", sd_key = "natal_temp_sd",
       xlab = "Natal temperature (\u00b0C)",          colour = "#ea5f94", fill = "#ea5f94", family = "temp"),
  list(id = "adult_temp",    raw_col = "adult_temp_ctx", b_idx = 4, quad_idx = 5,
       mean_key = "adult_temp_mean", sd_key = "adult_temp_sd",
       xlab = "Adult breeding temperature (\u00b0C)", colour = "#ea5f94", fill = "#ea5f94", family = "temp"),
  list(id = "natal_rain",    raw_col = "natal_rain_ctx", b_idx = 6, quad_idx = NA,
       mean_key = "natal_rain_mean", sd_key = "natal_rain_sd",
       xlab = "Natal rainfall (mm)",                  colour = "#2A7B9B", fill = "#3FA8CF", family = "rain"),
  list(id = "adult_rain",    raw_col = "adult_rain_ctx", b_idx = 7, quad_idx = NA,
       mean_key = "adult_rain_mean", sd_key = "adult_rain_sd",
       xlab = "Adult breeding rainfall (mm)",         colour = "#2A7B9B", fill = "#3FA8CF", family = "rain")
)

#helper functions for plots
N_GRID <- 80

make_std_seq <- function(x_raw) seq(min(scale(x_raw)), max(scale(x_raw)), length.out = N_GRID)
to_real      <- function(x_std, mean_val, sd_val) x_std * sd_val + mean_val

make_contrast <- function(t, b_idx, quad_idx, P_y) {
  v <- rep(0, P_y); v[1] <- 1; v[b_idx] <- t
  if (!is.na(quad_idx)) v[quad_idx] <- t^2
  v
}

make_pred <- function(seq_std, b_idx, quad_idx) {
  t(sapply(seq_std, function(t) {
    draws <- tanh(B_cpc %*% make_contrast(t, b_idx, quad_idx, ncol(B_cpc)))
    c(median = median(draws),
      q05 = quantile(draws, 0.05), q95 = quantile(draws, 0.95),
      q10 = quantile(draws, 0.10), q90 = quantile(draws, 0.90),
      q25 = quantile(draws, 0.25), q75 = quantile(draws, 0.75))
  }))
}

to_df <- function(x_real, pred_mat) {
  data.frame(
    x      = x_real,      median = pred_mat[, "median"],
    q05    = pred_mat[, "q05.5%"],  q95 = pred_mat[, "q95.95%"],
    q10    = pred_mat[, "q10.10%"], q90 = pred_mat[, "q90.90%"],
    q25    = pred_mat[, "q25.25%"], q75 = pred_mat[, "q75.75%"]
  )
}

#crn plot function
crn_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold")
  )

ribbon_plot <- function(df, xlab, colour, fill, rug_x,
                        ylab = "Among-individual correlation", ylim = c(0, 1)) {
  ggplot(df, aes(x = x)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = fill, alpha = 0.15) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = fill, alpha = 0.20) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = fill, alpha = 0.30) +
    geom_line(aes(y = median), colour = colour, linewidth = 1.2, lineend = "round") +
    geom_rug(data = data.frame(x = rug_x), aes(x = x), colour = colour, alpha = 0.4,
             linewidth = 0.3, sides = "b", length = unit(0.03, "npc")) +
    scale_y_continuous(limits = ylim, breaks = seq(0, 1, 0.2)) +
    labs(x = xlab, y = ylab) +
    crn_theme
}

#posterior density plot function
halfeye_theme <- theme_classic() +
  theme(
    axis.text.y      = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.text.x      = element_text(size = 11, colour = "grey20"),
    axis.title.x     = element_text(size = 12, colour = "grey20"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.margin      = margin(15, 15, 10, 15)
  )

make_single_halfeye <- function(values, colour, fill, xlab = "Posterior effect size") {
  ggplot(data.frame(value = values), aes(x = value, y = 1)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    stat_halfeye(.width = c(0.50, 0.90), point_colour = "black", point_size = 2.2, stroke = 1.2,
                 fill = fill, colour = colour, slab_alpha = 0.55, normalize = "panels") +
    scale_y_continuous(expand = expansion(add = c(0.04, 0.05))) +
    labs(x = xlab, y = NULL) +
    halfeye_theme
}

make_combined_halfeye <- function(linear_vals, quad_vals, colour, fill,
                                  xlab = "Posterior effect size") {
  df <- data.frame(
    value = c(linear_vals, quad_vals),
    term  = factor(rep(c("Linear", "Quadratic"), each = length(linear_vals)),
                   levels = c("Quadratic", "Linear"))
  )
  ggplot(df, aes(x = value, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    stat_halfeye(.width = c(0.50, 0.90), point_colour = "black", point_size = 2.2, stroke = 1.2,
                 fill = fill, colour = colour, slab_alpha = 0.55, normalize = "panels") +
    scale_y_discrete(expand = expansion(add = c(0.15, 0.55))) +
    labs(x = xlab, y = NULL) +
    halfeye_theme +
    theme(axis.text.y = element_text(size = 12, colour = "grey20"))
}

#__________________________CRN PLOTS FOR MAIN TEXT FIG 2 __________________________

ribbon_plots <- list()
for (spec in predictor_specs) {
  std_seq  <- make_std_seq(context_key[[spec$raw_col]])
  real_seq <- to_real(std_seq, scaling_part1[[spec$mean_key]], scaling_part1[[spec$sd_key]])
  pred_df  <- to_df(real_seq, make_pred(std_seq, spec$b_idx, spec$quad_idx))
  
  ribbon_plots[[spec$id]] <- ribbon_plot(pred_df, spec$xlab, spec$colour, spec$fill,
                                         rug_x = context_key[[spec$raw_col]])
}

# Adult | natal combined row, one per environmental variable
ribbon_rows <- list()
families <- unique(vapply(predictor_specs, `[[`, character(1), "family"))
for (fam in families) {
  ids      <- vapply(Filter(function(s) s$family == fam, predictor_specs), `[[`, character(1), "id")
  adult_id <- grep("^adult_", ids, value = TRUE)
  natal_id <- grep("^natal_", ids, value = TRUE)
  ribbon_rows[[fam]] <- ribbon_plots[[adult_id]] | ribbon_plots[[natal_id]]
}

# ── Print each ribbon panel ───────────────────────────────────────────────────
print(ribbon_plots$natal_density)
print(ribbon_plots$adult_density)
print(ribbon_plots$natal_temp)
print(ribbon_plots$adult_temp)
print(ribbon_plots$natal_rain)
print(ribbon_plots$adult_rain)

print(ribbon_rows$density)
print(ribbon_rows$temp)
print(ribbon_rows$rain)

#_______________POSTERIOR DENSITY PLOTS FOR MAIN TEXT FIG 2_______________

halfeyes <- list()
for (spec in predictor_specs) {
  values <- B_cpc[, spec$b_idx]
  
  if (is.na(spec$quad_idx)) {
    halfeyes[[spec$id]] <- make_single_halfeye(values, spec$colour, spec$fill)
  } else {
    quad_values <- B_cpc[, spec$quad_idx]
    halfeyes[[paste0(spec$id, "_linear")]]    <- make_single_halfeye(values, spec$colour, spec$fill)
    halfeyes[[paste0(spec$id, "_quadratic")]] <- make_single_halfeye(quad_values, spec$colour, spec$fill)
    halfeyes[[paste0(spec$id, "_combined")]]  <- make_combined_halfeye(values, quad_values, spec$colour, spec$fill)
  }
}

# ── Print each half-eye panel ─────────────────────────────────────────────────
print(halfeyes$natal_density)
print(halfeyes$adult_density)
print(halfeyes$natal_temp_linear)
print(halfeyes$natal_temp_quadratic)
print(halfeyes$natal_temp_combined)
print(halfeyes$adult_temp_linear)
print(halfeyes$adult_temp_quadratic)
print(halfeyes$adult_temp_combined)
print(halfeyes$natal_rain)
print(halfeyes$adult_rain)

