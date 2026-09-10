# ==============================================================================
# PART 2 CRN MODEL
#   Whether cross-life-stage mismatch (adult minus natal) in temperature,
#   rainfall, and population density modulates the carry-over correlation
#   between natal chick weight and adult body mass.
#
# Requires crn_wytham_pairedcontext.stan on STAN_MODEL_PATH
# ==============================================================================

library(cmdstanr)
library(dplyr)
library(posterior)
library(ggplot2)
library(ggdist)
library(tidyr)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load prepared data
# ══════════════════════════════════════════════════════════════════════════════

#Load stan_dl_part1 and context_key created and saved in DATA PREP.R
stan_dl_part1 <- readRDS(file.path(ENTERYOURDIRECTORY, "stan_dl_part1.RDS"))
context_key   <- readRDS(file.path(ENTERYOURDIRECTORY, "context_key.RDS"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Compute mismatch axes (adult - natal)
# ══════════════════════════════════════════════════════════════════════════════
mismatch_temp    <- context_key$adult_temp_ctx - context_key$natal_temp_ctx
mismatch_rain    <- context_key$adult_rain_ctx - context_key$natal_rain_ctx
mismatch_density <- context_key$adult_density  - context_key$natal_density

mismatch_vars <- list(temp = mismatch_temp, rain = mismatch_rain, density = mismatch_density)

cat("\n── Mismatch variable summaries ──\n")
for (nm in names(mismatch_vars)) {
  x <- mismatch_vars[[nm]]
  cat(sprintf("  %-8s mean = %6.3f  sd = %6.3f  range = [%.3f, %.3f]\n",
              nm, mean(x), sd(x), min(x), max(x)))
}

cat("\n── Correlation between mismatch axes ──\n")
print(round(cor(cbind(mismatch_temp, mismatch_rain, mismatch_density)), 3))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Standardize + collinearity/VIF check
# ══════════════════════════════════════════════════════════════════════════════
mismatch_temp_std    <- scale(mismatch_temp)[, 1]
mismatch_rain_std    <- scale(mismatch_rain)[, 1]
mismatch_density_std <- scale(mismatch_density)[, 1]

X_check <- cbind(
  temp_mm = mismatch_temp_std,    temp_mm2 = mismatch_temp_std^2,
  rain_mm = mismatch_rain_std,    rain_mm2 = mismatch_rain_std^2,
  dens_mm = mismatch_density_std, dens_mm2 = mismatch_density_std^2
)

cat("\n── Correlation matrix, mismatch predictor set ──\n")
print(round(cor(X_check), 3))

cat("\n── VIF, mismatch predictor set ──\n")
for (j in seq_len(ncol(X_check))) {
  r2 <- summary(lm(X_check[, j] ~ X_check[, -j]))$r.squared
  cat(sprintf("  %-10s VIF = %.2f\n", colnames(X_check)[j], 1 / (1 - r2)))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Build context matrix Xc: here X4_mismatch (P_y = 7: intercept + 3 mismatch axes, each
# linear + quadratic)
# ══════════════════════════════════════════════════════════════════════════════
X4_mismatch <- as.matrix(cbind(
  intercept   = 1,
  temp_mm     = mismatch_temp_std,
  temp_mm2    = mismatch_temp_std^2,
  rain_mm     = mismatch_rain_std,
  rain_mm2    = mismatch_rain_std^2,
  density_mm  = mismatch_density_std,
  density_mm2 = mismatch_density_std^2
))
rownames(X4_mismatch) <- context_key$context_id

stopifnot(nrow(X4_mismatch) == stan_dl_part1$C, ncol(X4_mismatch) == 7, sum(is.na(X4_mismatch)) == 0)

cat("\n── X4_mismatch ──\n")
cat("Rows:", nrow(X4_mismatch), "  Cols:", ncol(X4_mismatch), "  NAs:", sum(is.na(X4_mismatch)), "\n")

# Parameter coding for B_cpc / B_v (columns of X4_mismatch):
#   [2] temp_mm    [3] temp_mm2
#   [4] rain_mm    [5] rain_mm2
#   [6] density_mm [7] density_mm2

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Assemble stan_dl_mismatch with new context matrix with mismatch terms
# ══════════════════════════════════════════════════════════════════════════════
stan_dl_mismatch     <- stan_dl_part1
stan_dl_mismatch$X   <- X4_mismatch
stan_dl_mismatch$P_y <- ncol(X4_mismatch)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Save prepared model inputs
# ══════════════════════════════════════════════════════════════════════════════
scaling_mismatch <- list(
  temp_mean    = mean(mismatch_temp),    temp_sd    = sd(mismatch_temp),
  rain_mean    = mean(mismatch_rain),    rain_sd    = sd(mismatch_rain),
  density_mean = mean(mismatch_density), density_sd = sd(mismatch_density)
)

mismatch_df <- data.frame(
  context_id       = context_key$context_id,
  mismatch_temp    = mismatch_temp,
  mismatch_rain    = mismatch_rain,
  mismatch_density = mismatch_density
)

saveRDS(stan_dl_mismatch,  file.path(OUTPUT_DIR, "stan_dl_mismatch.RDS"))
saveRDS(X4_mismatch,       file.path(OUTPUT_DIR, "X4_mismatch.RDS"))
saveRDS(scaling_mismatch,  file.path(OUTPUT_DIR, "scaling_mismatch.RDS"))
saveRDS(mismatch_df,       file.path(OUTPUT_DIR, "mismatch_df.RDS"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Compile
# ══════════════════════════════════════════════════════════════════════════════
crn_mod_mismatch <- cmdstan_model(STAN_MODEL_PATH, stanc_options = list("O1"))
cat("Model compiled\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Run Model
# ══════════════════════════════════════════════════════════════════════════════
crn_fit_mismatch_full <- crn_mod_mismatch$sample(
  data            = stan_dl_mismatch,
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
saveRDS(crn_fit_mismatch_full, file.path(OUTPUT_DIR, "crn_fit_mismatch_full.RDS"))

crn_fit_mismatch_full$diagnostic_summary()

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Extract and save posteriors
# ══════════════════════════════════════════════════════════════════════════════
params_mismatch <- c("B_m_chick", "B_m_adult", "B_cpc", "B_v",
                     "sd_season_chick", "sd_season_adult",
                     "sd_breedbox", "sd_brood_chick")

post_draws_mismatch <- crn_fit_mismatch_full$draws(variables = params_mismatch, format = "df")
post_crn_mismatch   <- crn_fit_mismatch_full$draws(variables = c("B_cpc", "B_v"), format = "df")

saveRDS(post_draws_mismatch, file.path(OUTPUT_DIR, "posterior_draws_mismatch.RDS"))
saveRDS(post_crn_mismatch,   file.path(OUTPUT_DIR, "posterior_crn_mismatch.RDS"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 — Posterior summaries
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Part 2 — CRN slopes (B_cpc) ──\n")
slopes_mm <- list(
  list("B_cpc[1,1]", "Correlation intercept (baseline)"),
  list("B_cpc[2,1]", "Correlation x thermal mismatch (linear)"),
  list("B_cpc[3,1]", "Correlation x thermal mismatch2 (quadratic)"),
  list("B_cpc[4,1]", "Correlation x rain mismatch (linear)"),
  list("B_cpc[5,1]", "Correlation x rain mismatch2 (quadratic)"),
  list("B_cpc[6,1]", "Correlation x density mismatch (linear)"),
  list("B_cpc[7,1]", "Correlation x density mismatch2 (quadratic)")
)
for (s in slopes_mm) {
  draws <- post_crn_mismatch[[s[[1]]]]
  cat(sprintf("%-44s median = %6.3f  [%6.3f, %6.3f]  P(direction) = %.3f\n",
              s[[2]], median(draws), quantile(draws, 0.10), quantile(draws, 0.90),
              max(mean(draws > 0), mean(draws < 0))))
}

cat("\n── Part 2 — variance reaction norm slopes (B_v) ──\n")
slopes_v_mm <- list(
  list("B_v[1,1]", "Chick wt variance intercept"),
  list("B_v[2,1]", "Chick wt variance x thermal mismatch"),
  list("B_v[3,1]", "Chick wt variance x thermal mismatch2"),
  list("B_v[4,1]", "Chick wt variance x rain mismatch"),
  list("B_v[5,1]", "Chick wt variance x rain mismatch2"),
  list("B_v[6,1]", "Chick wt variance x density mismatch"),
  list("B_v[7,1]", "Chick wt variance x density mismatch2"),
  list("B_v[1,2]", "Adult wt variance intercept"),
  list("B_v[2,2]", "Adult wt variance x thermal mismatch"),
  list("B_v[3,2]", "Adult wt variance x thermal mismatch2"),
  list("B_v[4,2]", "Adult wt variance x rain mismatch"),
  list("B_v[5,2]", "Adult wt variance x rain mismatch2"),
  list("B_v[6,2]", "Adult wt variance x density mismatch"),
  list("B_v[7,2]", "Adult wt variance x density mismatch2")
)
for (s in slopes_v_mm) {
  draws <- post_crn_mismatch[[s[[1]]]]
  cat(sprintf("%-40s median = %6.3f  [%6.3f, %6.3f]\n",
              s[[2]], median(draws), quantile(draws, 0.10), quantile(draws, 0.90)))
}

cat("\n── Fixed effects on chick weight (B_m_chick) ──\n")
fixed_chick_mm <- lapply(seq_len(stan_dl_mismatch$P_chick), function(j) {
  draws <- post_draws_mismatch[[paste0("B_m_chick[", j, "]")]]
  c(median = median(draws), lo = quantile(draws, 0.10), hi = quantile(draws, 0.90))
})
names(fixed_chick_mm) <- colnames(stan_dl_mismatch$X_chick)
print(fixed_chick_mm)

cat("\n── Fixed effects on adult weight (B_m_adult) ──\n")
fixed_adult_mm <- lapply(seq_len(stan_dl_mismatch$P_adult), function(j) {
  draws <- post_draws_mismatch[[paste0("B_m_adult[", j, "]")]]
  c(median = median(draws), lo = quantile(draws, 0.10), hi = quantile(draws, 0.90))
})
names(fixed_adult_mm) <- colnames(stan_dl_mismatch$X_adult)
print(fixed_adult_mm)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 11 - Plots to visualise CRN slopes
# ══════════════════════════════════════════════════════════════════════════════

#Load what is needed for making plots if running on it's own

mismatch_df       <- readRDS(file.path(YOURDIRECTORY, "mismatch_df.RDS"))
scaling_mismatch  <- readRDS(file.path(YOURDIRECTORY, "scaling_mismatch.RDS"))
post_crn_mismatch <- readRDS(file.path(YOURDIRECTORY, "posterior_crn_mismatch.RDS"))


B_cpc <- post_crn_mismatch %>% select(starts_with("B_cpc[")) %>% as.matrix()
stopifnot(ncol(B_cpc) == 7)

# Predictor variables specifications ─────────────────────────
mismatch_specs <- list(
  list(id = "temp_mismatch",    raw_col = "mismatch_temp",    b_idx = 2, quad_idx = 3,
       mean_key = "temp_mean",    sd_key = "temp_sd",
       xlab = "Thermal mismatch (\u00b0C)\nadult - natal temperature",
       colour = "#ea5f94", fill = "#ea5f94"),
  list(id = "rain_mismatch",    raw_col = "mismatch_rain",    b_idx = 4, quad_idx = NA,
       mean_key = "rain_mean",    sd_key = "rain_sd",
       xlab = "Rainfall mismatch (mm)\nadult - natal rainfall",
       colour = "#2A7B9B", fill = "#3FA8CF"),
  list(id = "density_mismatch", raw_col = "mismatch_density", b_idx = 6, quad_idx = NA,
       mean_key = "density_mean", sd_key = "density_sd",
       xlab = "Density mismatch (pairs)\nadult - natal density",
       colour = "#3B6D11", fill = "#97C459")
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
    x   = x_real,    median = pred_mat[, "median"],
    q05 = pred_mat[, "q05.5%"],  q95 = pred_mat[, "q95.95%"],
    q10 = pred_mat[, "q10.10%"], q90 = pred_mat[, "q90.90%"],
    q25 = pred_mat[, "q25.25%"], q75 = pred_mat[, "q75.75%"]
  )
}

mismatch_ribbon_plot <- function(df, xlab, colour, fill, rug_x,
                                 ylab = "Among-individual correlation") {
  ggplot(df, aes(x = x)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted",  colour = "grey50", linewidth = 0.6) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = fill, alpha = 0.15) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = fill, alpha = 0.20) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = fill, alpha = 0.30) +
    geom_line(aes(y = median), colour = colour, linewidth = 1.2, lineend = "round") +
    geom_rug(data = data.frame(x = rug_x), aes(x = x), colour = colour, alpha = 0.4,
             linewidth = 0.3, sides = "b", length = unit(0.03, "npc")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    theme(
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 14),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}

# ── CRN plots ─────────────────────────
ribbon_plots_mm <- list()
for (spec in mismatch_specs) {
  std_seq  <- make_std_seq(mismatch_df[[spec$raw_col]])
  real_seq <- to_real(std_seq, scaling_mismatch[[spec$mean_key]], scaling_mismatch[[spec$sd_key]])
  pred_df  <- to_df(real_seq, make_pred(std_seq, spec$b_idx, spec$quad_idx))
  
  ribbon_plots_mm[[spec$id]] <- mismatch_ribbon_plot(pred_df, spec$xlab, spec$colour, spec$fill,
                                                     rug_x = mismatch_df[[spec$raw_col]])
}

# ── posterior density plots ─────────────────────────
mm_slopes_long <- data.frame(
  `Temperature mismatch (linear)`    = post_crn_mismatch$`B_cpc[2,1]`,
  `Temperature mismatch (quadratic)` = post_crn_mismatch$`B_cpc[3,1]`,
  `Rainfall mismatch (linear)`       = post_crn_mismatch$`B_cpc[4,1]`,
  `Rainfall mismatch (quadratic)`    = post_crn_mismatch$`B_cpc[5,1]`,
  `Density mismatch (linear)`        = post_crn_mismatch$`B_cpc[6,1]`,
  `Density mismatch (quadratic)`     = post_crn_mismatch$`B_cpc[7,1]`,
  check.names = FALSE
) %>%
  pivot_longer(everything(), names_to = "label", values_to = "value") %>%
  mutate(axis = case_when(
    grepl("Temperature", label) ~ "Temperature",
    grepl("Rainfall",    label) ~ "Rainfall",
    grepl("Density",     label) ~ "Density"
  ))

mm_label_order <- c(
  "Rainfall mismatch (quadratic)",    "Rainfall mismatch (linear)",
  "Temperature mismatch (quadratic)", "Temperature mismatch (linear)",
  "Density mismatch (quadratic)",     "Density mismatch (linear)"
)
mm_slopes_long$label <- factor(mm_slopes_long$label, levels = mm_label_order)
mm_slopes_long$axis  <- factor(mm_slopes_long$axis,  levels = c("Temperature", "Rainfall", "Density"))

mm_axis_colours <- c(Temperature = "#ea5f94", Rainfall = "#3FA8CF", Density = "#97C459")

p_halfeye_mm_combined <- ggplot(mm_slopes_long, aes(x = value, y = label, fill = axis, colour = axis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  stat_halfeye(.width = c(0.50, 0.90), point_colour = "black", point_size = 2.2, stroke = 1.2,
               slab_alpha = 0.55, normalize = "panels") +
  scale_fill_manual(values = mm_axis_colours) +
  scale_colour_manual(values = mm_axis_colours) +
  labs(x = "Posterior effect size (atanh scale)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y      = element_text(size = 13, colour = "grey20"),
    axis.ticks.y     = element_blank(), axis.line.y = element_blank(),
    axis.text.x      = element_text(size = 12.5, colour = "grey20"),
    axis.title.x     = element_text(size = 13.5, colour = "grey20"),
    legend.position  = "none",
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.margin      = margin(15, 15, 10, 15)
  )

# ── Print each panel ───────────────────────────────────────────────────────────
print(ribbon_plots_mm$temp_mismatch)
print(ribbon_plots_mm$rain_mismatch)
print(ribbon_plots_mm$density_mismatch)
print(p_halfeye_mm_combined)
