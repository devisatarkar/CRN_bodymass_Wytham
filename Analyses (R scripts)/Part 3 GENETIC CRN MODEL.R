# ==============================================================================
# PART 3 Genetic CRN MODEL
#   Bivariate animal model decomposing the phenotypic carry-over correlation
#   into an additive genetic component (r_A) and a non-genetic individual
#   component (r_E), to see how each component varies across temperature 
#   and population density during natal and adult life stages.
#
# Requires crn_wytham_genetic+nongenetic.stan on STAN_MODEL_PATH
# ==============================================================================

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Preparing the relatedness matrix
# ══════════════════════════════════════════════════════════════════════════════

library(nadiv)
library(Matrix)
library(dplyr)

# Load the pruned pedigree and the individual order to match
# ══════════════════════════════════════════════════════════════════════════════
ped_pruned <- readRDS(file.path(ENTERYOURDIRECTORY, "ped_pruned.RDS"))
df_chick_p <- readRDS(file.path(ENTERYOURDIRECTORY, "df_chick_p.RDS"))

# Build the pedigree index: row in ped_pruned for each id_num 1:I
# ══════════════════════════════════════════════════════════════════════════════
id_order <- df_chick_p$id
ped_idx  <- match(id_order, ped_pruned$id)

stopifnot(sum(is.na(ped_idx)) == 0)
stopifnot(min(ped_idx) >= 1, max(ped_idx) <= nrow(ped_pruned))

cat("\n── Pedigree index ──\n")
cat("NAs:", sum(is.na(ped_idx)), "— expected 0\n")
cat("Range:", min(ped_idx), "to", max(ped_idx), "\n")

# Spot check a few individuals map back to the ring number expected
for (i in c(1, round(nrow(df_chick_p) / 2), nrow(df_chick_p))) {
  ring    <- df_chick_p$id[i]
  id_back <- as.character(ped_pruned$id[ped_idx[i]])
  cat("id_num", i, "ring:", ring, "— match:", ring == id_back, "\n")
}

# Compute A (dense numerator relationship matrix) 
# ══════════════════════════════════════════════════════════════════════════════

A_full <- makeA(ped_pruned)
cat("A_full dimensions:", nrow(A_full), "x", ncol(A_full), "\n")

A_matrix <- as.matrix(A_full[ped_idx, ped_idx])

stopifnot(nrow(A_matrix) == nrow(df_chick_p), ncol(A_matrix) == nrow(df_chick_p))

cat("A_matrix dimensions:", nrow(A_matrix), "x", ncol(A_matrix), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Load data for CRN model
# ══════════════════════════════════════════════════════════════════════════════

library(cmdstanr)
library(dplyr)
library(posterior)

#Load stan_dl_part1 and context_key created and saved in DATA PREP.R
stan_dl_part1 <- readRDS(file.path(ENTERYOURDIRECTORY, "stan_dl_part1.RDS"))
context_key   <- readRDS(file.path(ENTERYOURDIRECTORY, "context_key.RDS"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Build context matrix Xc: here X4_genetic (P_y = 5)
# ══════════════════════════════════════════════════════════════════════════════
natal_temp_z    <- scale(context_key$natal_temp_ctx)[, 1]
adult_temp_z    <- scale(context_key$adult_temp_ctx)[, 1]
natal_density_z <- scale(context_key$natal_density)[, 1]
adult_density_z <- scale(context_key$adult_density)[, 1]

X4_genetic <- as.matrix(cbind(
  intercept       = 1,
  natal_temp_z    = natal_temp_z,
  adult_temp_z    = adult_temp_z,
  natal_density_z = natal_density_z,
  adult_density_z = adult_density_z
))
rownames(X4_genetic) <- context_key$context_id

cat("\n── VIF check ──\n")
X_check <- X4_genetic[, -1]
for (j in seq_len(ncol(X_check))) {
  r2 <- summary(lm(X_check[, j] ~ X_check[, -j]))$r.squared
  cat(sprintf("  %-16s VIF = %.2f\n", colnames(X_check)[j], 1 / (1 - r2)))
}
cat("\n── X4_genetic ──\n")
cat("Rows:", nrow(X4_genetic), "  Cols:", ncol(X4_genetic), "  NAs:", sum(is.na(X4_genetic)), "\n")

# Parameter coding for B_cpc_A / B_cpc_E / B_v_A / B_v_E (columns of
# X4_genetic):
#   [1] intercept  [2] natal_temp
#   [3] adult_temp   [4] natal_density 
#   [5] adult_density 

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Assemble stan_dl_genetic
# Everything except X / P_y is from stan_dl_part1; A is new here
# ══════════════════════════════════════════════════════════════════════════════
stan_dl_genetic     <- stan_dl_part1
stan_dl_genetic$X   <- X4_genetic
stan_dl_genetic$P_y <- ncol(X4_genetic)
stan_dl_genetic$A   <- A_matrix

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Compile Stan model
# ══════════════════════════════════════════════════════════════════════════════
crn_mod_genetic <- cmdstan_model(STAN_MODEL_PATH, stanc_options = list("O1"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Run model, 9,000 draws (3 chains x 3,000 sampling)
# ══════════════════════════════════════════════════════════════════════════════

crn_fit_genetic_full <- crn_mod_genetic$sample(
  data            = stan_dl_genetic,
  output_dir      = OUTPUT_DIR,
  iter_warmup     = 1500,
  iter_sampling   = 3000,
  chains          = 3,
  parallel_chains = 3,
  adapt_delta     = 0.97,
  max_treedepth   = 15,
  init            = 0.01,
  refresh         = 50,
  seed            = 1234
)
saveRDS(crn_fit_genetic_full, file.path(OUTPUT_DIR, "crn_fit_genetic_full.RDS"))

crn_fit_genetic_full$diagnostic_summary()
crn_fit_genetic_full$summary(
  variables = c("B_cpc_A", "B_v_A", "B_cpc_E", "B_v_E", "B_m_chick", "B_m_adult",
                "sd_season_chick", "sd_season_adult", "sd_breedbox", "sd_brood_chick")
) %>% print(n = 60)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Extract and save posteriors
# ══════════════════════════════════════════════════════════════════════════════
params_genetic <- c("B_cpc_A", "B_v_A", "B_cpc_E", "B_v_E", "B_m_chick", "B_m_adult",
                    "sd_season_chick", "sd_season_adult", "sd_breedbox", "sd_brood_chick")

post_draws_genetic <- crn_fit_genetic_full$draws(variables = params_genetic, format = "df")
post_crn_genetic   <- crn_fit_genetic_full$draws(
  variables = c("B_cpc_A", "B_v_A", "B_cpc_E", "B_v_E"), format = "df")

saveRDS(post_draws_genetic, file.path(YOURDIRECTORY, "posterior_draws_genetic.RDS"))
saveRDS(post_crn_genetic,   file.path(YOURDIRECTORY, "posterior_crn_genetic.RDS"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Posterior summaries
# ══════════════════════════════════════════════════════════════════════════════

summarise_param <- function(draws_df, param, label, width = 30) {
  d <- draws_df[[param]]
  cat(sprintf(paste0("%-", width, "s median = %7.3f  [%7.3f, %7.3f]\n"),
              label, median(d), quantile(d, 0.10), quantile(d, 0.90)))
}

predictor_labels <- c("Intercept (baseline, atanh)", "Natal temp", "Adult temp",
                      "Natal density", "Adult density")

print_reaction_norm <- function(draws_df, prefix, col = 1) {
  for (j in seq_along(predictor_labels)) {
    summarise_param(draws_df, sprintf("%s[%d,%d]", prefix, j, col), predictor_labels[j])
  }
}

cat("\n== r_A (genetic correlation) reaction norm ==\n")
print_reaction_norm(post_draws_genetic, "B_cpc_A")
cat("\n== r_E (non-genetic correlation) reaction norm ==\n")
print_reaction_norm(post_draws_genetic, "B_cpc_E")
cat("\n== VA (genetic variance) reaction norm - chick weight ==\n")
print_reaction_norm(post_draws_genetic, "B_v_A", col = 1)
cat("\n== VA (genetic variance) reaction norm - adult weight ==\n")
print_reaction_norm(post_draws_genetic, "B_v_A", col = 2)
cat("\n== VE (non-genetic variance) reaction norm - chick weight ==\n")
print_reaction_norm(post_draws_genetic, "B_v_E", col = 1)
cat("\n== VE (non-genetic variance) reaction norm - adult weight ==\n")
print_reaction_norm(post_draws_genetic, "B_v_E", col = 2)

cat("\n== Fixed effects - chick weight ==\n")
for (j in seq_len(stan_dl_genetic$P_chick)) {
  summarise_param(post_draws_genetic, paste0("B_m_chick[", j, "]"),
                  colnames(stan_dl_genetic$X_chick)[j])
}
cat("\n== Fixed effects - adult weight ==\n")
for (j in seq_len(stan_dl_genetic$P_adult)) {
  summarise_param(post_draws_genetic, paste0("B_m_adult[", j, "]"),
                  colnames(stan_dl_genetic$X_adult)[j])
}

# ── Variance decomposition and heritability at mean environment ──────────────

# Vr (residual variance) is fixed externally (from an lme4 repeatability model), and not estimated here
Vr_chick <- 0.10^2    # 0.01
Vr_adult <- 0.619^2   # 0.383

VA_chick <- exp(post_draws_genetic$`B_v_A[1,1]`)
VA_adult <- exp(post_draws_genetic$`B_v_A[1,2]`)
VE_chick <- exp(post_draws_genetic$`B_v_E[1,1]`)
VE_adult <- exp(post_draws_genetic$`B_v_E[1,2]`)

VP_chick <- VA_chick + VE_chick + Vr_chick
VP_adult <- VA_adult + VE_adult + Vr_adult
h2_chick <- VA_chick / VP_chick
h2_adult <- VA_adult / VP_adult

report_var <- function(x, label) {
  cat(sprintf("%-10s median = %6.3f  [%6.3f, %6.3f]\n", label,
              median(x), quantile(x, 0.10), quantile(x, 0.90)))
}

cat("\n== Variance decomposition at mean environment ==\n")
cat("\nChick weight:\n")
report_var(VA_chick, "VA"); report_var(VE_chick, "VE")
cat(sprintf("%-10s = %.3f (fixed)\n", "Vr", Vr_chick))
report_var(VP_chick, "VP"); report_var(h2_chick, "h2")

cat("\nAdult weight:\n")
report_var(VA_adult, "VA"); report_var(VE_adult, "VE")
cat(sprintf("%-10s = %.3f (fixed)\n", "Vr", Vr_adult))
report_var(VP_adult, "VP"); report_var(h2_adult, "h2")

cat("\n== Baseline correlations at mean environment ==\n")
report_var(tanh(post_draws_genetic$`B_cpc_A[1,1]`), "r_A")
report_var(tanh(post_draws_genetic$`B_cpc_E[1,1]`), "r_E")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Plots to visualise CRN slopes 
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(ggdist)

#Load what is needed for making plots if running on it's own
post_crn_genetic <- readRDS(file.path(YOURDIRECTORY, "posterior_crn_genetic.RDS")) 
context_key      <- readRDS(file.path(ENTERYOURDIRECTORY, "context_key.RDS")) #from DATA PREP.R
 
# Predictor variables specifications 
genetic_specs <- list(
  list(id = "adult_temp",    raw_col = "adult_temp_ctx", vary_col = 3,
       xlab = "Adult temperature (\u00b0C)",  label = "adult temp",    env_group = "Temperature"),
  list(id = "natal_temp",    raw_col = "natal_temp_ctx", vary_col = 2,
       xlab = "Natal temperature (\u00b0C)",  label = "natal temp",    env_group = "Temperature"),
  list(id = "adult_density", raw_col = "adult_density",  vary_col = 5,
       xlab = "Adult population density",     label = "adult density", env_group = "Population density"),
  list(id = "natal_density", raw_col = "natal_density",  vary_col = 4,
       xlab = "Natal population density",     label = "natal density", env_group = "Population density")
)
 
# helper functions for plots
N_PRED <- 100
 
make_std_seq <- function(x_raw, n = N_PRED) seq(min(scale(x_raw)), max(scale(x_raw)), length.out = n)
to_real      <- function(x_std, mean_val, sd_val) x_std * sd_val + mean_val
 
make_Xpred <- function(vary_col, seq_z, P_y = 5) {
  X <- matrix(0, nrow = length(seq_z), ncol = P_y)
  X[, 1] <- 1
  X[, vary_col] <- seq_z
  X
}
 
get_B <- function(prefix, draws, P_y = 5) {
  sapply(seq_len(P_y), function(j) {
    col <- paste0(prefix, "[", j, ",1]")
    if (col %in% names(draws)) draws[[col]] else rep(NA, nrow(draws))
  })
}
 
predict_corr <- function(B_mat, Xpred, x_raw) {
  preds <- tanh(B_mat %*% t(Xpred))
  data.frame(
    x    = x_raw,
    med  = apply(preds, 2, median),
    lo95 = apply(preds, 2, quantile, 0.025), hi95 = apply(preds, 2, quantile, 0.975),
    lo80 = apply(preds, 2, quantile, 0.10),  hi80 = apply(preds, 2, quantile, 0.90),
    lo50 = apply(preds, 2, quantile, 0.25),  hi50 = apply(preds, 2, quantile, 0.75)
  )
}
 
col_A <- "#7F77DD"   # genetic (r_A)
col_E <- "#3DB88F"   # non-genetic (r_E)
 
theme_crn <- theme_classic(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(colour = "grey30"),
    axis.ticks       = element_line(colour = "grey30"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.text.y      = element_text(margin = margin(r = 6)),
    axis.text.x      = element_text(margin = margin(t = 6)),
    plot.margin      = margin(t = 15, r = 18, b = 15, l = 15),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

#function for crn plots 
make_overlay_panel <- function(df_A, df_E, xlab, rug_vals, ylim = c(-0.75, 1.0)) {
  df_A$comp <- "Genetic (r_A)"
  df_E$comp <- "Non-genetic (r_E)"
  df_both <- bind_rows(df_A, df_E)
  df_both$comp <- factor(df_both$comp, levels = c("Genetic (r_A)", "Non-genetic (r_E)"))
 
  ggplot(df_both, aes(x = x, fill = comp, colour = comp)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.10, colour = NA) +
    geom_ribbon(aes(ymin = lo80, ymax = hi80), alpha = 0.18, colour = NA) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.28, colour = NA) +
    geom_line(aes(y = med), linewidth = 1.2) +
    geom_rug(data = data.frame(x = rug_vals, comp = "Genetic (r_A)"),
             aes(x = x), inherit.aes = FALSE,
             colour = "grey50", alpha = 0.5, linewidth = 0.3, length = unit(0.03, "npc")) +
    scale_colour_manual(values = c("Genetic (r_A)" = col_A, "Non-genetic (r_E)" = col_E)) +
    scale_fill_manual(values   = c("Genetic (r_A)" = col_A, "Non-genetic (r_E)" = col_E)) +
    scale_y_continuous(limits = ylim, breaks = seq(-0.5, 1.0, 0.25),
                        expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = xlab, y = "Correlation") +
    theme_crn
}
 
# crn plot panels
B_cpc_A <- get_B("B_cpc_A", post_crn_genetic)
B_cpc_E <- get_B("B_cpc_E", post_crn_genetic)
 
overlay_panels <- list()
for (spec in genetic_specs) {
  raw_vals <- context_key[[spec$raw_col]]
  mean_val <- mean(raw_vals)
  sd_val   <- sd(raw_vals)
 
  z_seq   <- make_std_seq(raw_vals)
  raw_seq <- to_real(z_seq, mean_val, sd_val)
  Xp      <- make_Xpred(spec$vary_col, z_seq)
 
  df_A <- predict_corr(B_cpc_A, Xp, raw_seq)
  df_E <- predict_corr(B_cpc_E, Xp, raw_seq)
 
  overlay_panels[[spec$id]] <- make_overlay_panel(df_A, df_E, spec$xlab, raw_vals)
}
 
# posterior density plots
slopes_B <- bind_rows(lapply(genetic_specs, function(spec) {
  data.frame(
    value     = c(post_crn_genetic[[sprintf("B_cpc_A[%d,1]", spec$vary_col)]],
                  post_crn_genetic[[sprintf("B_cpc_E[%d,1]", spec$vary_col)]]),
    parameter = rep(c(paste0("r_A \u00d7 ", spec$label), paste0("r_E \u00d7 ", spec$label)),
                     each = nrow(post_crn_genetic)),
    component = rep(c("Genetic (r_A)", "Non-genetic (r_E)"), each = nrow(post_crn_genetic)),
    env_group = spec$env_group
  )
}))
 
param_levels <- unlist(lapply(genetic_specs, function(spec)
  c(paste0("r_A \u00d7 ", spec$label), paste0("r_E \u00d7 ", spec$label))))
 
slopes_B$parameter <- factor(slopes_B$parameter, levels = rev(param_levels))
slopes_B$component <- factor(slopes_B$component, levels = c("Genetic (r_A)", "Non-genetic (r_E)"))
slopes_B$env_group  <- factor(slopes_B$env_group, levels = c("Temperature", "Population density"))

#POSTERIOR DENSITY PLOT FOR FIGURE IN MAIN TEXT
p_halfeye_genetic <- ggplot(slopes_B, aes(x = value, y = parameter, fill = component, colour = component)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  stat_halfeye(
    .width = c(0.50, 0.80, 0.95), point_interval = "median_qi",
    slab_alpha = 0.60, normalize = "panels", point_size = 2.5,
    point_colour = "black"
  ) +
  scale_fill_manual(values   = c("Genetic (r_A)" = col_A, "Non-genetic (r_E)" = col_E)) +
  scale_colour_manual(values = c("Genetic (r_A)" = col_A, "Non-genetic (r_E)" = col_E)) +
  facet_wrap(~ env_group, scales = "free_y", ncol = 1) +
  labs(x = "Slope (atanh scale)", y = NULL) +
  theme_crn +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(size = 12, face = "bold", colour = "grey20"),
    panel.spacing    = unit(0.8, "lines")
  )
 
# ── Print each panel ─────────────────────────────────────────────────────────
print(overlay_panels$adult_temp)
print(overlay_panels$natal_temp)
print(overlay_panels$adult_density)
print(overlay_panels$natal_density)
print(p_halfeye_genetic)
 
