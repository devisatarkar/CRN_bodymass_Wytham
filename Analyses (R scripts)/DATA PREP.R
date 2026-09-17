# ==============================================================================
# PREPARING DATA, MATRICES AND STAN LIST FOR MODEL
# ==============================================================================

library(dplyr)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load dataframes
# ══════════════════════════════════════════════════════════════════════════════

df_chick_p <- read.csv("ENTERYOURDIRECTORY/df_chick_p.csv")
df_adult_p <- read.csv("ENTERYOURDIRECTORY/df_adult_p.csv")

gtbt_density_peryear <- read.csv("ENTERYOURDIRECTORY/gtbt_density_peryear.csv")

I_p <- n_distinct(df_chick_p$id_num) #number of individuals
C_p <- n_distinct(df_adult_p$context_id) #number of contexts

cat("── Loaded data ──\n")
cat("df_chick_p rows (individuals):      ", nrow(df_chick_p), "\n")
cat("df_adult_p rows (breeding attempts):", nrow(df_adult_p), "\n")
cat("I (individuals):", I_p, "  C (contexts):", C_p, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Create context_key from df_chick_p + df_adult_p
# Since context-level temperature and rain must be averaged over every individual x breeding attempt row
# join the individual-level natal covariates from df_chick_p onto df_adult_p by id_num, then aggregate by context_id
# ══════════════════════════════════════════════════════════════════════════════
adult_with_natal <- df_adult_p %>%
  left_join(
    df_chick_p %>% select(id_num, birthyear, natal_meantemp, natal_meanrain),
    by = "id_num"
  )

context_key <- adult_with_natal %>%
  group_by(context_id) %>%
  summarise(
    n_ind          = n_distinct(id_num),
    birthyear      = first(birthyear),
    breedyear      = first(breedyear),
    natal_temp_ctx = mean(natal_meantemp, na.rm = TRUE),
    adult_temp_ctx = mean(breed_meantemp, na.rm = TRUE),
    natal_rain_ctx = mean(natal_meanrain, na.rm = TRUE),
    adult_rain_ctx = mean(breed_meanrain, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  arrange(context_id)

cat("Contexts:", nrow(context_key), "\n")
for (v in c("natal_temp_ctx", "adult_temp_ctx", "natal_rain_ctx", "adult_rain_ctx")) {
  cat(sprintf("  %-16s NAs = %d\n", v, sum(is.na(context_key[[v]]))))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Add population density
# ══════════════════════════════════════════════════════════════════════════════
combined_density <- gtbt_density_peryear %>%
  rename(combined_density = popdens_bt_gt)

context_years <- union(context_key$birthyear, context_key$breedyear) %>% unique() %>% sort()
missing_years <- setdiff(context_years, combined_density$year)
cat("\nContext years missing from density data:",
    if (length(missing_years) == 0) "none" else paste(missing_years, collapse = ", "), "\n")

context_key <- context_key %>%
  left_join(combined_density %>% select(year, combined_density) %>%
              rename(natal_density = combined_density),
            by = c("birthyear" = "year")) %>%
  left_join(combined_density %>% select(year, combined_density) %>%
              rename(adult_density = combined_density),
            by = c("breedyear" = "year"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Final context-level matrix Xc: here X4_part1, P_y = 9
# Temperature has quadratic terms (natal and adult); rainfall and density are linear.
# ══════════════════════════════════════════════════════════════════════════════
natal_temp_std <- scale(context_key$natal_temp_ctx)[, 1]
adult_temp_std <- scale(context_key$adult_temp_ctx)[, 1]
natal_rain_std <- scale(context_key$natal_rain_ctx)[, 1]
adult_rain_std <- scale(context_key$adult_rain_ctx)[, 1]
natal_dens_std <- scale(context_key$natal_density)[, 1]
adult_dens_std <- scale(context_key$adult_density)[, 1]

X_check <- cbind(
  natal_temp  = natal_temp_std, natal_temp2 = natal_temp_std^2,
  adult_temp  = adult_temp_std, adult_temp2 = adult_temp_std^2,
  natal_rain  = natal_rain_std, adult_rain  = adult_rain_std,
  natal_dens  = natal_dens_std, adult_dens  = adult_dens_std
)

#checking VIF and collinearity between variables
cat("\n── Correlation matrix, Part 1 predictor set ──\n")
print(round(cor(X_check), 3))

cat("\n── VIF, Part 1 predictor set ──\n")
for (j in seq_len(ncol(X_check))) {
  r2 <- summary(lm(X_check[, j] ~ X_check[, -j]))$r.squared
  cat(sprintf("  %-12s VIF = %.2f\n", colnames(X_check)[j], 1 / (1 - r2)))
}

X4_part1 <- as.matrix(cbind(
  intercept     = 1,
  natal_temp    = natal_temp_std,
  natal_temp2   = natal_temp_std^2,
  adult_temp    = adult_temp_std,
  adult_temp2   = adult_temp_std^2,
  natal_rain    = natal_rain_std,
  adult_rain    = adult_rain_std,
  natal_density = natal_dens_std,
  adult_density = adult_dens_std
))
rownames(X4_part1) <- context_key$context_id

cat("\n── X4_part1 ──\n")
cat("Rows:", nrow(X4_part1), "  Cols:", ncol(X4_part1), "  NAs:", sum(is.na(X4_part1)), "\n")

# Parameter coding for B_cpc / B_v (columns of X4_part1):
#   [2] natal_temp    [3] natal_temp2
#   [4] adult_temp    [5] adult_temp2
#   [6] natal_rain    [7] adult_rain
#   [8] natal_density [9] adult_density

#for back-transforming later
scaling_part1 <- list(
  natal_temp_mean = mean(context_key$natal_temp_ctx), natal_temp_sd = sd(context_key$natal_temp_ctx),
  adult_temp_mean = mean(context_key$adult_temp_ctx), adult_temp_sd = sd(context_key$adult_temp_ctx),
  natal_rain_mean = mean(context_key$natal_rain_ctx), natal_rain_sd = sd(context_key$natal_rain_ctx),
  adult_rain_mean = mean(context_key$adult_rain_ctx), adult_rain_sd = sd(context_key$adult_rain_ctx),
  natal_dens_mean = mean(context_key$natal_density),  natal_dens_sd = sd(context_key$natal_density),
  adult_dens_mean = mean(context_key$adult_density),  adult_dens_sd = sd(context_key$adult_density)
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Indexing for model in stan (cmat, corder, idc, c_id)
# ══════════════════════════════════════════════════════════════════════════════

# For each paired context, which unique individuals appear in df_adult_p?
c_list_p <- vector("list", C_p)
for (c in 1:C_p) {
  ctx <- context_key$context_id[c]
  c_list_p[[c]] <- sort(unique(df_adult_p$id_num[df_adult_p$context_id == ctx]))
}
cn_p <- sapply(c_list_p, length)

cat("─Individuals per paired context ──\n")
cat("Min:         ", min(cn_p), "\n") #10
cat("Median:      ", median(cn_p), "\n") #35
cat("Max:         ", max(cn_p), "\n") #234

cm_p <- max(cn_p)
cat("cm (max individuals in any context):", cm_p, "\n\n") #234

# Build cmat_p (C_p x cm_p)
cmat_p <- matrix(0L, nrow = C_p, ncol = cm_p)
for (c in 1:C_p) {
  ids <- c_list_p[[c]]
  cmat_p[c, 1:length(ids)] <- ids
}
stopifnot(sum(cmat_p > 0) == sum(cn_p))

temp_p <- t(cmat_p)
corder_p <- data.frame(
  id_num     = temp_p[temp_p > 0],
  context_id = rep(seq_len(C_p), times = cn_p)
)
corder_p$pos <- seq_len(nrow(corder_p))

# create idc_adult_p 
# For each row in df_adult_p, what is the position of that individual
# in that paired context's slot in cmat_p?

idc_adult_p <- corder_p$pos[
  match(paste(df_adult_p$id_num, df_adult_p$context_id, sep = "."),
        paste(corder_p$id_num, corder_p$context_id, sep = "."))
]
stopifnot(length(idc_adult_p) == nrow(df_adult_p))
stopifnot(max(idc_adult_p) == sum(cn_p))
stopifnot(sum(is.na(idc_adult_p)) == 0)

# create idc_chick_p 
first_appearance_p <- corder_p %>% group_by(id_num) %>% slice(1) %>% ungroup()
idc_chick_p <- first_appearance_p$pos[match(df_chick_p$id_num, first_appearance_p$id_num)]
stopifnot(length(idc_chick_p) == nrow(df_chick_p))
stopifnot(sum(is.na(idc_chick_p)) == 0)

# ──  Build c_id vectors ─────────
c_id_adult_p <- df_adult_p$context_id #paired context index for each adult observation

first_context_p <- corder_p %>% group_by(id_num) %>% slice(1) %>% ungroup() %>%
  select(id_num, context_id)
c_id_chick_p <- first_context_p$context_id[match(df_chick_p$id_num, first_context_p$id_num)] #paired context index for each chick observation

stopifnot(sum(is.na(c_id_adult_p)) == 0, sum(is.na(c_id_chick_p)) == 0)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Season (birth/breed year), nestbox, and brood random-effect indices
# ══════════════════════════════════════════════════════════════════════════════

#birthyear random effect
key_birthyr_p <- sort(unique(df_chick_p$birthyear))
season_id_chick_p <- match(df_chick_p$birthyear, key_birthyr_p)
n_birthyears_p <- length(key_birthyr_p)

#breedyear random effect
key_breedyr_p <- sort(unique(df_adult_p$breedyear))
season_id_adult_p <- match(df_adult_p$breedyear, key_breedyr_p)
n_breedyears_p <- length(key_breedyr_p)

stopifnot(min(season_id_chick_p) == 1, max(season_id_chick_p) == n_birthyears_p)
stopifnot(min(season_id_adult_p) == 1, max(season_id_adult_p) == n_breedyears_p)

#brood id and breeding nestbox random effects
breedbox_id_p    <- df_adult_p$breedbox_num
brood_id_chick_p <- df_chick_p$brood_num

stopifnot(min(breedbox_id_p) == 1, max(breedbox_id_p) == n_distinct(breedbox_id_p))
stopifnot(min(brood_id_chick_p) == 1, max(brood_id_chick_p) == n_distinct(brood_id_chick_p))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Fixed-effects matrices X1: here X_chick_p, and X2: here X_adult_p
# ══════════════════════════════════════════════════════════════════════════════
X_chick_p <- as.matrix(cbind(
  intercept            = 1,
  april_birthlaydate   = scale(df_chick_p$april_birthlaydate)[, 1],
  natalbroodsize       = scale(df_chick_p$natalbroodsize)[, 1],
  natal_meantemp       = scale(df_chick_p$natal_meantemp)[, 1],
  natal_meanrain       = scale(df_chick_p$natal_meanrain)[, 1],
  absnatal_hd_mismatch = scale(df_chick_p$absnatal_hd_mismatch)[, 1],
  birthyear            = scale(df_chick_p$birthyear)[, 1]
))

X_adult_p <- as.matrix(cbind(
  intercept      = 1,
  april_laydate  = scale(df_adult_p$april_laydate)[, 1],
  breed_meantemp = scale(df_adult_p$breed_meantemp)[, 1],
  breed_meanrain = scale(df_adult_p$breed_meanrain)[, 1],
  num_chicks     = scale(df_adult_p$num_chicks)[, 1],
  age_years      = scale(df_adult_p$age_years)[, 1],
  sex_num        = df_adult_p$sex_num,
  breedyear      = scale(df_adult_p$breedyear)[, 1]
))

stopifnot(nrow(X_chick_p) == nrow(df_chick_p), sum(is.na(X_chick_p)) == 0)
stopifnot(nrow(X_adult_p) == nrow(df_adult_p), sum(is.na(X_adult_p)) == 0)

cat("\n── Fixed-effects matrices ──\n")
cat("X_chick_p:", nrow(X_chick_p), "x", ncol(X_chick_p),
    " X_adult_p:", nrow(X_adult_p), "x", ncol(X_adult_p), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Scale response variables
# ══════════════════════════════════════════════════════════════════════════════
chick_mean_p <- mean(df_chick_p$chickweight)
chick_sd_p   <- sd(df_chick_p$chickweight)
adult_mean_p <- mean(df_adult_p$adultweight)
adult_sd_p   <- sd(df_adult_p$adultweight)

chick_weight_scaled_p <- (df_chick_p$chickweight - chick_mean_p) / chick_sd_p
adult_weight_scaled_p <- (df_adult_p$adultweight - adult_mean_p) / adult_sd_p

scaling_params_p <- list(
  chick_mean = chick_mean_p, chick_sd = chick_sd_p,
  adult_mean = adult_mean_p, adult_sd = adult_sd_p
)

cat("Chick weight — mean:", round(chick_mean_p, 3),
    "sd:", round(chick_sd_p, 3), "\n")
cat("Adult weight — mean:", round(adult_mean_p, 3),
    "sd:", round(adult_sd_p, 3), "\n")
cat("Any NAs in chick_weight_scaled_p:", sum(is.na(chick_weight_scaled_p)), "\n")
cat("Any NAs in adult_weight_scaled_p:", sum(is.na(adult_weight_scaled_p)), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Assemble the final, complete Stan data list
# ══════════════════════════════════════════════════════════════════════════════
stan_dl_part1 <- list(
  N           = nrow(df_chick_p), #number of individuals
  M           = nrow(df_adult_p), #number of breeding attempts
  C           = C_p,              #number of contexts
  C_birth     = n_birthyears_p,   #number of birthyears 
  C_breed     = n_breedyears_p,   #number of breedyears
  I           = I_p,              #number of individuals
  D           = 2L,               # number of response variables
  P_y         = ncol(X4_part1),   # number of context variables (incl intercept)
  P_chick     = ncol(X_chick_p),  # number of fixed effects for nestling mass
  P_adult     = ncol(X_adult_p),  # number of fixed effects for adult mass
  N_breedbox  = n_distinct(df_adult_p$breedbox_num), #number of breeding nest boxes
  N_broods    = n_distinct(df_chick_p$brood_num), #number of unique broods
  
  # indexing over paired contexts
  
  c_id_chick  = c_id_chick_p, 
  c_id_adult  = c_id_adult_p,
  idc_chick   = idc_chick_p,
  idc_adult   = idc_adult_p,
  id_chick_lm = seq_len(nrow(df_chick_p)),
  id_adult_lm = seq_len(nrow(df_adult_p)),
  
  #individual indices
  id_chick    = df_chick_p$id_num,
  id_adult    = df_adult_p$id_num,
  
  #year random effects
  season_id_chick = season_id_chick_p,
  season_id_adult = season_id_adult_p,
  
  #nestbox and brood identity RE
  breedbox_id    = breedbox_id_p,
  brood_id_chick = brood_id_chick_p,
  
  #predictor matrices
  X       = X4_part1,
  X_chick = X_chick_p,
  X_adult = X_adult_p,
  
  #cmat structure
  cmat = cmat_p,
  cm   = cm_p,
  cn   = cn_p,
  cnt  = sum(cn_p),
  
  #response variables
  chick_weight = chick_weight_scaled_p,
  adult_weight = adult_weight_scaled_p
)

stopifnot(nrow(stan_dl_part1$X) == stan_dl_part1$C, ncol(stan_dl_part1$X) == stan_dl_part1$P_y)
stopifnot(nrow(stan_dl_part1$X_chick) == stan_dl_part1$N)
stopifnot(nrow(stan_dl_part1$X_adult) == stan_dl_part1$M)

cat("\n── stan_dl_part1 NA check ──\n")
any_na <- FALSE
for (nm in names(stan_dl_part1)) {
  n_na <- sum(is.na(stan_dl_part1[[nm]]))
  if (n_na > 0) { cat("  WARNING —", nm, "has", n_na, "NAs\n"); any_na <- TRUE }
}
if (!any_na) cat("  All elements NA-free\n")

cat("\n── Final dimensions ──\n")
cat("N:", stan_dl_part1$N, " M:", stan_dl_part1$M, " C:", stan_dl_part1$C,
    " I:", stan_dl_part1$I, " P_y:", stan_dl_part1$P_y,
    " P_chick:", stan_dl_part1$P_chick, " P_adult:", stan_dl_part1$P_adult, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 — Save
# ══════════════════════════════════════════════════════════════════════════════

saveRDS(stan_dl_part1, "stan_dl_part1.RDS")
saveRDS(X4_part1, "X4_part1.RDS")
saveRDS(scaling_part1, "scaling_part1.RDS")
saveRDS(context_key, "context_key.RDS")
saveRDS(scaling_params_p, "scaling_params_p.RDS")
