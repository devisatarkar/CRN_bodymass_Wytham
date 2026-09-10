functions {

//functions are used fom prior work by
//Dan Schrage (https://gitlab.com/dschrage/rcovreg)
  
  real sum_square_x(matrix x, int i, int j) {
    int j_prime;
    real sum_x = 0;
    if (j == 1) return(sum_x);
    j_prime = 1;
    while (j_prime < j) {
      sum_x += x[i, j_prime]^2;
      j_prime += 1;
    }
    return(sum_x);
  }

  matrix lkj_to_chol_corr(row_vector constrained_reals, int ntrait) {
    int z_counter;
    matrix[ntrait, ntrait] x;
    z_counter = 1;
    x[1, 1] = 1;
    for (j in 2:ntrait) {
      x[1, j] = 0;
    }
    for (i in 2:ntrait) {
      for (j in 1:ntrait) {
        if (i == j) {
          x[i, j] = sqrt(1 - sum_square_x(x, i, j));
        } else if (i > j) {
          x[i, j] = constrained_reals[z_counter] *
                    sqrt(1 - sum_square_x(x, i, j));
          z_counter += 1;
        } else {
          x[i, j] = 0;
        }
      }
    }
    return(x);
  }
}

data {
  int<lower=1> N;        // Number of chick-weight observations
  int<lower=1> M;        // Number of adult-weight observations
  int<lower=1> C;        // Number of paired contexts (birthyr × breedyr = 113)
  int<lower=1> C_birth;  // Number of unique birth years 
  int<lower=1> C_breed;  // Number of unique breeding years 
  int<lower=1> I;        // Number of unique individuals
  int<lower=1> D;        // Number of traits (D = 2: chick weight, adult weight)
  int<lower=1> P_y;      // Number of environmental predictors on the CRN (incl. intercept)
  int<lower=1> P_chick;  // Number of fixed-effect predictors for chick weight mean
  int<lower=1> P_adult;  // Number of fixed-effect predictors for adult weight mean

  int<lower=1> N_breedbox; // Number of unique nestboxes (adult breeding location RE)
  int<lower=1> N_broods;   // Number of unique broods (chick clutch RE)
// using two independent univariate REs

  // ── Index arrays linking observations to contexts ─────────────────────────
  array[N] int<lower=1> c_id_chick;   // Context (paired) index for each chick obs
  array[M] int<lower=1> c_id_adult;   // Context (paired) index for each adult obs

  // ── Index arrays linking observations to positions in cmat ───────────────
  array[N] int<lower=1> idc_chick;    // Position in mat_G column 1 for chick obs
  array[M] int<lower=1> idc_adult;    // Position in mat_G column 2 for adult obs

  // ── Index arrays linking observations to fixed-effect linear predictors ──
  array[N] int<lower=1> id_chick_lm;  // Row index into mu_chick_fe for each chick obs
  array[M] int<lower=1> id_adult_lm;  // Row index into mu_adult_fe for each adult obs

  // ── Index arrays linking observations to individuals ─────────────────────
  array[N] int<lower=1> id_chick;     // Individual ID for each chick obs
  array[M] int<lower=1> id_adult;     // Individual ID for each adult obs

  // ── Season random effect indices ──────────────────────────────────────────
  array[N] int<lower=1> season_id_chick; // Birth year index (1:C_birth) for each chick obs
  array[M] int<lower=1> season_id_adult; // Breed year index (1:C_breed) for each adult obs

  // ── Grouping indices for nestbox and brood REs ────────────────────────────
  array[M] int<lower=1> breedbox_id;  // Nestbox ID for each adult obs
  array[N] int<lower=1> brood_id_chick; // Brood ID for each chick obs
                                      // seperate univariate REs for nestbox (adult) and brood (chick)

  // ── Environmental predictor matrices ─────────────────────────────────────
  matrix[C, P_y]     X;       // CRN predictors (intercept, natal temp, adult temp)
  matrix[N, P_chick] X_chick; // Fixed-effect predictors for chick weight mean
  matrix[M, P_adult] X_adult; // Fixed-effect predictors for adult weight mean

  // ── Context structure for mat_G construction ─────────────────────────────
  int<lower=1>       cm;           // Max individuals observed in any single context
  array[C, cm] int   cmat;         // Matrix of individual IDs per context (C rows, cm cols)
  array[C]     int   cn;           // Count of individuals in each context
  int<lower=1>       cnt;          // Total individual-context observations (sum of cn)

  // ── Response variables ────────────────────────────────────────────────────
  vector[N] chick_weight; // Day-15 chick mass (standardised)
  vector[M] adult_weight; // Adult body mass (standardised)
}

transformed data {
   // No LA = cholesky_decompose(A) here. removed identity matrix
  int ncor = (D * (D - 1)) %/% 2; // Number of unique correlation parameters
                                    // For D=2 traits: ncor = 1 (one correlation)
  real sd_E_chick = 0.10; // sd_E_chick = 0.10:  as measurement error only (one obs per individual)
  real sd_E_adult = 0.619; // sd_E_adult = 0.619: from lme4 REML on the paired-context subset; standardised = 0.619

  matrix[C, P_y]     Q      = qr_thin_Q(X) * sqrt(C - 1);
  matrix[P_y, P_y]   R      = qr_thin_R(X) / sqrt(C - 1);
  matrix[P_y, P_y]   R_inv  = inverse(R);

  matrix[N, P_chick]       Q_chick     = qr_thin_Q(X_chick) * sqrt(N - 1);
  matrix[P_chick, P_chick] R_chick     = qr_thin_R(X_chick) / sqrt(N - 1);
  matrix[P_chick, P_chick] R_inv_chick = inverse(R_chick);

  matrix[M, P_adult]       Q_adult     = qr_thin_Q(X_adult) * sqrt(M - 1);
  matrix[P_adult, P_adult] R_adult     = qr_thin_R(X_adult) / sqrt(M - 1);
  matrix[P_adult, P_adult] R_inv_adult = inverse(R_adult);
}

parameters {
  // ── Fixed effects  ───────
  vector[P_chick] B_mq_chick; // Chick weight mean reaction norm 
  vector[P_adult] B_mq_adult; // Adult weight mean reaction norm 

  // ── CRN parameters ─────────────────────
  matrix[P_y, ncor] B_cpcq; // Slopes of canonical partial correlations on environment (QR)
                              // After tanh back-transform: how r(chick,adult) varies with
                              // natal and adult env
  matrix[P_y, D]    B_vq;   // Slopes of log among-individual variances on environment (QR)
                              // After exp back-transform: how SD(chick) and SD(adult) vary

  // ── Individual-level random effects 
  matrix[cnt, D] Z_G; 
  // ── Season random effects (non-centred, split by life stage) ─────────────
  vector[C_birth] z_season_chick; // birth-year season effects for chick weight
  real<lower=0>   sd_season_chick; // SD of birth-year season effects
  vector[C_breed] z_season_adult;  // breed-year season effects for adult weight
  real<lower=0>   sd_season_adult; // SD of breed-year season effects
 //natal vs breeding environments separate

  // ── Nestbox and brood random effects ──────
  vector[N_breedbox] z_breedbox;    // nestbox effects on adult weight
  real<lower=0>      sd_breedbox;   // SD of nestbox effects
  vector[N_broods] z_brood_chick;   //brood effects on chick weight
  real<lower=0>    sd_brood_chick;  // SD of brood effects
}

model {
  // ── Fixed-effect linear predictors ────────────────────
  vector[N] mu_chick_fe = Q_chick * B_mq_chick; 
  vector[M] mu_adult_fe = Q_adult * B_mq_adult; 

  // ── CRN: context-specific correlations and SDs ───────────────────────────
  matrix[C, ncor] cpc_G = tanh(Q * B_cpcq); 
  matrix[C, D]    sd_G  = sqrt(exp(Q * B_vq)); 

  // ── Scale season random effects ───────────────────────────────────────────
  vector[C_birth]  re_season_chick = z_season_chick * sd_season_chick;
  vector[C_breed]  re_season_adult = z_season_adult * sd_season_adult;

  // ── Scale nestbox and brood random effects ────────────────────────────────
  vector[N_breedbox] re_breedbox    = z_breedbox     * sd_breedbox;
  vector[N_broods]   re_brood_chick = z_brood_chick * sd_brood_chick;


  vector[N] mu_chick = mu_chick_fe[id_chick_lm]
                       + re_season_chick[season_id_chick]  // birth-year season
                       + re_brood_chick[brood_id_chick];   // brood/clutch effect

  vector[M] mu_adult = mu_adult_fe[id_adult_lm]
                       + re_season_adult[season_id_adult]  // breed-year season
                       + re_breedbox[breedbox_id];         // nestbox effect

  matrix[cnt, D] mat_G;
  int pos = 1;
  for (c in 1:C) {
    mat_G[pos:(pos + cn[c] - 1)] =
      Z_G[pos:(pos + cn[c] - 1)]
      * diag_pre_multiply(sd_G[c],
          lkj_to_chol_corr(cpc_G[c], D))';
    pos = pos + cn[c];
  }

  // ── Add individual effects to linear predictors ───────────────────────────
  for (n in 1:N) {
    mu_chick[n] += col(mat_G, 1)[idc_chick[n]]; // Column 1 = chick weight individual effect
  }
  for (m in 1:M) {
    mu_adult[m] += col(mat_G, 2)[idc_adult[m]]; // Column 2 = adult weight individual effect
  }

  // ── Likelihoods ───────────────────────────────────────────────────────────
  // Both traits Gaussian with externally fixed residual SDs.
 
  chick_weight ~ normal(mu_chick, sd_E_chick);
  adult_weight ~ normal(mu_adult, sd_E_adult);

  // ── Priors ────────────────────────────────────────────────────────────────
  to_vector(B_mq_chick) ~ normal(0, 1); 
  to_vector(B_mq_adult) ~ normal(0, 1);


  B_cpcq[1, ]          ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_cpcq[2:P_y, ]) ~ normal(0, 0.5);
  B_vq[1, ]            ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_vq[2:P_y, ])   ~ normal(0, 1); 

  to_vector(Z_G) ~ std_normal(); 

  // Season REs
  z_season_chick  ~ std_normal();
  z_season_adult  ~ std_normal();
  sd_season_chick ~ exponential(2); 
  sd_season_adult ~ exponential(2);

  // Nestbox and brood REs
  z_breedbox     ~ std_normal();
  sd_breedbox    ~ exponential(2);
  z_brood_chick  ~ std_normal();
  sd_brood_chick ~ exponential(2);

}

generated quantities {
  // ── Back-transform fixed effects from QR to original scale ───────────────
  vector[P_chick] B_m_chick = R_inv_chick * B_mq_chick;
  vector[P_adult] B_m_adult = R_inv_adult * B_mq_adult;

  // ── Back-transform CRN parameters from QR to original scale ─────────────

  matrix[P_y, ncor] B_cpc;
  matrix[P_y, D]    B_v;

  for (d in 1:ncor) {
    B_cpc[, d] = R_inv * B_cpcq[, d];
  }
  for (d in 1:D) {
    B_v[, d] = R_inv * B_vq[, d];
  }

  // ── Posterior predictive check quantities ─────────────────────────────────

  matrix[C, ncor] cpc_G_bis = tanh(Q * B_cpcq);
  matrix[C, D]    sd_G_bis  = sqrt(exp(Q * B_vq));

  vector[C_birth]  re_season_chick_bis = z_season_chick * sd_season_chick;
  vector[C_breed]  re_season_adult_bis = z_season_adult * sd_season_adult;

  vector[N_breedbox] re_breedbox_bis    = z_breedbox     * sd_breedbox;
  vector[N_broods]   re_brood_chick_bis = z_brood_chick * sd_brood_chick;

  vector[N] mu_chick_bis = (Q_chick * B_mq_chick)[id_chick_lm]
                           + re_season_chick_bis[season_id_chick]
                           + re_brood_chick_bis[brood_id_chick];
  vector[M] mu_adult_bis = (Q_adult * B_mq_adult)[id_adult_lm]
                           + re_season_adult_bis[season_id_adult]
                           + re_breedbox_bis[breedbox_id];

  matrix[cnt, D] mat_G_bis;
  {
    int pos_bis = 1;
    for (c in 1:C) {
      mat_G_bis[pos_bis:(pos_bis + cn[c] - 1)] =
        Z_G[pos_bis:(pos_bis + cn[c] - 1)]
        * diag_pre_multiply(sd_G_bis[c],
            lkj_to_chol_corr(cpc_G_bis[c], D))';
      pos_bis = pos_bis + cn[c];
    }
  }

  // ── Draw posterior predictive replicates ──────────────────────────────────
  array[N] real y_rep_chick;
  array[M] real y_rep_adult;

  for (n in 1:N) {
    mu_chick_bis[n] += col(mat_G_bis, 1)[idc_chick[n]];
    y_rep_chick[n]   = normal_rng(mu_chick_bis[n], sd_E_chick);
  }
  for (m in 1:M) {
    mu_adult_bis[m] += col(mat_G_bis, 2)[idc_adult[m]];
    y_rep_adult[m]   = normal_rng(mu_adult_bis[m], sd_E_adult);
  }

}