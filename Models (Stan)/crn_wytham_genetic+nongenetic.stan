

functions {
  // ── Unchanged from phenotypic model  ────────────────────
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
    for (j in 2:ntrait) x[1, j] = 0;
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
  
  int<lower=1> N;          // chick observations = 4345
  int<lower=1> M;          // adult observations = 6362
  int<lower=1> C;          // paired contexts = 113
  int<lower=1> C_birth;    // birth year levels
  int<lower=1> C_breed;    // breed year levels
  int<lower=1> I;          // analysis individuals = 4345
  int<lower=1> D;          // traits = 2
  int<lower=1> P_y;        // context matrix cols = 3
  int<lower=1> P_chick;    // chick fixed effect cols = 7
  int<lower=1> P_adult;    // adult fixed effect cols = 8
  int<lower=1> N_breedbox;
  int<lower=1> N_broods;

  //  indexing ─────────────────────────────────────────────────────
  array[N] int<lower=1> c_id_chick;
  array[M] int<lower=1> c_id_adult;
  array[N] int<lower=1> idc_chick;
  array[M] int<lower=1> idc_adult;
  array[N] int<lower=1> id_chick_lm;
  array[M] int<lower=1> id_adult_lm;
  array[N] int<lower=1> id_chick;
  array[M] int<lower=1> id_adult;
  array[N] int<lower=1> season_id_chick;
  array[M] int<lower=1> season_id_adult;
  array[M] int<lower=1> breedbox_id;
  array[N] int<lower=1> brood_id_chick;

  // ── predictor matrices ──────────────────────────────────────────────
  matrix[C, P_y]     X;
  matrix[N, P_chick] X_chick;
  matrix[M, P_adult] X_adult;

  //  ────────────────────────────────────────────
  int<lower=1>       cm;
  array[C, cm] int   cmat;
  array[C]     int   cn;
  int<lower=1>       cnt;

  // ──  response variables ────────────────────────────────────────────────────
  vector[N] chick_weight;
  vector[M] adult_weight;

  // ── Relatedness matrix ────────────────────────────────────────────────
  matrix[I, I] A;
}

transformed data {
  int ncor = (D * (D - 1)) %/% 2;  // = 1 for D=2

  // ── Fixed residual SDs — unchanged ────────────────────────────────────────
  real sd_E_chick = 0.10;
  real sd_E_adult = 0.619;

  // ── new: Cholesky factor of A ──────────────────────────────────────────────

  matrix[I, I] LA = cholesky_decompose(A);

  // ── QR reparameterisation  ─────────────────────────────────────
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
  // ── Fixed effects — unchanged ──────────────────────────────────────────────
  vector[P_chick] B_mq_chick;
  vector[P_adult] B_mq_adult;

  // ── Genetic CRN parameters ─────────────────────────────────────────────────
  // These replace the old single B_cpcq and B_vq
  // B_cpcq_A: reaction norm slopes for genetic correlation (QR scale)
  // B_vq_A:   reaction norm slopes for genetic variance (QR scale)
  // Back-transformed to B_cpc_A and B_v_A in generated quantities
  matrix[P_y, ncor] B_cpcq_A;
  matrix[P_y, D]    B_vq_A;

  // ── Non-genetic individual CRN parameters ──────────────────────────────────
  // New parameters — capture environmental carry-over beyond genetic effects
  // B_cpcq_E: reaction norm slopes for non-genetic correlation (QR scale)
  // B_vq_E:   reaction norm slopes for non-genetic individual variance (QR scale)
  matrix[P_y, ncor] B_cpcq_E;
  matrix[P_y, D]    B_vq_E;

  // ── Genetic random effects ──────────────────────
  // Z_G: cnt x D matrix of standardised genetic deviations
  // Same dimensions and indexing as old Z_G in phenotypic model
  // Prior is std_normal — pedigree structure enters through LA multiplication
  // in the model block, not through the prior
  matrix[cnt, D] Z_G;

  // ── Non-genetic individual random effects ──────────────────────────────────
  // Z_E: cnt x D matrix — same structure as Z_G but NO pedigree constraint
  // Independent across individuals — captures stable non-genetic differences
  // Renamed from old Z_G to distinguish from genetic component
  matrix[cnt, D] Z_E;

  // ── Season, breedbox, brood random effects — unchanged ────────────────────
  vector[C_birth] z_season_chick;
  real<lower=0>   sd_season_chick;
  vector[C_breed] z_season_adult;
  real<lower=0>   sd_season_adult;
  vector[N_breedbox] z_breedbox;
  real<lower=0>      sd_breedbox;
  vector[N_broods] z_brood_chick;
  real<lower=0>    sd_brood_chick;
}

model {
  // ── Fixed effect means — unchanged ────────────────────────────────────────
  vector[N] mu_chick_fe = Q_chick * B_mq_chick;
  vector[M] mu_adult_fe = Q_adult * B_mq_adult;

  // ── Context-specific genetic CRN predictions ───────────────────────────────
  // sd_G[c]: genetic SDs for each trait in context c
  // cpc_G[c]: genetic correlation between traits in context c
  matrix[C, D]    sd_G  = sqrt(exp(Q * B_vq_A));
  matrix[C, ncor] cpc_G = tanh(Q * B_cpcq_A);

  // ── Context-specific non-genetic individual CRN predictions ───────────────
  // sd_Ei[c]: non-genetic individual SDs for each trait in context c
  // cpc_Ei[c]: non-genetic individual correlation between traits in context c
  matrix[C, D]    sd_Ei  = sqrt(exp(Q * B_vq_E));
  matrix[C, ncor] cpc_Ei = tanh(Q * B_cpcq_E);

  // ── Season/breedbox/brood random effects — unchanged ──────────────────────
  vector[C_birth]    re_season_chick = z_season_chick * sd_season_chick;
  vector[C_breed]    re_season_adult = z_season_adult * sd_season_adult;
  vector[N_breedbox] re_breedbox     = z_breedbox     * sd_breedbox;
  vector[N_broods]   re_brood_chick  = z_brood_chick  * sd_brood_chick;

  // ── Linear predictors — unchanged ─────────────────────────────────────────
  vector[N] mu_chick = mu_chick_fe[id_chick_lm]
                       + re_season_chick[season_id_chick]
                       + re_brood_chick[brood_id_chick];
  vector[M] mu_adult = mu_adult_fe[id_adult_lm]
                       + re_season_adult[season_id_adult]
                       + re_breedbox[breedbox_id];

  // ── Genetic random effects 
  
  matrix[cnt, D] mat_G;
  {
    int pos = 1;
    for (c in 1:C) {
      mat_G[pos:(pos + cn[c] - 1)] =
        LA[cmat[c, 1:cn[c]], cmat[c, 1:cn[c]]]   // pedigree 
        * Z_G[pos:(pos + cn[c] - 1)]             
        * diag_pre_multiply(sd_G[c],              
            lkj_to_chol_corr(cpc_G[c], D))';     
      pos = pos + cn[c];
    }
  }

  // ── Non-genetic individual effects — same loop structure, no LA ─────────────
  // Z_E is NOT multiplied by LA — these effects are independent across
  // individuals, so there is no pedigree constraint.
  // Otherwise identical structure to mat_G loop.
  matrix[cnt, D] mat_E;
  {
    int pos = 1;
    for (c in 1:C) {
      mat_E[pos:(pos + cn[c] - 1)] =
        Z_E[pos:(pos + cn[c] - 1)]                
        * diag_pre_multiply(sd_Ei[c],            
            lkj_to_chol_corr(cpc_Ei[c], D))';    
      pos = pos + cn[c];
    }
  }

  // ── genetic + non-genetic effects to linear predictors ────────────────
  
  for (n in 1:N) {
    mu_chick[n] += col(mat_G, 1)[idc_chick[n]]   // genetic chick effect
                 + col(mat_E, 1)[idc_chick[n]];   // non-genetic chick effect
  }
  for (m in 1:M) {
    mu_adult[m] += col(mat_G, 2)[idc_adult[m]]   // genetic adult effect
                 + col(mat_E, 2)[idc_adult[m]];   // non-genetic adult effect
  }

  // ── Likelihood — unchanged ─────────────────────────────────────────────────
  chick_weight ~ normal(mu_chick, sd_E_chick);
  adult_weight ~ normal(mu_adult, sd_E_adult);

  // ── Priors — genetic CRN parameters ───────────────────────────────────────
  // Same prior structure as old B_cpcq and B_vq in phenotypic model
  B_cpcq_A[1, ]                  ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_cpcq_A[2:P_y, ]) ~ normal(0, 0.5);
  B_vq_A[1, ]                    ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_vq_A[2:P_y, ])   ~ normal(0, 1);

  // ── Priors — non-genetic individual CRN parameters ────────────────────────
  // same prior structure as genetic component
  B_cpcq_E[1, ]                  ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_cpcq_E[2:P_y, ]) ~ normal(0, 0.5);
  B_vq_E[1, ]                    ~ normal(0, 0.5);
  if (P_y > 1)
    to_vector(B_vq_E[2:P_y, ])   ~ normal(0, 1);

  // ── Priors — random effects ────────────────────────────────────────────────

  to_vector(Z_G) ~ std_normal();

  // Z_E: std_normal prior — same as old Z_G prior in phenotypic model
  to_vector(Z_E) ~ std_normal();

  // ── Priors — fixed effects — unchanged ────────────────────────────────────
  to_vector(B_mq_chick) ~ normal(0, 1);
  to_vector(B_mq_adult) ~ normal(0, 1);

  // ── Priors — season/breedbox/brood — unchanged ────────────────────────────
  z_season_chick  ~ std_normal();
  z_season_adult  ~ std_normal();
  sd_season_chick ~ exponential(2);
  sd_season_adult ~ exponential(2);
  z_breedbox      ~ std_normal();
  sd_breedbox     ~ exponential(2);
  z_brood_chick   ~ std_normal();
  sd_brood_chick  ~ exponential(2);
}

generated quantities {
  // ── Fixed effects back-transformed — unchanged ─────────────────────────────
  vector[P_chick] B_m_chick = R_inv_chick * B_mq_chick;
  vector[P_adult] B_m_adult = R_inv_adult * B_mq_adult;

  // ── Genetic CRN slopes back-transformed from QR scale ─────────────────────
  // B_cpc_A: genetic correlation reaction norm — analogue of old B_cpc
  // B_v_A:   genetic variance reaction norm — analogue of old B_v
  // Interpretation identical to phenotypic model but now genetic only
  matrix[P_y, ncor] B_cpc_A;
  matrix[P_y, D]    B_v_A;

  // ── Non-genetic individual CRN slopes back-transformed ────────────────────
  // B_cpc_E: non-genetic correlation reaction norm — new parameter
  // B_v_E:   non-genetic variance reaction norm — new parameter
  matrix[P_y, ncor] B_cpc_E;
  matrix[P_y, D]    B_v_E;

  for (d in 1:ncor) {
    B_cpc_A[, d] = R_inv * B_cpcq_A[, d];
    B_cpc_E[, d] = R_inv * B_cpcq_E[, d];
  }
  for (d in 1:D) {
    B_v_A[, d] = R_inv * B_vq_A[, d];
    B_v_E[, d] = R_inv * B_vq_E[, d];
  }
}