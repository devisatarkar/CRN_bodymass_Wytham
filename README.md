# **Strength of developmental carry-over varies continuously with environmental conditions in a wild bird population**

Data and R/Stan code for reproducing the analyses and figures in "***Strength of developmental carry-over varies continuously with environmental conditions in a wild bird population"***, examining developmental carry-over effects of body mass in great tits (*Parus major*) at Wytham Woods, Oxford, using the Covariance Reaction Norm (CRN) framework.

All models run using:
R v.4.3.3 https://www.r-project.org/,
CmdStanR v.0.9.0 https://mc-stan.org/cmdstanr/

## File list

-   `Source Data:`contains files needed to run the R scripts for analysis. Main files include:

1.  df_chick_p.csv — 4345 x 13 table, each row is an individual nestling, with natal information
2.  df_adult_p.csv — 6362 x 15 table, each row is a breeding attempt by an individual, with breeding attempt information
3.  gtbt_density_peryear.csv — population density (number of breeding pairs of great tits and blue tits) per year
4.  ped_pruned.csv — Pruned pedigree used to create relatedness matrix, with id, dam, and sire columns

-   `Analyses (R scripts):`contains R scripts needed to run the CRN models:

1.  DATA PREP.R — Uses `df_chick_p.RDS` / `df_adult_p.RDS` as inputs and builds the context key, context matrix, and stan data list for the CRN model in Stan

2.  Part 1 CRN MODEL.R — Main analysis R script for Part 1. Uses Stan file `crn_wytham_pairedcontext.stan`. Output (posterior draws) stored in `posterior_draws_part1.RDS`. Also contains code to reproduce Figure 2 using posterior draws.

3.  Part 2 CRN MODEL.R — Main analysis R script for Part 2. Uses Stan file `crn_wytham_pairedcontext.stan`. Output (posterior draws) stored in `posterior_draws_mismatch.RDS`. Also contains code to reproduce Figure 3 using posterior draws.

4.  Part 3 GENETIC CRN MODEL.R — Main analysis R script for Part 3. Uses Stan file `crn_wytham_genetic+nongenetic.stan`. Output (posterior draws) stored in `posterior_draws_genetic.RDS`. Also contains code to reproduce Figure 4 using posterior draws, create relatedness matrix using `ped_pruned.RDS`, and calculate heritability at mean environmental conditions.

-   `Models (Stan):` contains the main CRN statistical models:

1.  `crn_wytham_pairedcontext.stan` — Stan model shared by Parts 1-2.
2.  `crn_wytham_genetic+nongenetic.stan` — Stan model for Part 3; adds the genetic/non-genetic (`_A`/`_E`) decomposition and the pedigree relatedness matrix `A` as a data input, on top of the same phenotypic structure as the Part 1-2 model.

-   `Outputs:` contains the posterior draws generated after running all models. Also contains the Stan data list and context key generated in DATA PREP.R, which is used as an input for the Analyses scripts. Following are the files in this folder:

1.  `stan_dl_part1.RDS` — generated in DATA PREP.R
2.  `context_key.RDS` — generated in DATA PREP.R
3.  `posterior_draws_part1.RDS` — output (posterior draws) of Part 1 model
4.  `posterior_crn_part1.RDS`— output (posterior draws) of Part 1 model; CRN specific parameters
5.  `posterior_draws_mismatch.RDS` — output (posterior draws) of Part 2 model
6.  `posterior_crn_mismatch.RDS`— output (posterior draws) of Part 2 model; CRN specific parameters
7.  `posterior_draws_genetic.RDS` — output (posterior draws) of Part 3 model
8.  `posterior_crn_genetic.RDS`— output (posterior draws) of Part 3 model; CRN specific parameters

Posterior draws can be used to reproduce figures (Code available at the end of each analysis script).

## Data Information

**Following are descriptions of the data files provided with the code. Only data explicitly required for running all the models have been provided. These data come from the long-term individual-based study of great tits in Wytham Woods, Oxfordshire, UK.**

1.  df_chick_p.RDS (N = 4345 observations)

- id_num: unique numeric index for each individual (runs from 1 to N)
- id: unique identifier for each individual
- birthyear: year in which individual hatched
- birthyr_num: birthyear-specific numeric index
- broodid: unique identifier for the brood an individual belongs to
- brood_num: brood-specific numeric index
- chickweight: weight in grams of individual at 15 days old (response variable)
- april_birthlaydate: lay date (clutch initiation date) in april days (number of days since april 1st)
- natalbroodsize: number of hatched chicks in brood
- natal_meantemp: Average daily temperature in ºC during developmental period (0-15 days old)
- natal_meanrain: Average daily rainfall in mm during developmental period (0-15 days old)
- absnatal_hd_mismatch: absolute difference between the half-fall date of winter moth larvae and the tenth day post-hatching of the focal individual
- context_id: which context (birthyear x breeding year) the individual belongs to

2.  df_adult_p.RDS (6362 observations)

- id_num: unique numeric index for each individual from df_chick_p
- id: unique identifier for each individual (these ids were nestlings in df_chick_p, and each row in df_adult_p corresponds to their subsequent breeding attempts)
- breedyear: year in which breeding attempt has occurred
- breedyr_num: breeding year-specific numeric index
- Pnum: unique identifier for each breeding attempt
- breed_nestbox: nestbox in which breeding attempt has occurred
- breedbox_num: nestbox-specific numeric index
- adultweight: weight in grams of individual during breeding attempt (response variable)
- april_laydate: lay date (clutch initiation date of brood produced by focal individual) in april days (number of days since april 1st)
- breed_meantemp: Average daily temperature in ºC during developmental period of offspring (0-15 days old)
- breed_meanrain: Average daily rainfall in mm during developmental period of offspring (0-15 days old)
- num_chicks: number of hatched chicks in brood produced in breeding attempt
- age_years: age of breeding individual
- sex_num: sex of breeding individual (0 if female, 1 if male)
- context_id: which context (birthyear x breeding year) the individual belongs to

3.  gtbt_density_peryear.RDS (66 observations)

- year: year of monitoring
- popdens_bt_gt: number of breeding pairs of great tits and blue tits

4.  ped_pruned.RDS (8633 observations)

- id: unique identifier for each individual
- dam: unique identifier for the mother
- sire: unique identifier for the social father
