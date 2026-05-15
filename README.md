<img width="751" height="516" alt="dens_overlay plot" src="https://github.com/user-attachments/assets/cbaa4d8d-f084-45d7-b5d6-c3128afd6fb2" />
# Bayesian Regression on Poverty Income Ratio

*Lui Yiu Wa — 2026-01-06*

## Introduction

The data used in this analysis comes from the **National Health and Nutrition Examination Survey (NHANES)** conducted from August 2021 to August 2023.

This project aims to run a regression to model how demographic features affect the family monthly poverty income ratio — defined as "the ratio of monthly family income to the HHS poverty guidelines specific to family size" (NHANES). Specifically, it models the household reference person's age, gender, and education's effects on the poverty income ratio, adjusted for the number of people in the household.

> **Disclaimer:** I self-taught most of this and there is still an ocean of things to learn, so if you spot any errors along the way — that's probably the reason. With that in mind, we shall begin!

---

## Loading the Libraries and the Data

```r
library(haven)
library(tidyverse)
library(mice)
library(miceadds)
library(brms)
library(posterior)
library(bayestestR)
library(broom.mixed)
library(bayesplot)
```

```r
demographic <- read_xpt("DEMO_L.xpt")
income      <- read_xpt("INQ_L.xpt")

df <- inner_join(demographic, income)
## Joining with `by = join_by(SEQN)`

df <- df %>% select(
  ID            = SEQN,
  num_ppl       = DMDHHSIZ,
  HH_gender     = DMDHRGND,
  HH_age        = DMDHRAGZ,
  HH_edu        = DMDHREDZ,
  poverty_index = INDFMMPI
)

df <- df %>% mutate(
  HH_gender = as.factor(HH_gender),
  HH_edu    = as.factor(HH_edu)
)
```

---

## Handling Missing Data

We check for the presence of missing data:

```r
ini <- mice(df, maxit = 0)
ini$nmis
```

```
##           ID     num_ppl   HH_gender      HH_age      HH_edu poverty_index
##            0           0        7818        7809        8187          2944
```

```r
md.pattern(df)
```

```
##      ID num_ppl poverty_index HH_age HH_gender HH_edu  
## 3012  1       1             1      1         1      1   0
##   76  1       1             1      1         1      0   1
##    7  1       1             1      1         0      1   1
## 5894  1       1             1      0         0      0   3
##  727  1       1             0      1         1      1   1
##  300  1       1             0      1         1      0   2
##    2  1       1             0      1         0      0   3
## 1915  1       1             0      0         0      0   4
##       0       0          2944   7809      7818   8187 26758
```

And of course there is missing data — can't have it too easy!

With missing data, we have 2 popular approaches: multiple imputation and complete case analysis. If we can assume a **"Missing at Random" (MAR)** structure, multiple imputation is preferred as it avoids bias and discards less information. So, let us boldly assume MAR for a second (we will come back to this) and perform multiple imputation:

```r
if (FALSE) {
  df_imp <- futuremice(df, m = 20, parallelseed = 555, method = "rf",
                       maxit = 10, ntree = 25, verbose = T, n.core = 14)
  saveRDS(df_imp, "imputed_data.RDS")
}

df_imp      <- readRDS("imputed_data.RDS")
complete_df <- complete(df_imp, "long", include = T)
```

We can compare the distributions of `poverty_index` before and after imputation as a crude sanity check. They looked good — the imputed values were not very off.

---

## Model: Zero-One Inflated Beta

The distribution of `poverty_index` is highly irregular. NHANES censors all values beyond 5, and a ratio cannot drop below 0, so there are two point masses at 0 and 5. After scaling (`poverty_index / 5`), it is reasonable to use a **zero-one inflated beta (ZOIB) model**:

1. Each data point has a probability $\pi$ to be at the point masses (0 or 1) — the **zero-one inflated part** (zoi).
2. Given it is at a point mass, it has a probability $\lambda$ to be exactly 1 — the **conditional-one inflated part** (coi).
3. Otherwise, it follows a **beta distribution** with location and scale parameters $(\mu, \phi)$.

```r
complete_df <- complete_df %>% mutate(
  poverty_index_scaled = poverty_index / 5,
  HH_edu    = as.factor(HH_edu),
  HH_gender = as.factor(HH_gender)
)

complete_df <- as.mids(complete_df)
```

---

## Model Fitting with Hamiltonian Monte Carlo

We define the full model formula for `brms`:

```r
formulas <- bf(
  poverty_index_scaled ~ num_ppl + HH_gender + HH_age + HH_edu,
  phi ~ num_ppl + HH_gender + HH_age + HH_edu,
  zoi ~ num_ppl + HH_gender + HH_age + HH_edu,
  coi ~ num_ppl + HH_gender + HH_age + HH_edu
)

if (FALSE) {
  # The old gods have forsaken these chains.
  # The frogs are leaping and leaping and leaping...
  imp_fit <- brm_multiple(formulas, data = complete_df,
                          family = zero_one_inflated_beta(),
                          chains = 4, iter = 2000, cores = 4,
                          seed = 555)
  saveRDS(imp_fit, "impfit.RDS")
}

imp_fit <- read_rds("impfit.RDS")
```

---

## MCMC Diagnostics

After the excruciating Hamiltonian physics-fueled black magic, we perform diagnostics to ensure chains converged and we have enough valid samples. This requires careful handling as the model was fitted on multiply imputed datasets.

```r
# Modified from brms documentation
draws <- as_draws_array(imp_fit)
m     <- 20

draws_per_dat <- lapply(1:m, \(i) subset_draws(draws, chain = i))
fit_check     <- lapply(draws_per_dat, summarise_draws, default_convergence_measures())

problematic_list <- list()

for (i in seq_along(fit_check)) {
  problematic_rows <- fit_check[[i]] %>%
    filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400)

  if (nrow(problematic_rows) > 0) {
    problematic_list[[paste0("dataset_", i)]] <- problematic_rows %>%
      mutate(dataset = i) %>%
      select(dataset, everything())
  }
}

all_problematic <- bind_rows(problematic_list)
```

**Result:** 15 problematic parameters were found across 11 datasets, but they are either borderline acceptable (Rhat ≈ 1.01) or are `lp__`, the log-posterior, which is known to have difficult geometry for HMC to sample from. We can safely proceed.

---

## Posterior Predictive Checks

We check whether our model can capture the distribution of `poverty_index_scaled` via simulation (using one imputed dataset due to computational constraints):

```r
if (FALSE) {
  pp_check(imp_fit, type = "dens_overlay", ndraws = 500, set.seed(555))
  pp_check(imp_fit, type = "stat", stat = "mean", ndraws = 500, set.seed(555))
  pp_check(imp_fit, type = "ecdf_overlay", ndraws = 500, set.seed(555))
}
```

**Density overlay** — comparing observed $y$ vs. simulated $y_{\text{rep}}$:

![Density overlay posterior predictive check](https://raw.githubusercontent.com/cipher6120-bit/Bayesian-regression-on-family-income/main/dens_overlay-plot.jpg)

**Empirical CDF overlay:**

![ECDF overlay posterior predictive check](https://raw.githubusercontent.com/cipher6120-bit/Bayesian-regression-on-family-income/main/ecdf_overlay-plot.jpg)

**Simulated mean statistic** — the true mean $T(y)$ falls comfortably within the predictive distribution $T(y_{\text{rep}})$:

![Mean statistic posterior predictive check](https://raw.githubusercontent.com/cipher6120-bit/Bayesian-regression-on-family-income/main/mean_stat-plot.jpg)

Not too horrible if I do say so myself! The model captures the dynamic of the response distribution well.

---

## Model Summary and Interpretation

```r
imp_fit_broom <- tidyMCMC(imp_fit, conf.int = T)
print(imp_fit_broom, n = Inf)
```

```
## # A tibble: 29 x 5
##    term                 estimate std.error conf.low conf.high
##    <chr>                   <dbl>     <dbl>    <dbl>     <dbl>
##  1 b_Intercept            -0.563   0.0931   -0.760    -0.387
##  2 b_phi_Intercept         1.49    0.142     1.24      1.75
##  3 b_zoi_Intercept        -1.80    0.279    -2.30     -1.21
##  4 b_coi_Intercept        -0.938   0.621    -2.22      0.227
##  5 b_num_ppl              -0.128   0.00902  -0.144    -0.109
##  6 b_HH_gender2           -0.193   0.0338   -0.256    -0.126
##  7 b_HH_age                0.139   0.0330    0.0808    0.208
##  8 b_HH_edu2               0.262   0.0523    0.161     0.355
##  9 b_HH_edu3               1.01    0.0829    0.795     1.15
## 10 b_phi_num_ppl           0.0940  0.0127    0.0692    0.118
## ... (29 total rows)
```

We focus on the effect of `HH_gender2` (female reference person), `HH_age`, and `HH_edu3` (college graduate or above) on $\mu$ and `coi`:

```r
if (FALSE) {
  mcmc_areas(imp_fit, pars = c("b_HH_age", "b_coi_HH_age",
                               "b_HH_gender2", "b_coi_HH_gender2",
                               "b_HH_edu3", "b_coi_HH_edu3"),
             prob = 0.95) +
    ggtitle("Posterior Distributions of selected model parameters")
}
```

![Posterior Distributions of selected model parameters](https://raw.githubusercontent.com/cipher6120-bit/Bayesian-regression-on-family-income/main/Posterior-Distributions-of-selected-model-parameters.jpg)

Key findings:

1. **Higher age** of the household reference person is associated with a **higher** scaled poverty ratio.
2. **Female** reference persons are associated with a **lower** scaled poverty ratio.
3. **College graduates or above** are associated with a **higher** scaled poverty ratio.

It is also worth noting that the posterior distribution of `b_HH_edu3` is bi- or even tri-modal — a subtlety that frequentist modeling would fail to capture, as it relies heavily on asymptotic normality assumptions.

---

## Sensitivity Analysis on the MAR Assumption

The MAR assumption necessary for valid multiple imputation cannot be verified statistically — we can only argue for it using domain knowledge. I suspect respondents with lower poverty ratios may be more likely to not respond (possibly due to embarrassment), suggesting the mechanism may actually be **Missing Not At Random (MNAR)**.

We simulate what would happen under MNAR using the **delta adjustment method**: artificially adding a bias term $\delta$ to missing values during imputation, adapted from Stef Van Buuren's *Flexible Imputation of Missing Data*:

```r
delta <- c(-0.5)
post  <- ini$post

if (FALSE) {
  for (i in 1:length(delta)) {
    d   <- delta[i]
    cmd <- paste("imp[[j]][,i] <- imp[[j]][,i] +", d)
    post["poverty_index"] <- cmd
    sens_imp <- mice(df, method = "rf", m = 3, post = post,
                     maxit = 10, seed = 555)
  }
  saveRDS(sens_imp, "sens_imp.RDS")
}

sens_imp <- readRDS("sens_imp.RDS")
```

After re-fitting and comparing estimates side by side:

```
## # A tibble: 29 x 3
##    term               estimate_imp estimate_delta_sens
##  1 b_Intercept              -0.563              -0.791
##  2 b_HH_gender2             -0.193              -0.219
##  3 b_HH_age                  0.139               0.170
##  4 b_HH_edu3                 1.01                1.27
##  ...
## 18 b_zoi_edu2                0.232              -0.0602   # direction changed!
## 23 b_coi_HH_edu2            -0.0918              1.20    # direction changed!
```

The directions of most model parameters are consistent, implying that if the systematic bias is $\leq -0.5$ on the unscaled `poverty_index`, our conclusions will remain roughly the same. However, `b_zoi_edu2` and `b_coi_HH_edu2` reversed direction, so our conclusions may **not** be entirely robust under a misspecified missing mechanism.

> Ideally, we would explore more extreme bias values with more imputed datasets — but HMC takes a long time to run and I was getting grumpier and grumpier as it ran, so I'll call it here.

---

## Informative Prior Sensitivity (Supplementary)

We also explored the robustness of the results to an informative prior specification:

```r
priors <- c(set_prior("normal(-1, 5)", class = "Intercept", dpar = "mu"),
            set_prior("normal(-1, 5)", class = "Intercept", dpar = "coi"))
```

Fitting on a subset of imputed datasets (2 out of 20 for speed), the direction of all effects remained consistent and estimates were very close to those from the non-informative prior model — confirming that the prior choice did not substantively influence inference, which is reassuring given the large sample size.

---

## Conclusion

In this project, I performed **Bayesian regression on a multiply imputed NHANES dataset** using `brms` (and Stan under the hood) with:

- A zero-one inflated beta model to handle the censored, bounded response
- Full MCMC diagnostics (Rhat, ESS) adapted for multiply imputed data
- Posterior predictive checks for model fit assessment
- Sensitivity analysis for the MAR assumption via delta adjustment
- Informative prior sensitivity analysis

I am still working on my skills in summarizing complex models and data visualization — so pardon if it's not the best. But this basically wraps up the project!

---

*I am not a statistician — at least, not yet. Take all of this with appropriate skepticism and a grain of salt.*
