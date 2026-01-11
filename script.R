library(haven)
library(tidyverse)
library(mice)
library(miceadds)
library(brms)
library(posterior)
library(bayestestR)
library(broom.mixed)
library(bayesplot)
demographic <- read_xpt("DEMO_L.xpt")
demographic
income <- read_xpt("INQ_L.xpt")
income

df <- inner_join(demographic, income)
df

df <- df %>% select(
 ID = SEQN,
 num_ppl = DMDHHSIZ,
 HH_gender = DMDHRGND,
 HH_age = DMDHRAGZ,
 HH_edu = DMDHREDZ,
 poverty_index = INDFMMPI
)
df %>% mutate(
  HH_gender = as.factor(HH_gender),
  HH_edu = as.factor(HH_edu)
)
#df_imp <- futuremice(df, m = 20, parallelseed = 555, method = "rf", 
#                     maxit = 10, ntree = 25, 
#                     verbose = T, n.core = 14)
#saveRDS(df_imp, "imputed_data.RDS")
df_imp <- readRDS("imputed_data.RDS")


complete_df <- complete(df_imp, "long", include = T)

#plot poverty_index from complete_df
ggplot(complete_df, aes(x = poverty_index)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Poverty Index",
       x = "Poverty Index",
       y = "Frequency") +
  theme_minimal()

#plot poverty_index from df
ggplot(df, aes(x = poverty_index)) +
  geom_histogram(binwidth = 0.1, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Poverty Index (Original Data)",
       x = "Poverty Index",
       y = "Frequency") +
  theme_minimal()


complete_df <- complete_df %>% mutate(
  poverty_index_scaled = poverty_index / 5,
  HH_edu = as.factor(HH_edu),
  HH_gender = as.factor(HH_gender)
)

complete_df <- as.mids(complete_df)

formulas <- bf(poverty_index_scaled ~ num_ppl + 
                 HH_gender + HH_age + HH_edu,
               phi ~ num_ppl + HH_gender + HH_age + HH_edu,
               zoi ~ num_ppl + HH_gender + HH_age + HH_edu,
               coi ~ num_ppl + HH_gender + HH_age + HH_edu
               )

#imp_fit <- brm_multiple(formulas, data = complete_df,
#                       family = zero_one_inflated_beta(),
#                        chains = 4, iter = 2000, cores = 4,
#                        seed = 555)
                        

#The old gods have forsaken these chains. 
#The frogs are leaping and leaping and leaping and leaping... 
#All they have got is the remembrance of the conjugate priors. 
#The Hamiltonian paths of the spheres tumbling down the wicked hills of the Log-posterior-Stan
#fluttering in numerical fluctuations.
#All is lost, all is lost...

#saveRDS(imp_fit, "impfit.RDS")
imp_fit <- read_rds("impfit.RDS")

draws <- as_draws_array(imp_fit)
m <- 20

draws_per_dat <- lapply(1:m, \(i) subset_draws(draws, chain = i))
fit_check <- lapply(draws_per_dat, summarise_draws, default_convergence_measures())

#check if any rhat >1.01 or ess_bulk < 400 or ess_tail < 400 in fit_check, return which row is problematic
# Initialize a list to store problematic rows for each dataset
problematic_list <- list()

# Check each dataset
for(i in seq_along(fit_check)) {
  problematic_rows <- fit_check[[i]] %>%
    filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400)
  
  if(nrow(problematic_rows) > 0) {
    problematic_list[[paste0("dataset_", i)]] <- problematic_rows %>%
      mutate(dataset = i) %>%
      select(dataset, everything())
  }
}

all_problematic <- bind_rows(problematic_list)

if(nrow(all_problematic) > 0) {
  cat("Found", nrow(all_problematic), "problematic parameters across", 
      length(problematic_list), "datasets:\n\n")
  print(all_problematic, n = nrow(all_problematic))
} else {
  cat("No problematic parameters found in any dataset 
      (using criteria: rhat > 1.01 or ess_bulk < 400 or ess_tail < 400).\n")
}

pp_check(imp_fit, type = "dens_overlay", ndraws = 500, set.seed(555))
pp_check(imp_fit, type = "stat", stat = "mean", ndraws = 500, set.seed(555))
pp_check(imp_fit, type = "ecdf_overlay", ndraws = 500, set.seed(555))

mcmc_areas(imp_fit, pars = c("b_HH_age", "b_coi_HH_age",
                             "b_HH_gender2","b_coi_HH_gender2",
                             "b_HH_edu3", "b_coi_HH_edu3"),
            prob = 0.95) +
  ggtitle("Posterior Distributions of selected model parameters")

#not horrible!

imp_fit_broom <- tidyMCMC(imp_fit, conf.int = T)
print(imp_fit_broom, n = Inf)
ini <- mice(df, maxit = 0)
ini$nmis
md.pattern(df)

#sensitivity checking for priors


priors <- c(set_prior("normal(-1,5)", class = "Intercept", dpar = "mu"),
            set_prior("normal(-1,5)", class = "Intercept", dpar = "coi"))


partial_complete_df <- subset_datlist(
  datlist = datlist_create(complete_df),  
  index = c(4,6)
) %>% datlist2mids()


#prior_sens_fit <- brm_multiple(formulas, data = partial_complete_df,
#                       family = zero_one_inflated_beta(),
#                       prior = priors,
#                        chains = 4, iter = 2000, cores = 4,
#                        seed = 555)
#saveRDS(prior_sens_fit, "prior_sens_fit.RDS")

prior_sens_fit <- readRDS("prior_sens_fit.RDS")

draws <- as_draws_array(prior_sens_fit)
m <- 5

draws_per_dat <- lapply(1:m, \(i) subset_draws(draws, chain = i))
prior_sens_fit_check <- lapply(draws_per_dat, summarise_draws, default_convergence_measures())
imp_fit_broom <- tidyMCMC(imp_fit)
sens_fit_broom <- tidyMCMC(prior_sens_fit)


imp_fit_renamed <- imp_fit_broom %>%
  rename_with(~ paste0(., "_imp"), -term)

sens_fit_renamed <- sens_fit_broom %>%
  rename_with(~ paste0(., "_sens"), -term)

combined_table <- imp_fit_renamed %>%
  left_join(sens_fit_renamed, by = "term")


print(combined_table, n = Inf)

#Stef Van Burren code
delta <- c(-0.5)

post <- ini$post
imp.all.undamped <- vector("list", length(delta))

#for (i in 1:length(delta)) {
#  d <- delta[i]
#  cmd <- paste("imp[[j]][,i] <- imp[[j]][,i] +", d)
#  post["poverty_index"] <- cmd
#  sens_imp <- mice(df,  method = "rf", m = 3, post = post, maxit = 10,
#              seed = 555)
#}

#saveRDS(sens_imp, "sens_imp.RDS")

sens_imp <- readRDS("sens_imp.RDS")
sens_complete_df <- complete(sens_imp, "long", include = T)

sens_complete_df <- sens_complete_df %>%
  mutate(poverty_index = ifelse(.imp != 0 & poverty_index < 0, 0, poverty_index))

hist(sens_complete_df$poverty_index)

sens_complete_df <- sens_complete_df %>% mutate(
  poverty_index_scaled = poverty_index / 5,
  HH_edu = as.factor(HH_edu),
  HH_gender = as.factor(HH_gender)
)

sens_complete_df <- as.mids(sens_complete_df)

#delta_sens_fit <- brm_multiple(formulas, data = sens_complete_df,
#                       family = zero_one_inflated_beta(),
#                        chains = 4, iter = 2000, cores = 4,
#                        seed = 555)

#saveRDS(delta_sens_fit, "delta_sens_fit.RDS")
delta_sens_fit <- readRDS("delta_sens_fit.RDS")

draws <- as_draws_array(delta_sens_fit)
m <- 3

draws_per_dat <- lapply(1:m, \(i) subset_draws(draws, chain = i))
delta_sens_fit_check <- lapply(draws_per_dat, summarise_draws, default_convergence_measures())

delta_sens_fit_check
view(delta_sens_fit_check)

delta_sens_fit_broom <- tidyMCMC(delta_sens_fit)



delta_sens_fit_renamed <- delta_sens_fit_broom %>%
  rename_with(~ paste0(., "_delta_sens"), -term)

combined_table <- imp_fit_renamed %>%
  left_join(sens_fit_renamed, by = "term")
print(combined_table, n = Inf)


