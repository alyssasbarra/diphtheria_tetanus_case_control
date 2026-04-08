#### diphtheria OG model best

library(tidyverse)
library(ggplot2)
library(readxl)
library(lubridate)
library(haven)
library(janitor)
library(sf)
library(data.table)
library(epitools)
library(survival)
library(rstan)
library(ggbeeswarm)
library(scales)
library(rstan)
library(loo)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


`%!in%` = Negate(`%in%`)

## -----------------------------------------------------------
## Read in PIRMZ final data 
dat_pirmz <- read_dta("PIRMZ data with wts.dta")

dat_pirmz %>%
  mutate(province = as_factor(province)) %>%
  mutate(district = as_factor(district)) %>%
  mutate(agecat = as_factor(agecat)) %>%
  mutate(gender = as_factor(gender)) %>%
  mutate(hivstatusfinal = as_factor(hivstatusfinal)) %>%
  mutate(immunity_m = as_factor(immunity_m)) %>%
  mutate(immunity_m2 = as_factor(immunity_m2)) %>%
  mutate(immunity_r = as_factor(immunity_r)) %>%
  mutate(sampletype = as_factor(sampletype)) -> dat_pirmz

dat_pirmz <- dat_pirmz %>% dplyr::select(-contains("wt_pirmz_jkn"))


levels(dat_pirmz$immunity_m2) <- levels(dat_pirmz$immunity_m)
dat_pirmz$immunity_m2[which(is.na(dat_pirmz$final_quant_m))] <- "missing"

## -----------------------------------------------------------
## Read in data to link PIRMZ to samples tested 

dat_link <- read_xlsx('Updated list of retrieved DipTox Samples _2024-01_29-928samples.xlsx')
dat_link <- dat_link %>% dplyr::select("ZM (ID 1)", "ID (Global Spec)")
colnames(dat_link) <- c('Sample', 'pirmz_id')

## doesn't have barcodes yet -- to discuss with Andrea / Saki
dat_link2 <- fread('67 of the 168 proposed_samples_DipTox_secondary_list_20240219.csv')
dat_link3 <- fread('96 of 168 proposed_samples_DipTox_secondary_list_20240215_retrieved.csv')

dat_link_all <- fread('~/Downloads/PIMRV sample status report (1)_v2.csv', header = F)
colnames(dat_link_all) <- c('error_pirmz_id', 'Sample')
 
dat_link_all$pirmz_id <- gsub("\036", "-", dat_link_all$error_pirmz_id)
dat_link_all$error_pirmz_id <- NULL

## -----------------------------------------------------------
## Read in results from diphtheria testing
dat.conc.d <- data.table()

for(i in 1:17){
  load(paste0('/Plate',i,'_Diphtheria_flexfit_conc.RData'))
  dat.all.conc_copy <- copy(dat.all.conc)
  dat.all.conc_copy$batch <- i
  dat.conc.d <- rbind(dat.conc.d, dat.all.conc_copy)
}

dat.conc.d <- subset(dat.conc.d, Sample %!in% c('Cal 1', 'Cal 2', 'Cal 3', 'Cal 4',  
                                                'Pos Ctrl', 'Neg Ctrl', 'BLANK', 
                                                'IC High (IC4)', 'IC Low (IC2)', 'IC1', 'IC5', 'IC7'))

dat.conc.d <- subset(dat.conc.d, select=c("Sample", "OD", "concentration", "batch"))
colnames(dat.conc.d) <- c("Sample", "OD_d", "concentration_d", "batch")

dat.conc.d.final <- data.table()
dat.conc.d.final_dup <- data.table()
for(i in 1:length(unique(dat.conc.d$Sample))){
  
  dat.conc.d_i <- subset(dat.conc.d, Sample == unique(dat.conc.d$Sample)[i])
  
  if(dim(dat.conc.d_i)[1] >1){
    message()
    dat.conc.d_i <- subset(dat.conc.d_i, batch == min(dat.conc.d_i$batch))
    dat.conc.d_i_dup <- subset(dat.conc.d_i, batch == max(dat.conc.d_i$batch))
    
    dat.conc.d.final <- rbind(dat.conc.d.final, dat.conc.d_i)
    dat.conc.d.final_dup <- rbind(dat.conc.d.final_dup, dat.conc.d_i_dup)
    
    
  }else{
    dat.conc.d.final <- rbind(dat.conc.d.final, dat.conc.d_i)
  }
}



## -----------------------------------------------------------
## Read in results from tetanus testing
dat.conc.t <- data.table()

for(i in 1:18){
  load(paste0('/Plate',i,'_Tetanus_flexfit_conc.RData'))
  dat.all.conc_copy <- copy(dat.all.conc)
  dat.all.conc_copy$batch <- i
  dat.conc.t <- rbind(dat.conc.t, dat.all.conc_copy)
  
}

dat.conc.t <- subset(dat.conc.t, Sample %!in% c('Cal 1', 'Cal 2', 'Cal 3', 'Cal 4', 'Cal 5', 'Cal 6', 'Cal 7', 
                                                'High Ctrl', 'Low Ctrl', 'BLANK', 
                                                'IC High (IC3 old)', 'IC Low (IC6)'))

dat.conc.t <- subset(dat.conc.t, select=c("Sample", "OD", "concentration", "batch"))
colnames(dat.conc.t) <- c("Sample", "OD_t", "concentration_t", "batch")

dat.conc.t.final <- data.table()
dat.conc.t.final_dup <- data.table()
for(i in 1:length(unique(dat.conc.t$Sample))){
  
  dat.conc.t_i <- subset(dat.conc.t, Sample == unique(dat.conc.t$Sample)[i])
  
  if(dim(dat.conc.t_i)[1] >1){
    message(unique(dat.conc.t$Sample)[i])
    dat.conc.t_i_dup <- subset(dat.conc.t_i, batch == max(dat.conc.t_i$batch))
    dat.conc.t_i <- subset(dat.conc.t_i, batch == min(dat.conc.t_i$batch))

    
    dat.conc.t.final <- rbind(dat.conc.t.final, dat.conc.t_i)
    dat.conc.t.final_dup <- rbind(dat.conc.t.final_dup, dat.conc.t_i_dup)
    
  }else{
    dat.conc.t.final <- rbind(dat.conc.t.final, dat.conc.t_i)
  }
}

## -----------------------------------------------------------
## Check out replicates

head(dat.conc.t.final_dup)

dat.conc.t.final_dup$concentration_t_rep <- dat.conc.t.final_dup$concentration_t
dat.conc.t.final_dup$OD_t_rep <- dat.conc.t.final_dup$OD_t
dat.conc.t.final_dup$batch_rep <- dat.conc.t.final_dup$batch

dat.conc.t.final_dup$concentration_t <- NULL
dat.conc.t.final_dup$OD_t <- NULL
dat.conc.t.final_dup$batch <- NULL


dat.conc.t.final_dup_test <- merge(dat.conc.t.final_dup, dat.conc.t.final, all.x=T, all.y=F)

head(dat.conc.t.final_dup_test)

## ------------------------------------------------------
## Link results to PIRMZ samples
dat.conc <- merge(dat.conc.t.final, dat.conc.d.final, by='Sample', all.x=T, all.y=T)

df <- merge(dat_link_all, dat.conc, by='Sample', all.y=T, all.x=F)
df$id <- gsub("\\-.*","",df$pirmz_id)
df$pirmz_id <- NULL
df <- unique(df)
df <- merge(df, dat_pirmz, by='id', all.x=T)

missing <- subset(df, is.na(immunity_m2))
missing[1:10,1:5] # 2 barcodes missing results TODO:: fix

df <- subset(df, !is.na(immunity_m2))

df$immunity_m_binary <- ifelse(df$immunity_m2 %in% c('negative'), 'negative', 'positive')

## -----------------------------------------------------------
## Use diphtheria thresholds to determine seropositivity  

df$immunity_d <- ifelse(df$concentration_d< 0.01, 'negative',
                        ifelse(df$concentration_d >= 0.01 & df$concentration_d < 0.1, 'uncertain',
                               ifelse(df$concentration_d >= 0.1 & df$concentration_d < 1, 'positive',
                               ifelse(df$concentration_d >= 1, 'high positive', NA))))

#df$immunity_d_binary <- ifelse(df$immunity_d %in% c('high positive', 'positive'), 'positive', 'negative')
df$immunity_d_binary <- ifelse(df$immunity_d %in% c('high positive', 'positive', 'uncertain'), 'positive', 'negative')

## -----------------------------------------------------------
## Use tetanus thresholds to determine seropositivity  

df$immunity_t_kit <- ifelse(df$concentration_t < 0.01, 'negative',
                            ifelse(df$concentration_t >= 0.01 & df$concentration_t < 0.15, 'border',
                        ifelse(df$concentration_t >= 0.15, 'positive', NA)))

df$immunity_t <- ifelse(df$concentration_t < 0.2, 'negative',
                        ifelse(df$concentration_t >= 0.2, 'positive', NA))

df$immunity_t_binary <- df$immunity_t

## -----------------------------------------------------------
## Determine how many samples tested so far are in matched pairs 

dat_pairs <- fread('proposed_samples_DipTox_20231103_rearranged-1326samples.csv')
dat_pairs <- subset(dat_pairs, select = c('id', 'pair_number'))

df <- merge(df, dat_pairs, by='id', all.x=T)
df <- subset(df, !is.na(pair_number) & !is.na(concentration_t) & !is.na(concentration_d))

complete_tested_pairs <- df$pair_number[which(duplicated(df$pair_number))]
df_complete <- subset(df, pair_number %in% c(complete_tested_pairs))

df_not_complete <- subset(df, pair_number %!in% c(complete_tested_pairs))

ftable(df_complete$immunity_t_binary ~ df_complete$immunity_m_binary)
ftable(df_complete$immunity_d_binary ~ df_complete$immunity_m_binary)

ftable(df_complete$immunity_d_binary ~ df_complete$immunity_t_binary)

df_complete$immunity_d_binary2 <- ifelse(df_complete$immunity_d_binary == 'positive', 1, 0)

df_complete %>% group_by(age) %>% summarise(mean=mean(immunity_d_binary2, na.rm=T)) -> test_data_for_stan


##### --------
#### stan model

df_complete$A_6yo <- ifelse(df_complete$age >= 6,1,0)
df_complete$A_7yo <- ifelse(df_complete$age >= 7,1,0)
df_complete$A_8yo <- ifelse(df_complete$age >= 8,1,0)
df_complete$A_9yo <- ifelse(df_complete$age >= 9,1,0)
df_complete$A_10yo <- ifelse(df_complete$age >= 10,1,0)


stan_code2 <- "
data {
  int<lower=1> N;                 // Number of observations
  vector[N] log_concentration_d;   // Response variable (log-transformed)
  vector[N] log_age;               // Predictor: log(age)
  vector[N] A_6yo;                 // Dummy variable for age 6
  vector[N] A_7yo;                 // Dummy variable for age 7
  vector[N] A_8yo;                 // Dummy variable for age 8
  vector[N] A_9yo;                 // Dummy variable for age 9
  vector[N] A_10yo;                // Dummy variable for age 10
}

parameters {
  real alpha;                      // Intercept
  real beta_age;                    // Slope for log(age)
  real beta_6yo;                    // Coefficient for age 6
  real beta_7yo;                    // Coefficient for age 7
  real beta_8yo;                    // Coefficient for age 8
  real beta_9yo;                    // Coefficient for age 9
  real beta_10yo;                   // Coefficient for age 10
  real<lower=0> sigma;              // Standard deviation of residuals
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta_age ~ normal(0, 10);
  beta_6yo ~ normal(0, 5);
  beta_7yo ~ normal(0, 5);
  beta_8yo ~ normal(0, 5);
  beta_9yo ~ normal(0, 5);
  beta_10yo ~ normal(0, 5);
  sigma ~ normal(0, 5);
  
  // Likelihood
  log_concentration_d ~ normal(
    alpha + beta_age * log_age + 
      beta_6yo * A_6yo + 
      beta_7yo * A_7yo + 
      beta_8yo * A_8yo + 
      beta_9yo * A_9yo + 
      beta_10yo * A_10yo, sigma);
}

generated quantities {
  vector[N] log_lik;         // Log-likelihood for each observation

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(log_concentration_d[i] | (alpha + beta_age * log_age[i] + beta_6yo * A_6yo[i] + beta_7yo * A_7yo[i] + 
                                      beta_8yo * A_8yo[i] + beta_9yo * A_9yo[i] + beta_10yo * A_10yo[i]), sigma);
  }
}

"

# Prepare data for Stan
stan_data <- list(
  N = nrow(df_complete),
  log_concentration_d = log(df_complete$concentration_d),
  log_age = log(df_complete$age),
  A_6yo = df_complete$A_6yo,
  A_7yo = df_complete$A_7yo,
  A_8yo = df_complete$A_8yo,
  A_9yo = df_complete$A_9yo,
  A_10yo = df_complete$A_10yo
)

set.seed(123)
# Compile the Stan model
stan_model <- stan(model_code = stan_code2, data = stan_data, 
                   iter = 10000, chains = 4)

log_lik1 <- extract_log_lik(stan_model)
waic1 <- waic(log_lik1)


# Extract log-likelihood samples
log_lik_samples <- extract(stan_model, pars = "log_lik")$log_lik  # Matrix of size (iterations x N)

# Step 1: Compute \bar{D}
mean_log_lik <- colMeans(log_lik_samples)  # Mean log-likelihood per observation
D_bar <- -2 * sum(mean_log_lik)

# Step 2: Compute D(\hat{\theta})
# Extract posterior means of parameters
posterior_means <- extract(stan_model, pars = c("alpha", "beta_age", "beta_6yo", "beta_7yo", "beta_8yo", "beta_9yo", "beta_10yo", "sigma"))
alpha_hat <- mean(posterior_means$alpha)
beta_age_hat <- mean(posterior_means$beta_age)
beta_6yo_hat <- mean(posterior_means$beta_6yo)
beta_7yo_hat <- mean(posterior_means$beta_7yo)
beta_8yo_hat <- mean(posterior_means$beta_8yo)
beta_9yo_hat <- mean(posterior_means$beta_9yo)
beta_10yo_hat <- mean(posterior_means$beta_10yo)
sigma_hat <- mean(posterior_means$sigma)

# Compute log-likelihood at posterior mean
log_lik_hat <- dnorm(stan_data$log_concentration_d, 
                     mean = alpha_hat + beta_age_hat * stan_data$log_age +
                       beta_6yo_hat * stan_data$A_6yo +
                       beta_7yo_hat * stan_data$A_7yo +
                       beta_8yo_hat * stan_data$A_8yo +
                       beta_9yo_hat * stan_data$A_9yo +
                       beta_10yo_hat * stan_data$A_10yo, 
                     sd = sigma_hat, 
                     log = TRUE)

D_hat <- -2 * sum(log_lik_hat)

# Step 3: Compute DIC
DIC <- 2 * D_bar - D_hat

print(DIC)


posterior_samples <- extract(stan_model)

n_samples <- length(posterior_samples$alpha)
new_age <- c(2,3,4,5,6,7,8,9,10)
# Initialize matrix to store predictions
predictions <- matrix(NA, n_samples, 9)

A_6yo <- ifelse(new_age >= 6,1,0)
A_7yo <- ifelse(new_age >= 7,1,0)
A_8yo <- ifelse(new_age >= 8,1,0)
A_9yo <- ifelse(new_age >= 9,1,0)
A_10yo <- ifelse(new_age >= 10,1,0)

# Loop through each posterior sample and calculate predicted y values
for (i in 1:n_samples) {
  for(a in 1:9){
    predictions[i, a] <- (posterior_samples$alpha[i] + (posterior_samples$beta_age[i] * log(new_age[a])) +
                            ((posterior_samples$beta_6yo[i]) * A_6yo[a]) +
                            ((posterior_samples$beta_7yo[i]) * A_7yo[a]) +
                            ((posterior_samples$beta_8yo[i]) * A_8yo[a]) +
                            ((posterior_samples$beta_9yo[i]) * A_9yo[a]) +
                            ((posterior_samples$beta_10yo[i]) * A_10yo[a]))
  }
}

quantile(posterior_samples$alpha, probs=c(0.5))
quantile(posterior_samples$beta_age, probs=c(0.5))
quantile(posterior_samples$beta_6yo, probs=c(0.5))
quantile(posterior_samples$beta_7yo, probs=c(0.5))
quantile(posterior_samples$beta_8yo, probs=c(0.5))
quantile(posterior_samples$beta_9yo, probs=c(0.5))
quantile(posterior_samples$beta_10yo, probs=c(0.5))

quantile(posterior_samples$alpha, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_age, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_6yo, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_7yo, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_8yo, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_9yo, probs=c(0.025, 0.975))
quantile(posterior_samples$beta_10yo, probs=c(0.025, 0.975))


predictions_counterfactual <- matrix(NA, n_samples, 9)

# Loop through each posterior sample and calculate predicted y values
for (i in 1:n_samples) {
  for(a in 1:9){
    predictions_counterfactual[i, a] <- (posterior_samples$alpha[i] + (posterior_samples$beta_age[i] * log(new_age[a]))) 
  }
}

# Mean predictions for each new x
predictions <- exp(predictions)
mean_predictions <- apply(predictions, 2, mean)

predictions_counterfactual <- exp(predictions_counterfactual)
mean_predictions_counterfactual <- apply(predictions_counterfactual, 2, mean)

# 95% credible interval (2.5th and 97.5th percentiles)
lower_bound <- apply(predictions, 2, quantile, 0.025)
upper_bound <- apply(predictions, 2, quantile, 0.975)

# Display predictions with intervals
pred_results <- data.frame(
  new_age = new_age,
  mean_prediction = (mean_predictions),
  lower_bound = (lower_bound),
  upper_bound = (upper_bound)
)

lower_bound_counterfactual <- apply(predictions_counterfactual, 2, quantile, 0.025)
upper_bound_counterfactual <- apply(predictions_counterfactual, 2, quantile, 0.975)

# Display predictions with intervals
pred_results_counterfactual <- data.frame(
  new_age = new_age,
  mean_prediction = (mean_predictions_counterfactual),
  lower_bound = (lower_bound_counterfactual),
  upper_bound = (upper_bound_counterfactual)
)

print(pred_results)


## actual scenario
predictions_binary <- ifelse(predictions > 0.10,1,0)
seroprev_by_age <- colMeans(predictions_binary)
titer_by_age <- colMeans(predictions)

df_age_seroprev_stan <- data.table(cbind(c(2:10), seroprev_by_age, titer_by_age))
colnames(df_age_seroprev_stan) <- c('age', 'seroprevalence', 'titer')

## counterfactual scenario
predictions_binary_counterfactual <- ifelse(predictions_counterfactual > 0.10,1,0)
seroprev_by_age <- colMeans(predictions_binary_counterfactual)
titer_by_age <- colMeans(predictions_counterfactual)

df_age_seroprev_stan_counterfactual <- data.table(cbind(c(2:10), seroprev_by_age, titer_by_age))
colnames(df_age_seroprev_stan_counterfactual) <- c('age', 'seroprevalence_counterfactual', 'titer_counterfactual')

## combine
df_age_seroprev_stan_all <- merge(df_age_seroprev_stan, df_age_seroprev_stan_counterfactual, by='age')

df_complete$new_test_d <- ifelse(df_complete$concentration_d >=0.1, 1,0)
df_complete %>% group_by(age) %>% summarise(mean_seroprev = mean(new_test_d), mean_titer=mean(concentration_d), mean_titer_log = exp(mean(log(concentration_d)))) -> test_d

test_d_mean <- merge(df_age_seroprev_stan_all, test_d, by = 'age')
colnames(test_d_mean) <- c('age', 'model_seroprev', 'model_titer','model_seroprev_counterfactual', 'model_titer_counterfactual','observed_seroprev', 'observed_titer', 'observed_titer_log')

(test_d_mean$model_seroprev - test_d_mean$model_seroprev_counterfactual)
(test_d_mean$model_titer - test_d_mean$model_titer_counterfactual)

p <- ggplot() + scale_y_log10() + 
  geom_line(aes(new_age, (mean_predictions_counterfactual)), color = "#D55E00") +
  geom_line(aes(new_age, (mean_predictions)), color = "#0072B2") +
  geom_hline(yintercept=0.1, lty=2) +
  geom_ribbon(aes(x = (new_age), ymin = (lower_bound), ymax = (upper_bound)), alpha = 0.2, fill="#0072B2") +
  geom_ribbon(aes(x = (new_age), ymin = (lower_bound_counterfactual), ymax = (upper_bound_counterfactual)), alpha = 0.2, fill="#D55E00") +
  # labs(title = "Predicted Values with 95% Credible Interval") +
  xlab("Age (years)") + ylab("Diphtheria antitoxin IgG antibody concentration (IU/mL)") + 
  # geom_point(data = test_d, aes(age, observed_titer),  color = "orange", size = 2) +
  geom_point(data = test_d_mean, aes(age, observed_titer_log), size =2) +
  theme_minimal()

png(file = paste0('/continuous_model_additive_stan_log_age_best_negative_paper.png'),
    width = 6,
    height = 4.5,
    units = "in",
    res = 600)
print(p)
dev.off()



## -----------------
## can we get uncertainty intervals around "boosting rates"??

differences <- predictions - predictions_counterfactual

round(quantile(differences[,5], probs=c(0.5, 0.025, 0.975)),3)
round(quantile(differences[,6], probs=c(0.5, 0.025, 0.975)),3)
round(quantile(differences[,7], probs=c(0.5, 0.025, 0.975)),3)
round(quantile(differences[,8], probs=c(0.5, 0.025, 0.975)),3)
round(quantile(differences[,9], probs=c(0.5, 0.025, 0.975)),3)


seroprev_by_age_ci <- round(apply(differences, 2, quantile, probs=c(0.025,0.5, 0.975)),3)
seroprev_by_age_ci
seroprev_by_age_mean <- round(apply(differences, 2, mean),3)
seroprev_by_age_mean

test_df <- data.table(cbind(c(2:10),seroprev_by_age_mean))

annual_boosting <- ggplot() + 
  geom_segment(aes(x = 5, xend=5.9999, y = seroprev_by_age_mean[4], yend= seroprev_by_age_mean[4]), color = '#0072B2') +
  geom_segment(aes(x = 6, xend=6.9999, y = seroprev_by_age_mean[5], yend= seroprev_by_age_mean[5]), color = '#0072B2') +
  geom_segment(aes(x = 7, xend=7.9999, y = seroprev_by_age_mean[6], yend= seroprev_by_age_mean[6]), color = '#0072B2') +
  geom_segment(aes(x = 8, xend=8.9999, y = seroprev_by_age_mean[7], yend= seroprev_by_age_mean[7]), color = '#0072B2') +
  geom_segment(aes(x = 9, xend=9.9999, y = seroprev_by_age_mean[8], yend= seroprev_by_age_mean[8]), color = '#0072B2') +
  geom_segment(aes(x = 10, xend=10.9999, y = seroprev_by_age_mean[9], yend= seroprev_by_age_mean[9]), color = '#0072B2') +
  geom_segment(aes(x = 6, xend=6, y = seroprev_by_age_mean[4], yend= seroprev_by_age_mean[5]), color = '#0072B2') +
  geom_segment(aes(x = 7, xend=7, y = seroprev_by_age_mean[5], yend= seroprev_by_age_mean[6]), color = '#0072B2') +
  geom_segment(aes(x = 8, xend=8, y = seroprev_by_age_mean[6], yend= seroprev_by_age_mean[7]), color = '#0072B2') + 
  geom_segment(aes(x = 9, xend=9, y = seroprev_by_age_mean[7], yend= seroprev_by_age_mean[8]), color = '#0072B2') + 
  geom_segment(aes(x = 10, xend=10, y = seroprev_by_age_mean[8], yend= seroprev_by_age_mean[9]), color = '#0072B2') + 
  annotate("text", x=6.28, y=mean(c(seroprev_by_age_mean[5],seroprev_by_age_mean[4])), label= seroprev_by_age_mean[5] - seroprev_by_age_mean[4], color = 'grey40') + 
  annotate("text", x=7.28, y=mean(c(seroprev_by_age_mean[5],seroprev_by_age_mean[6])), 
           label= (seroprev_by_age_mean[6] - seroprev_by_age_mean[5]- seroprev_by_age_mean[4]), color = 'grey40') + 
  annotate("text", x=8.28, y=mean(c(seroprev_by_age_mean[6],seroprev_by_age_mean[7])), 
           label= (seroprev_by_age_mean[7] - seroprev_by_age_mean[6] ), color = 'grey40') + 
  annotate("text", x=9.28, y=mean(c(seroprev_by_age_mean[7],seroprev_by_age_mean[8])), 
           label= (seroprev_by_age_mean[8] - seroprev_by_age_mean[7] ), color = 'grey40') + 
  annotate("text", x=10.28, y=mean(c(seroprev_by_age_mean[8],seroprev_by_age_mean[9])), 
           label= (seroprev_by_age_mean[9] - seroprev_by_age_mean[8] ), color = 'grey40') + 
  theme_bw() + 
  xlab("Age (years)") + ylab("Cummulative annual boosting rate (IU/mL)")

png(file = paste0('/additive_boosting_rate_titer_best_negative.png'),
    width = 6,
    height = 4,
    units = "in",
    res = 600)
print(annual_boosting)
dev.off()



differences_binary <- predictions_binary - predictions_binary_counterfactual


seroprev_by_age_ci <- apply(differences_binary, 2, quantile, probs=c(0.025,0.5, 0.975))
seroprev_by_age_mean <- round(apply(differences_binary, 2, mean),3)
seroprev_by_age_sd <- apply(differences_binary, 2, sd)

annual_seroconversion <- ggplot() + 
  geom_segment(aes(x = 5, xend=5.9999, y = seroprev_by_age_mean[4], yend= seroprev_by_age_mean[4]), color = '#0072B2') +
  geom_segment(aes(x = 6, xend=6.9999, y = seroprev_by_age_mean[5], yend= seroprev_by_age_mean[5]), color = '#0072B2') +
  geom_segment(aes(x = 7, xend=7.9999, y = seroprev_by_age_mean[6], yend= seroprev_by_age_mean[6]), color = '#0072B2') +
  geom_segment(aes(x = 8, xend=8.9999, y = seroprev_by_age_mean[7], yend= seroprev_by_age_mean[7]), color = '#0072B2') +
  geom_segment(aes(x = 9, xend=9.9999, y = seroprev_by_age_mean[8], yend= seroprev_by_age_mean[8]), color = '#0072B2') +
  geom_segment(aes(x = 10, xend=10.9999, y = seroprev_by_age_mean[9], yend= seroprev_by_age_mean[9]), color = '#0072B2') +
  geom_segment(aes(x = 6, xend=6, y = seroprev_by_age_mean[4], yend= seroprev_by_age_mean[5]), color = '#0072B2') +
  geom_segment(aes(x = 7, xend=7, y = seroprev_by_age_mean[5], yend= seroprev_by_age_mean[6]), color = '#0072B2') +
  geom_segment(aes(x = 8, xend=8, y = seroprev_by_age_mean[6], yend= seroprev_by_age_mean[7]), color = '#0072B2') + 
  geom_segment(aes(x = 9, xend=9, y = seroprev_by_age_mean[7], yend= seroprev_by_age_mean[8]), color = '#0072B2') + 
  geom_segment(aes(x = 10, xend=10, y = seroprev_by_age_mean[8], yend= seroprev_by_age_mean[9]), color = '#0072B2') + 
  annotate("text", x=6.28, y=mean(c(seroprev_by_age_mean[5],seroprev_by_age_mean[4])), label= seroprev_by_age_mean[5] - seroprev_by_age_mean[4], color = 'grey40') + 
  annotate("text", x=7.28, y=mean(c(seroprev_by_age_mean[5],seroprev_by_age_mean[6])), 
           label= (seroprev_by_age_mean[6] - seroprev_by_age_mean[5]- seroprev_by_age_mean[4]), color = 'grey40') + 
  annotate("text", x=8.28, y=mean(c(seroprev_by_age_mean[6],seroprev_by_age_mean[7])), 
           label= (seroprev_by_age_mean[7] - seroprev_by_age_mean[6] ), color = 'grey40') + 
  annotate("text", x=9.28, y=mean(c(seroprev_by_age_mean[7],seroprev_by_age_mean[8])), 
           label= (seroprev_by_age_mean[8] - seroprev_by_age_mean[7] ), color = 'grey40') + 
  annotate("text", x=10.28, y=mean(c(seroprev_by_age_mean[8],seroprev_by_age_mean[9])), 
           label= (seroprev_by_age_mean[9] - seroprev_by_age_mean[8] ), color = 'grey40') + 
  theme_bw() + 
  xlab("Age (years)") + ylab("Cummulative annual seroconversion rate (% points)")

png(file = paste0('/additive_seroconversion_rate_titer_best_negative.png'),
    width = 6,
    height = 4,
    units = "in",
    res = 600)
print(annual_seroconversion)
dev.off()
