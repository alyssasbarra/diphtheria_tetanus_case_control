
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

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

`%!in%` = Negate(`%in%`)


## -----------------------------------------------------------
## Read in final dataset

df_complete <- fread('df_complete.csv')


df_complete$immunity_t_binary <- ifelse(df_complete$concentration_t >=0.1, 1, 0)


df_complete$immunity_d_binary <- ifelse(df_complete$concentration_d >=0.1, 1, 0)


## -----------------------------------------------------------
## Conditional logistic regression set up


#### Run conditional logistic regresion -- tetanus 
datlist <- list(N=nrow(df_complete),
                n_grp=max(df_complete$pair_number),
                n_coef=1,
                x=df_complete[,c("immunity_m_binary")],
                y=df_complete$immunity_t_binary,
                grp=df_complete$pair_number)

clogit_stan2 <- stan("/clogit2.stan", data=datlist, seed=1001)

posterior_samples <- extract(clogit_stan2)
oddsratio_samples_tetanus_sensitivity_analysis_agegroup <- posterior_samples$oddsratio

#### Calculate median and 95% credible intervals for odds ratios
oddsratio_median_tetanus_sensitivity_analysis_agegroup <- apply(oddsratio_samples_tetanus_sensitivity_analysis_agegroup, c(2),FUN= median)
oddsratio_ci_tetanus_sensitivity_analysis_agegroup <- apply(oddsratio_samples_tetanus_sensitivity_analysis_agegroup, c(2), function(x) quantile(x, c(0.5, 0.025, 0.975)))
colnames(oddsratio_ci_tetanus_sensitivity_analysis_agegroup) <- c("measles")

round(t(oddsratio_ci_tetanus_sensitivity_analysis_agegroup),3)


#### Run conditional logistic regresion -- diphtheria
datlist <- list(N=nrow(df_complete),
                n_grp=max(df_complete$pair_number),
                n_coef=1,
                x=df_complete[,c("immunity_m_binary")],
                y=df_complete$immunity_d_binary,
                grp=df_complete$pair_number)

clogit_stan2 <- stan("~/Downloads/clogit2.stan", data=datlist, seed=1001)

posterior_samples <- extract(clogit_stan2)
oddsratio_samples_tetanus_sensitivity_analysis_agegroup <- posterior_samples$oddsratio

#### Calculate median and 95% credible intervals for odds ratios
oddsratio_median_tetanus_sensitivity_analysis_agegroup <- apply(oddsratio_samples_tetanus_sensitivity_analysis_agegroup, c(2),FUN= median)
oddsratio_ci_tetanus_sensitivity_analysis_agegroup <- apply(oddsratio_samples_tetanus_sensitivity_analysis_agegroup, c(2), function(x) quantile(x, c(0.5, 0.025, 0.975)))
colnames(oddsratio_ci_tetanus_sensitivity_analysis_agegroup) <- c("measles")

round(t(oddsratio_ci_tetanus_sensitivity_analysis_agegroup),3)


## -----------------------------------------------------------
## Compare mean titer differences across measles group (by antigen)

#### diphtheria
t.test(concentration_d ~ immunity_m_binary, data = df_complete)
wilcox.test(concentration_d ~ immunity_m_binary, data = df_complete)

var.test(concentration_d ~ immunity_m_binary, data = df_complete)

gg_distribution_d <- ggplot(data = df_complete) + theme_bw() + scale_x_log10() + 
  xlab("Diphtheria antitoxin IgG antibody concentration (IU/mL)") + labs(fill = "Measles serostatus") + 
  geom_density(aes(x=concentration_d, fill=immunity_m_binary, group=immunity_m_binary), alpha = 0.5) + 
  scale_fill_manual(values=c("#0072B2", "#D55E00")) + ggtitle("B. ")

png(file = paste0('/diphtheria_distribution_by_measles.png'),
    width = 8,
    height = 5,
    units = "in",
    res = 600)
print(gg_distribution_d)
dev.off()


ggplot(df_complete) + theme_bw() +
  geom_boxplot(aes(x=immunity_m_binary, y=concentration_t)) + scale_y_log10()


#### tetanus 
t.test(concentration_t ~ immunity_m_binary, data = df_complete)
wilcox.test(concentration_t ~ immunity_m_binary, data = df_complete)


var.test(concentration_t ~ immunity_m_binary, data = df_complete)

gg_distribution_t <- ggplot(data = df_complete) + theme_bw() + scale_x_log10() + xlab("Tetanus antitoxin IgG antibody concentration (IU/mL)") + labs(fill = "Measles serostatus") + 
  geom_density(aes(x=concentration_t, fill=immunity_m_binary, group=immunity_m_binary), alpha = 0.5) + 
  scale_fill_manual(values=c("#0072B2", "#D55E00"))+ ggtitle("A. ")

png(file = paste0('/tetanus_distribution_by_measles.png'),
    width = 8,
    height = 5,
    units = "in",
    res = 600)
print(gg_distribution_t)
dev.off()




