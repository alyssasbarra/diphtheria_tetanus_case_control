
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


## -----------------------------------------------------------
## Preliminary plots to visualize d & t

gg1 <- ggplot() + geom_point(data=df_complete, aes(x=(concentration_d), y = (concentration_t), color=immunity_m_binary), alpha = 0.2) + 
  theme_bw() + xlab("Diphtheria IU/mL") + ylab("Tetanus IU/mL") + labs(color="Measles serostatus") +
  geom_vline(xintercept = (0.01), lty = 2) + scale_y_log10() + scale_x_log10()+
  geom_vline(xintercept = (0.1), lty = 2) + 
  geom_vline(xintercept = (1), lty = 2) + 
  geom_hline(yintercept = (0.2), lty = 2) 

gg2_d <- ggplot() + geom_histogram(data=df_complete, aes(x=concentration_d), fill='darkgrey') + theme_bw() + 
  geom_vline(xintercept=0.01, color='#0072B2', lty=2) + 
  geom_vline(xintercept=0.1, color='#0072B2', lty=2) + scale_x_log10()+
  geom_vline(xintercept=1, color='#0072B2', lty=2) + 
  ggtitle("Diphtheria") + xlab("IU/mL") + ylab("Frequency")

gg2_t <- ggplot() + geom_histogram(data=df_complete, aes(x=concentration_t), fill='darkgrey') + theme_bw() + scale_x_log10() +
  geom_vline(xintercept=0.2, color='#0072B2', lty=2) + ggtitle("Tetanus") + xlab("IU/mL") + ylab("Frequency")


png(file = paste0('~/Dropbox/d&t/analysis_final_dataset/diphtheria_titers.png'),
    width = 8,
    height = 4,
    units = "in", 
    res = 600)
print(gg2_d)
dev.off()

png(file = paste0('~/Dropbox/d&t/analysis_final_dataset/tetanus_titers.png'),
    width = 8,
    height = 4,
    units = "in", 
    res = 600)
print(gg2_t)
dev.off()



## -----------------------------------------------------------
## Overall crude 2x2 tables

#### diptheria 2xmultiple
df_overall <- subset(df_complete, select= c("immunity_m_binary", "immunity_t_binary", "immunity_d"))

df_overall_d = df_overall %>% 
  group_by(immunity_m_binary, immunity_d) %>% 
  count() %>% as.tibble

overall_table_d <- xtabs(n ~ immunity_m_binary + immunity_d, data = df_overall_d)
print(overall_table_d)

#### diptheria 2x2
df_overall <- subset(df_complete, select= c("immunity_m_binary", "immunity_t_binary", "immunity_d_binary"))

df_overall_d = df_overall %>% 
  group_by(immunity_m_binary, immunity_d_binary) %>% 
  count() %>% as.tibble

overall_table_d <- xtabs(n ~ immunity_m_binary + immunity_d_binary, data = df_overall_d)
print(overall_table_d)


#### tetanus 2xkit cuts
df_overall <- subset(df_complete, select= c("immunity_m_binary", "immunity_t_kit", "immunity_d"))

df_overall_t = df_overall %>% 
  group_by(immunity_m_binary, immunity_t_kit) %>% 
  count() %>% as.tibble

overall_table_t <- xtabs(n ~ immunity_m_binary + immunity_t_kit, data = df_overall_t)
print(overall_table_t)

df_complete$immunity_d_binary2 <- ifelse(df_complete$immunity_d_binary == 'positive', 1, 0)
df_complete$immunity_t_binary2 <- ifelse(df_complete$immunity_t_binary == 'positive', 1, 0)


#### tetanus 2x2
df_overall <- subset(df_complete, select= c("immunity_m_binary", "immunity_t", "immunity_d"))

df_overall_t = df_overall %>% 
  group_by(immunity_m_binary, immunity_t) %>% 
  count() %>% as.tibble

overall_table_t <- xtabs(n ~ immunity_m_binary + immunity_t, data = df_overall_t)
print(overall_table_t)

df_complete$immunity_d_binary2 <- ifelse(df_complete$immunity_d_binary == 'positive', 1, 0)
df_complete$immunity_t_binary2 <- ifelse(df_complete$immunity_t_binary == 'positive', 1, 0)


## -----------------------------------------------------------
## Diphtheria and tetanus correlation

#### Test correlation between diphtheria and tetanus
cor.test(df_complete$concentration_t, df_complete$concentration_d)
ftable(df_complete$immunity_d ~ df_complete$immunity_t)

#### who are diphtheria positive and tetanus negative?
t_neg_d_pos <- subset(df_complete, immunity_t == 'negative'  & immunity_d %in% c('positive', 'high positive'))


## -----------------------------------------------------------
## Conditional logistic regression set up


#### Run conditional logistic regresion -- tetanus 
datlist <- list(N=nrow(df_complete),
                n_grp=max(df_complete$pair_number),
                n_coef=1,
                x=df_complete[,c("immunity_m_binary")],
                y=df_complete$immunity_t_binary2,
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
                y=df_complete$immunity_d_binary2,
                grp=df_complete$pair_number)

clogit_stan2 <- stan("/clogit2.stan", data=datlist, seed=1001)

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
  geom_density(aes(x=(concentration_d), fill=immunity_m_binary, group=immunity_m_binary), alpha = 0.5) + 
  scale_fill_manual(values=c("#0072B2", "#D55E00")) + ggtitle("B. ") + geom_vline(aes(xintercept = (0.01)), lty= 2)  + geom_vline(aes(xintercept = 0.1), lty= 2)

gg_distribution_d <- ggplot(data = df_complete) +
  theme_bw() +
  scale_x_log10(limits = c(0.003, NA)) +   # extend lower limit
  xlab("Diphtheria antitoxin IgG antibody concentration (IU/mL)") +
  labs(fill = "Measles serostatus") +
  geom_density(aes(x = concentration_d,
                   fill = immunity_m_binary,
                   group = immunity_m_binary),
               alpha = 0.5) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  ggtitle("B.") +
  geom_vline(xintercept = 0.01, lty = 2) +
  geom_vline(xintercept = 0.1, lty = 2)


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
  scale_fill_manual(values=c("#0072B2", "#D55E00"))+ ggtitle("A. ") + geom_vline(aes(xintercept = 0.2), lty= 2)

png(file = paste0('/tetanus_distribution_by_measles.png'),
    width = 8,
    height = 5,
    units = "in",
    res = 600)
print(gg_distribution_t)
dev.off()




