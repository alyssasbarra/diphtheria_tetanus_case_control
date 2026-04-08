
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
library(scales)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

`%!in%` = Negate(`%in%`)

## -----------------------------------------------------------
## Read in final dataset

df_complete <- fread('df_complete.csv')


## -----------------------------------------------------------
## Compare age distribution across concentrations -- diphtheria -- is proportion uncertain protection increasing as age is increasing?

df_complete$immunity_d_uncertain <- ifelse(df_complete$immunity_d == "uncertain", 1, 0)
df_complete$immunity_d_uncertain_n <- ifelse(df_complete$immunity_d == "uncertain", 1, 1)

df_complete$immunity_d_negative <- ifelse(df_complete$immunity_d == "negative", 1, 0)
df_complete$immunity_d_positive <- ifelse(df_complete$immunity_d %in% c("positive", "high positive"), 1, 0)

df_complete_uncertain <- aggregate(df_complete$immunity_d_uncertain, by=list(df_complete$age), FUN='sum')
colnames(df_complete_uncertain) <- c('age', 'uncertain')
df_complete_uncertain_n <- aggregate(df_complete$immunity_d_uncertain_n, by=list(df_complete$age), FUN='sum')
colnames(df_complete_uncertain_n) <- c('age', 'n')
df_complete_negative <- aggregate(df_complete$immunity_d_negative, by=list(df_complete$age), FUN='sum')
colnames(df_complete_negative) <- c('age', 'negative')
df_complete_positive <- aggregate(df_complete$immunity_d_positive, by=list(df_complete$age), FUN='sum')
colnames(df_complete_positive) <- c('age', 'positive')

df_age_uncertain <- merge(df_complete_uncertain, df_complete_uncertain_n, by='age')
df_age_uncertain <- merge(df_age_uncertain, df_complete_negative, by='age')
df_age_uncertain <- merge(df_age_uncertain, df_complete_positive, by='age')

df_age_uncertain$prop_uncertain <- df_age_uncertain$uncertain / df_age_uncertain$n
df_age_uncertain$prop_negative <- df_age_uncertain$negative / df_age_uncertain$n
df_age_uncertain$prop_positive <- df_age_uncertain$positive / df_age_uncertain$n

gg4 <- ggplot(data=df_age_uncertain) + theme_bw() + xlab("Age") + ylab("Proportion") + labs(color = "Serostatus") +
  geom_point(aes(x=age, y=prop_uncertain), color='#56B4E9') + geom_smooth(aes(x=age, y=prop_uncertain), color='#56B4E9') +
  geom_point(aes(x=age, y=prop_negative), color='#D55E00') + geom_smooth(aes(x=age, y=prop_negative), color='#D55E00') + #ggtitle("Diphtheria")+
  geom_point(aes(x=age, y=prop_positive), color='#0072B2') + geom_smooth(aes(x=age, y=prop_positive), color='#0072B2')

png(file = paste0('/diphtheria_serostatus_by_age_poster.png'),
    width = 5,
    height = 2,
    units = "in", 
    res = 600)
print(gg4)
dev.off()




df_age_uncertain$prop_uncertain <- df_age_uncertain$uncertain / df_age_uncertain$n
df_age_uncertain$prop_negative <- df_age_uncertain$negative / df_age_uncertain$n
df_age_uncertain$prop_positive <- df_age_uncertain$positive / df_age_uncertain$n

df_age_uncertain$prop_pos_given_not_neg <- df_age_uncertain$prop_positive / (df_age_uncertain$prop_uncertain + df_age_uncertain$prop_positive)
gg4_paper_d <- ggplot(data = df_age_uncertain) +
  theme_bw() + ylim(0,1) + 
  ggtitle("Diphtheria antitoxin IgG serostatus")+
  xlab("Age (years)") +
  ylab("Proportion") +
  labs(color = "Serostatus") +
  
  geom_point(aes(x = age, y = prop_uncertain, color = "Uncertain")) +
  geom_smooth(aes(x = age, y = prop_uncertain, color = "Uncertain")) +
  
  geom_point(aes(x = age, y = prop_negative, color = "Negative")) +
  geom_smooth(aes(x = age, y = prop_negative, color = "Negative")) +
  
  geom_point(aes(x = age, y = prop_positive, color = "Positive")) +
  geom_smooth(aes(x = age, y = prop_positive, color = "Positive")) +
  
  scale_color_manual(
    values = c(
      "Uncertain" = "#56B4E9",
      "Negative"  = "#D55E00",
      "Positive"  = "#0072B2"
    )
  )

gg4_paper_d_pos <- ggplot(data = df_age_uncertain) +
  theme_bw() + ylim(0,1) + 
  ggtitle("Proportion diphtheria seropositive among those not seronegative")+
  xlab("Age (years)") +
  ylab("Proportion") +
  labs(color = "") +
  
  geom_point(aes(x = age, y = prop_pos_given_not_neg, color='Proportion')) +
  geom_smooth(aes(x = age, y = prop_pos_given_not_neg, color='Proportion')) +
  
  scale_color_manual(
    values = c(
      "Proportion" = "darkgreen"
    )
  )

#### Test regressions
summary(lm(data = df_age_uncertain, prop_negative ~ age))
summary(lm(data = df_age_uncertain, prop_uncertain ~ age))
summary(lm(data = df_age_uncertain, prop_positive ~ age))



## -----------------------------------------------------------
## Compare age distribution across concentrations -- tetanus -- is proportion uncertain protection increasing as age is increasing?

df_complete$immunity_t_negative <- ifelse(df_complete$immunity_t == "negative", 1, 0)
df_complete$immunity_t_positive <- ifelse(df_complete$immunity_t %in% c("positive"), 1, 0)

df_complete$immunity_t_positive_n <- ifelse(df_complete$immunity_t == "positive", 1, 1)

df_complete_positive_n <- aggregate(df_complete$immunity_t_positive_n, by=list(df_complete$age), FUN='sum')
colnames(df_complete_positive_n) <- c('age', 'n')
df_complete_negative <- aggregate(df_complete$immunity_t_negative, by=list(df_complete$age), FUN='sum')
colnames(df_complete_negative) <- c('age', 'negative')
df_complete_positive <- aggregate(df_complete$immunity_t_positive, by=list(df_complete$age), FUN='sum')
colnames(df_complete_positive) <- c('age', 'positive')

df_age_uncertain <- merge(df_complete_positive_n, df_complete_negative, by='age')
df_age_uncertain <- merge(df_age_uncertain, df_complete_positive, by='age')

df_age_uncertain$prop_negative <- df_age_uncertain$negative / df_age_uncertain$n
df_age_uncertain$prop_positive <- df_age_uncertain$positive / df_age_uncertain$n


gg4_paper_t <- ggplot(data = df_age_uncertain) +
  theme_bw() +
  ggtitle("Tetanus antitoxin IgG serostatus")+
  xlab("Age (years)") +
  ylab("Proportion") +
  labs(color = "Serostatus") +
  
  ylim(0,1) +
  
  geom_point(aes(x = age, y = prop_negative, color = "Negative")) +
  geom_smooth(aes(x = age, y = prop_negative, color = "Negative")) +
  
  geom_point(aes(x = age, y = prop_positive, color = "Positive")) +
  geom_smooth(aes(x = age, y = prop_positive, color = "Positive")) +
  
  scale_color_manual(
    values = c(
      "Uncertain" = "#56B4E9",
      "Negative"  = "#D55E00",
      "Positive"  = "#0072B2"
    )
  )


library(cowplot)
gg4_paper_t_noleg <- gg4_paper_t + theme(legend.position = "none")

combined_plot <- plot_grid(
  gg4_paper_t,
  gg4_paper_d,
  gg4_paper_d_pos,
  ncol = 1,
  labels = c("a.", "b.", "c."),
  label_fontface = "bold",
  align = "v"
)


png(file = paste0('/diphtheria_tetanus_serostatus_by_age_paper.png'),
    width = 6,
    height = 6,
    units = "in", 
    res = 600)
print(combined_plot)
dev.off()




#### Test regressions
summary(lm(data = df_age_uncertain, prop_negative ~ age))
summary(lm(data = df_age_uncertain, prop_positive ~ age))

gg5_ta <- ggplot(data = df_complete, aes(x=age, y=concentration_t)) + geom_boxplot(aes(x=age, y=concentration_t, group=age)) + scale_y_log10()+ 
  theme_bw() + #ggtitle("Tetanus") + 
  scale_x_continuous( breaks=pretty_breaks()) +
  xlab("Age") + ylab('IU/mL') + geom_hline(aes(yintercept=0.2), color='red', lty=2, lwd=1)


png(file = paste0('/tetanus_conc_by_age_epi_slides.png'),
    width = 5,
    height = 2,
    units = "in", 
    res = 600)
print(gg5_ta)
dev.off()




gg5_ta <- ggplot(data = df_complete, aes(x=age, y=concentration_t)) + geom_boxplot(aes(x=age, y=concentration_t, group=age)) + scale_y_log10()+ 
  theme_bw() + ggtitle("Tetanus") + xlab("Age") + ylab('IU/mL') + geom_hline(aes(yintercept=0.2), color='red', lty=2) + facet_wrap(~immunity_m_binary)


png(file = paste0('/tetanus_conc_by_age_by_measles.png'),
    width = 14,
    height = 6,
    units = "in", 
    res = 600)
print(gg5_ta)
dev.off()


## -----------------------------------------------------------
## Regress diphtheria titer by age

summary(glm(data = df_complete, log(concentration_t) ~ age))
summary(glm(data = df_complete, log(concentration_t) ~ age + immunity_m_binary))

summary(glm(data = df_complete, log(concentration_d) ~ age))
summary(glm(data = df_complete, log(concentration_d) ~ age + immunity_m_binary))

gg_t_age <- ggplot(data = df_complete) + facet_wrap(~age, ncol=1)+ scale_x_log10()+ 
  geom_density(aes(x=concentration_t, group = immunity_m_binary, fill = immunity_m_binary), alpha = 0.5)+
  geom_vline(aes(xintercept = 0.2), lty=2) + labs(fill = 'Measles serostatus') + xlab("Tetanus IU/mL")


gg_d_age <- ggplot(data = df_complete) + facet_wrap(~age, ncol=1)+ scale_x_log10()+ 
  geom_density(aes(x=concentration_d, group = immunity_m_binary, fill = immunity_m_binary), alpha = 0.5)+
  geom_vline(aes(xintercept = 0.2), lty=2) + labs(fill = 'Measles serostatus')+ xlab("Diphtheria IU/mL")

png(file = paste0('/tetanus_by_age_measles.png'),
    width = 4,
    height = 10,
    units = "in",
    res = 600)
print(gg_t_age)
dev.off()

png(file = paste0('/diphtheria_by_age_measles.png'),
    width = 4,
    height = 10,
    units = "in",
    res = 600)
print(gg_d_age)
dev.off()


## -----------------------------------------------------------    
## diphtheria and tetanus correlation by age

gg_corr <- ggplot(data=df_complete) + #subset(df_complete, age %in% c(2,5,10))) + 
  geom_point(aes(x=concentration_d, y=concentration_t,color=as.factor(age)),alpha=0.2) + geom_vline(xintercept = 0.1, lty=2)+
  geom_hline(yintercept = 0.2, lty=2)+ guides(color = 'none') + 
  scale_x_log10() + scale_y_log10() + theme_bw() + labs(color='Age') + xlab("Diphtheria IU/mL") + ylab("Tetanus IU/mL") +
  facet_wrap(~age) + geom_smooth(aes(x=concentration_d, y=concentration_t), color='black', method = 'lm') + 
  scale_color_manual(values=c("#0072B2","#0072B2","#0072B2", "#D55E00",  "#D55E00",  "#D55E00", "#0072B2",  "#D55E00",  "#D55E00"))


png(file = paste0('/diphtheria_tetanus_corr_by_age.png'),
    width = 6,
    height = 4,
    units = "in",
    res = 600)
print(gg_corr)
dev.off()


df_complete$age_paper <- paste0(df_complete$age, '-year-olds')
df_complete$age_paper <- factor(df_complete$age_paper , levels=c("2-year-olds", "3-year-olds", "4-year-olds", "5-year-olds",
                                                                 "6-year-olds", "7-year-olds", "8-year-olds", "9-year-olds", "10-year-olds"))
gg_corr <- ggplot(data=df_complete) + #subset(df_complete, age %in% c(2,5,10))) + 
  geom_point(aes(x=concentration_d, y=concentration_t,color=as.factor(age)),alpha=0.2) + geom_vline(xintercept = 0.1, lty=2)+
  geom_hline(yintercept = 0.2, lty=2)+ guides(color = 'none') + 
  scale_x_log10() + scale_y_log10() + theme_bw() + labs(color='Age') + xlab("Diphtheria antitoxin IgG antibody concentration (IU/mL)") + ylab("Tetanus antitoxin IgG antibody concentration (IU/mL)") +
  facet_wrap(~age_paper) + geom_smooth(aes(x=concentration_d, y=concentration_t), color='black', method = 'lm') + 
  scale_color_manual(values=c("#0072B2","#0072B2","#0072B2", "#D55E00",  "#D55E00",  "#D55E00", "#0072B2",  "#D55E00",  "#D55E00"))


png(file = paste0('/diphtheria_tetanus_corr_by_age_paper.png'),
    width = 6,
    height = 4,
    units = "in",
    res = 600)
print(gg_corr)
dev.off()


cor.test(subset(df_complete, age ==2)$concentration_d, subset(df_complete, age ==2)$concentration_t)
cor.test(subset(df_complete, age ==3)$concentration_d, subset(df_complete, age ==3)$concentration_t)
cor.test(subset(df_complete, age ==4)$concentration_d, subset(df_complete, age ==4)$concentration_t)
cor.test(subset(df_complete, age ==5)$concentration_d, subset(df_complete, age ==5)$concentration_t)
cor.test(subset(df_complete, age ==6)$concentration_d, subset(df_complete, age ==6)$concentration_t)
cor.test(subset(df_complete, age ==7)$concentration_d, subset(df_complete, age ==7)$concentration_t)
cor.test(subset(df_complete, age ==8)$concentration_d, subset(df_complete, age ==8)$concentration_t)
cor.test(subset(df_complete, age ==9)$concentration_d, subset(df_complete, age ==9)$concentration_t)
cor.test(subset(df_complete, age ==10)$concentration_d, subset(df_complete, age ==10)$concentration_t)



