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


gg_rep <- ggplot(data = dat.conc.t.final_dup_test) + coord_equal() +
  scale_x_log10() + scale_y_log10() + theme_bw() + xlab("Replicate 1") + ylab("Replicate 2") + 
  geom_vline(xintercept=0.2, color='red', lty=2) + geom_hline(yintercept=0.2, color='red', lty=2) +
  geom_vline(xintercept=0.1, color='blue', lty=2) + geom_hline(yintercept=0.1, color='blue', lty=2) +
  geom_point(aes(x = OD_t, y = OD_t_rep),alpha=0.5) + ggtitle("Tetanus")


png(file = paste0('/tetanus_replicates.png'),
    width = 6,
    height = 4,
    units = "in", 
    res = 600)
print(gg_rep)
dev.off()


mod_rep <- lm(data =dat.conc.t.final_dup_test, OD_t_rep ~ OD_t)

gg_rep <- ggplot(data = dat.conc.t.final_dup_test) + coord_equal() +
  scale_x_log10() + scale_y_log10() + theme_bw() + xlab("Replicate 1") + ylab("Replicate 2") + 
  geom_vline(xintercept=0.2, color='red', lty=2) + geom_hline(yintercept=0.2, color='red', lty=2) +
  geom_vline(xintercept=0.1, color='blue', lty=2) + geom_hline(yintercept=0.1, color='blue', lty=2) +
  geom_point(aes(x = OD_t, y = OD_t_rep),alpha=0.5) + ggtitle("Tetanus") + 
  geom_smooth(aes(x = OD_t, y = OD_t_rep), method = 'lm')



## -----------------------------------------------------------
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

df$immunity_d_binary <- ifelse(df$immunity_d %in% c('high positive', 'positive', 'uncertain'), 'positive', 'negative')

## -----------------------------------------------------------
## Use tetanus thresholds to determine seropositivity  

df$immunity_t_kit <- ifelse(df$concentration_t < 0.01, 'negative',
                            ifelse(df$concentration_t >= 0.01 & df$concentration_t < 0.15, 'border',
                        ifelse(df$concentration_t >= 0.15, 'positive', NA)))

df$immunity_t <- ifelse(df$concentration_t < 0.2, 'negative',
                        ifelse(df$concentration_t >= 0.2, 'positive', NA))

df$immunity_t_binary <- df$immunity_t

df$immunity_t_sensitivity <- ifelse(df$concentration_t < 0.1, 'negative',
                                    ifelse(df$concentration_t >= 0.1, 'positive', NA))

## -----------------------------------------------------------
## Determine how many samples tested so far are in matched pairs 

dat_pairs <- fread('proposed_samples_DipTox_20231103_rearranged-1326samples.csv')
dat_pairs <- subset(dat_pairs, select = c('id', 'pair_number'))

df <- merge(df, dat_pairs, by='id', all.x=T)
df <- subset(df, !is.na(pair_number) & !is.na(concentration_t) & !is.na(concentration_d))

complete_tested_pairs <- df$pair_number[which(duplicated(df$pair_number))]
df_complete <- subset(df, pair_number %in% c(complete_tested_pairs))

df_not_complete <- subset(df, pair_number %!in% c(complete_tested_pairs))


## -----------------------------------------------------------
## Sample demographics

df_complete %>% group_by(immunity_m_binary, age) %>% summarise(n = n()) -> df_age

df_complete %>% group_by(province) %>% summarise(n = n()) -> df_province

df_complete %>% group_by(immunity_m_binary, hivstatusfinal) %>% summarise(n = n()) -> df_hiv

fwrite(df_complete, '~/Dropbox/d&t/analysis_final_dataset/df_complete.csv')



ftable(df_complete$immunity_m_binary ~ df_complete$gender)

ftable(df_complete$immunity_m_binary ~ df_complete$immunity_r)

ftable(df_complete$immunity_m_binary ~ df_complete$immunity_t)

ftable(df_complete$immunity_m_binary ~ df_complete$immunity_d_binary)


