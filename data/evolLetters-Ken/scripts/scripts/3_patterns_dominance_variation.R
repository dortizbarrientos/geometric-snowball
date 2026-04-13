# PCD: Phenotypic consequences of divergence
# Script author: Ken A. Thompson & Avneet K. Chhina
# SCRIPT 3: Dominance analysis
# This script takes derived phenotype data and uses it to evaluate dominance

# PACKAGES s####
library(magicfor)
library(tidyverse) # to do everything

# Helpers
# do a test given data in wide format
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  p = 2*pt(-abs(t),df)
  return(p) 
}

# overwrite select to make sure it's working from dplyr
select <- dplyr::select

# LOAD DATA ####
standardized_trait_df <- read_csv('data/derived_data/standardized_trait_df.csv')
parent_dists_df <- read_csv("data/parent_dists_df.csv")

# SUMMARY STATS ####
pooled_mean_all_pops_traits.df <- standardized_trait_df %>%
  select(Cross.Type, FishID, Cross, category, FW_pop, trait, scaled_val) %>% 
  group_by(category, trait, Cross, Cross.Type, FW_pop) %>% 
  summarise(mean_scaled_val = mean(scaled_val, na.rm = TRUE)) %>% # finding the average, standard deviation, and N (sample size) for traits for P1
  group_by(category, trait, FW_pop) %>%  # now we want to pool across families
  summarise(poolmean = mean(mean_scaled_val))  # adding_P2_n represents the number of individuals for each FW population
    
pooled_sd_n_all_pops_traits.df <- standardized_trait_df %>%
    group_by(category, trait, FW_pop) %>%  # now we want to pool across families
    summarise(poolsd = sd(scaled_val),
              totaln = n())

pooled_mean_sd_n_all_pops_traits.df <- pooled_mean_all_pops_traits.df %>% 
  left_join(., pooled_sd_n_all_pops_traits.df, by = c("category", "trait", "FW_pop"))

# DF of just Little Campbell River
LCR_pooled <- pooled_mean_sd_n_all_pops_traits.df %>% 
  ungroup() %>% 
  filter(category == "P1") %>% 
  select(trait, totaln, poolmean, poolsd) %>% 
  rename(LCR_n = totaln, 
         LCR_mean = poolmean,
         LCR_sd = poolsd) 

# two quick dfs of parent means for later use
LCR_trait_means <- LCR_pooled %>% 
  select(trait, LCR_mean)

FW_trait_means <- pooled_mean_sd_n_all_pops_traits.df %>% 
  filter(category == "P2") %>% 
  ungroup() %>% 
  select(trait, FW_pop, poolmean) %>% 
  rename(FW_mean = poolmean)

# P1: P2 DIVERGENCE TESTS ####
t.test_dataset <- pooled_mean_sd_n_all_pops_traits.df %>% 
  filter(category == "P2") %>% 
  left_join(.,LCR_pooled, by = c("trait")) %>%  # need  to bring the LCR values in line
  rowwise() %>% 
  mutate(t_test_p = t.test2(m1 = poolmean, m2 = LCR_mean, 
                            s1 = poolsd, s2 = LCR_sd,
                            n1 = totaln, n2 = LCR_n),
         diff_stats = t_test_p < 0.05, # a true or false of whether it's diff at P = 0.05
         parent_lcr_sd_diff = abs(poolmean - LCR_mean) / LCR_sd) # divergence in units of LCR SDs

traits_fw_pop_diff <- t.test_dataset %>% 
  ungroup() %>% 
  select(FW_pop, trait, diff_stats) 


scaled_traits_for_dominance.df <- standardized_trait_df %>%
  filter(Parent_Hybrid == "Hybrid") %>% 
  left_join(., traits_fw_pop_diff, by = c("FW_pop", "trait")) %>% 
  filter(diff_stats == T) %>% 
  # take the average within family, then across families (unweighted)
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(scaled_val = mean(scaled_val)) %>% 
  group_by(FW_pop, trait, category) %>% 
  summarise(scaled_val = mean(scaled_val)) %>% 
  left_join(., LCR_trait_means, by = "trait") %>% 
  left_join(., FW_trait_means, by = c("trait", "FW_pop")) %>% 
  mutate(z_dom_scaled_val = (scaled_val - LCR_mean) / (FW_mean - LCR_mean)) %>% 
  left_join(., parent_dists_df, by = "FW_pop")


name_vec <- sort(unique(standardized_trait_df$trait))
code_vec <- c("#AFR", "#DFR",  "GW", 
              "#GR", "HD", "#LAP",
              "BD", "BW",  "FDS", 
              "GRL",   "PF", "PG", 
              "PSL", "SDS", "SNT")
names_codes <- tibble(trait = name_vec, code = code_vec)

trait_mean_scaled_traits_for_dominance.df <- scaled_traits_for_dominance.df %>% 
  group_by(trait, category) %>% 
  summarise(n = n(), 
            sd = sd(z_dom_scaled_val),
            se = sd/sqrt(n),
            z_dom_scaled_val = mean(z_dom_scaled_val)) %>% 
  left_join(., names_codes, "trait")


# mean_dom dataset for Fig. 6B. 
raw_scaled_traits_for_dominance.df_F1 <- scaled_traits_for_dominance.df %>% filter(category == "F1")

scaled_traits_for_dominance.df_F1 <- scaled_traits_for_dominance.df %>% 
  filter(category == "F1") %>% 
  group_by(FW_pop, trait) %>% 
  summarise(mean_dom_scaled_val = mean(z_dom_scaled_val)) %>% 
  group_by(trait) %>% 
  summarise(mean_dom_scaled_val_across_pops = mean(mean_dom_scaled_val),
            n_pops = n(),
            se = sd(mean_dom_scaled_val) / sqrt(n_pops)) %>% 
  left_join(., names_codes , by = "trait")

scaled_traits_for_dominance.df_F2 <- scaled_traits_for_dominance.df %>% 
  filter(category == "F2") %>% 
  group_by(FW_pop, trait) %>% 
  summarise(mean_dom_scaled_val = mean(z_dom_scaled_val)) %>% 
  group_by(trait) %>% 
  summarise(mean_dom_scaled_val_across_pops = mean(mean_dom_scaled_val),
            n_pops = n(),
            se = sd(mean_dom_scaled_val) / sqrt(n_pops)) %>% 
  left_join(., names_codes , by = "trait")

# Dominance statistical analysis DF ####
traits_by_pops_df <- standardized_trait_df %>% 
  select(FW_pop, trait) %>% 
  unique() %>% 
  left_join(., t.test_dataset, by = c("FW_pop", "trait")) %>%
  filter(FW_pop != "LCR",
         diff_stats == T)

# Loop the model across each trait and population
magic_for(progress = T)

for (i in 1:nrow(traits_by_pops_df)) {
  
  df_in_loop = traits_by_pops_df[i,]
  
  dominance_df_in_loop <- standardized_trait_df %>% 
    filter(FW_pop %in% c("LCR", df_in_loop$FW_pop),
           trait == df_in_loop$trait) %>% 
    mutate(dominance = ifelse(category == "F1", 1, 
                              ifelse(category == "F2", 0.5, 0)),
           additive = ifelse(category == "P1", 0, 
                             ifelse(category == "P2", 1, 0.5))) 
  
  dom_mod_in_loop <- lmer(scaled_val ~ dominance + additive + (1|Cross/FW_pop), data = dominance_df_in_loop)
  
  anova_mod <- anova(dom_mod_in_loop, type = 3)
  summary_mod <- summary(dom_mod_in_loop)
  
  dom_p <- anova_mod$`Pr(>F)`[2]
  dom_est <- summary_mod$coefficients[2]
  add_est <- summary_mod$coefficients[3]
  
  FW_pop <- df_in_loop$FW_pop
  trait <- df_in_loop$trait
  print(i)
  put(FW_pop, trait, dom_p, dom_est, add_est)
}

dominance_results <- magic_result_as_dataframe() %>% 
  left_join(., traits_by_pops_df, by = c("FW_pop", "trait")) %>% 
  mutate(fresh_smaller = ifelse(add_est < 0, T, F), 
         dom_smaller = ifelse(dom_est < 0, T, F))

mean(dominance_results$dom_p < 0.05)
sum(dominance_results$dom_p < 0.05)
sum(dominance_results$dom_p < 0.05) / length(dominance_results$dom_p)


sig_dom_direction <- dominance_results %>% 
  filter(dom_p < 0.05) %>% 
  mutate(dom_fresh = ifelse(fresh_smaller == dom_smaller, T, F))
  
mean(sig_dom_direction$dom_fresh == T)

# Analysis of trait variation #### 


# * Sum the variances of individual traits
mean_variances <- filtered_traits_df %>% 
  mutate(scaled_val = as.numeric(scaled_val)) %>% 
  # calculate the variance for each trait within family
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(var = var(scaled_val, na.rm = T),
            family_size = n()) %>% 
  filter(family_size >= 5) %>% 
  # then sum the variances across traits
  group_by(FW_pop, Cross, category) %>% 
  summarise(mean_variance = mean(var, na.rm = T)) %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)) %>% 
  ungroup() %>% 
  filter(category %in% c("F1", "F2"))

individual_trait_variances.df <- filtered_traits_df %>% 
  mutate(scaled_val = as.numeric(scaled_val)) %>% 
  # calculate the variance for each trait within family
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(var = var(scaled_val, na.rm = T),
            family_size = n()) %>% 
  filter(family_size >= 5) %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)) %>% 
  ungroup() %>% 
  filter(category %in% c("F1", "F2"))

population_trait_variances_df <- individual_trait_variances.df %>% 
  group_by(FW_pop, trait, category, parent_trait_divergence) %>% 
  summarise(var = mean(var))





# Save Files ####
write_csv(scaled_traits_for_dominance.df_F1, 'data/derived_data/scaled_traits_for_dominance.df_F1.csv')
write_csv(scaled_traits_for_dominance.df_F2, 'data/derived_data/scaled_traits_for_dominance.df_F2.csv')
