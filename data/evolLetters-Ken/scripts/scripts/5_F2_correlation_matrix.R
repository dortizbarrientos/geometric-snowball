# Trait correlation matrix

# load packages ####
library(correlation)
library(tidyverse)

# load data ####
mouth_correction_data_log <- read_csv('data/derived_data/size_mouth_corrected_pheno_data.csv')

# the correlation matrix is meant to be something similar to a GENETIC CORRELATION matrix. We measure correlations within F2 families.
F2_correlation_matrix <- mouth_correction_data_log %>% 
  filter(category == "F2") %>% # only retain F2 hybrids
  dplyr::select(Cross, N.plates:GillRakerNo, size_correction_Body.Dpt.log.res:size_correction_Second.Dor.Spine.Lt.log.res, size_correction_Body.Wt.log.res:head.lt.correct.log) %>% 
  group_by(Cross) %>% 
  na.omit() %>% 
  mutate(n_in_fam = n()) %>% # how many individuals in each family
  filter(n_in_fam > 9) %>%  # retain famiiles with 10 or more individuals
  # filter(Cross %in% c("LCR-BUL_F2.1", "LCR-PXB_F2.1")) %>%  # checking to see if it works for just a couple exmaple families
  # dplyr::select(Cross, N.plates, Dorsal.Ray, GillRakerNo) %>%  # checking to see if it works for just a couple example traits
  mutate(group_id = group_indices()) %>% 
  group_by(group_id) %>% 
  dplyr::select(-Cross, -n_in_fam) %>% 
  correlation(p_adjust = "none")

correlation_coeffs_r <- F2_correlation_matrix$r
correlation_coeffs_p <- F2_correlation_matrix$p
correlation_coeffs_name1 <- F2_correlation_matrix$Parameter1
correlation_coeffs_name2 <- F2_correlation_matrix$Parameter2
correlation_tibble <- tibble(correlation_coeffs_r, correlation_coeffs_p, correlation_coeffs_name1, correlation_coeffs_name2) %>% 
  mutate(sig = ifelse(correlation_coeffs_p < 0.05, T, F))

# summary stats for reporting in manuscript ####
# mean_r
median(abs(correlation_coeffs_r))
# fraction significant
1- sum(correlation_coeffs_p < 0.05)/length(correlation_coeffs_p)

# write file
write_csv(correlation_tibble, 'data/derived_data/correlation_matrix_summaries.csv')

