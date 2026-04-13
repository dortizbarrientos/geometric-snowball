# PCD: Phenotypic consequences of divergence
# Script author: Ken A. Thompson & Avneet K. Chhina
# SCRIPT 2: Trait mismatch
# This script takes derived phenotype data and uses it to evaluate trait mismatch

# loading packages####
library(tidyverse) # to do everything

# Load data ####
LCR_scaled_parent_means.df <- read_csv("data/derived_data/LCR_scaled_parent_means.df.csv")
standardized_trait_df <- read_csv('data/derived_data/standardized_trait_df.csv')
category_df <- read_csv('data/derived_data/category_df.csv')
Cross_FW_pop_df <- read_csv('data/derived_data/Cross_FW_pop_df.csv')
dominance_dataset_diff_sd_or_stats <- read_csv('data/derived_data/dominance_dataset_diff_sd_or_stats.csv')
hybrid_scaled_traits_ind.df <- read_csv('data/derived_data/hybrid_scaled_traits_ind.df.csv')


# Helper functions ####
# pw dist function
dist2d <- function(p1_t1, p1_t2, p2_t1, p2_t2, h_t1, h_t2) {
  # define points
  a <- c(h_t1, h_t2) # hybrid trait 1 and hybrid trait 2
  b <- c(p1_t1, p1_t2) # anadromous parent trait 1 and 2
  c <- c(p2_t1, p2_t2) # freshwater parent trait 1 and 2
  
  # calculate the distance between point 'a' and line defined by two points 'b' and 'c'
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
}

# fimctopm to group each trait into a pair with another
make_pairwise_groups <- function(n) {  # number of traits
  combn(1:n, 2) %>% # make pairwise comparisons
    t(.) %>% # transpose the array
    as_tibble() %>% # convert to DF to use with a pipe
    mutate(group = row_number()) %>% # add a column for grouping variable; this will be a concatenated string # so this pairs different traits into two and creates a group. so, for example, 1 and 2 are paired and it is group #1.
    tidyr::gather(c(V1,V2), key = "ColNo", value = "TraitNo") %>%  
    group_by(TraitNo) %>% # group by trait No for concatenation
    mutate(Traits_by_group = paste0(group, collapse = " ")) %>%
    ungroup() %>%
    select(TraitNo, Traits_by_group) %>%
    unique() # removing any duplicates
}

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

# create vectors for naming crosses
f1_CrossName_vec <- c("LCR-LQL", "LCR-PLB", "LCR-PLL", 
                      "LCR-KLN", "LCR-CRN", "LCR-LQB",
                      "LCR-BUL", "LCR-NOR", "LCR-PAQ",
                      "LCR-PCH", "BUL-LCR", "LCR-LQB",
                      "LCR-PXL", "LCR-PXB")

p2_CrossName_vec <- c("LQL-LQL", "LQB-LQB", "PCH-PCH",
                      "PXL-PXL", "PLB-PLB", "PXB-PXB",
                      "KLN-KLN", "PAQ-PAQ", "LQL-LQL", 
                      "BUL-BUL", "PLL-PLL", "NOR-NOR",
                      "CRN-CRN", "PLL-PLL_F2")

f2_CrossName_vec <- c("LCR-LQB_F2", "LCR-CRN_F2", "LCR-KLN_F2",
                      "LCR-PCH_F2", "LCR-PAQ_F2", "LCR-NOR_F2",
                      "LCR-PLB_F2", "LCR-LQL_F2", "LCR-BUL_F2",
                      "LCR-PLL_F2", "LCR-PXB_F2", "LCR-PXL_F2")


#////////////////////#
#### MULTIVARIATE ####
#////////////////////#

# * All traits ####
# *  Hybrids ####
hyb_p1_mv_mm <- hybrid_scaled_traits_ind.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  filter(trait != "Gape.Wt.correct.log") %>% 
  mutate(ha = scaled_hybrid_value - P1) %>% # hybrid value minus anadromous parent value
  select(-scaled_hybrid_value, - P1) %>%  # remove hybrid value and anadromous parent
  mutate(ha = as.numeric(ha))

p1_p2_mv_mm <- FW_scaled_parent_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ba = P2 - P1) %>% # ba = freshwater popuation parent value minus anadrmous parent value
  select(-P2, -P1)

# final step puts it together to calculate mismatch
hybrid_mismatch.df <- hyb_p1_mv_mm %>% 
  left_join(p1_p2_mv_mm, by = c('FW_pop', 'trait')) %>% # joining 'ac' dataset by common columns 'FW_pop' and 'trait'
  group_by(FishID, Cross, FW_pop, category) %>% 
  mutate(t.across = (ha %*% ba) / (ba %*% ba)) %>% # each fish id has same t.across value
  mutate(d = ha - t.across * ba) %>% 
  summarize(hybrid_mismatch = sqrt(sum(d^2))) # dist is mismatch here

# * Parent 'mismatch' ####
parent_mismatch_input_df1 <- fw_parent_dataset_for_p_mismatch %>% 
  left_join(LCR_scaled_parent_means.df, by = "trait") %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ha = scaled_hybrid_value - P1) %>% # hybrid value minus anadromous parent value
  select(-scaled_hybrid_value, - P1) %>%   # remove hybrid value and anadromous parent
  mutate(ha = as.numeric(ha))
  
# calculating multivariate mismatch for hybrids ####
parent_mismatch_input_df2 <- parent_mismatch_input_df1 %>% 
  left_join(p1_p2_mv_mm, by = c('FW_pop', 'trait')) %>% # joining 'ac' dataset by common columns 'FW_pop' and 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>%
  filter(trait != "Gape.Wt.correct.log") %>% 
  group_by(FishID, Cross, FW_pop, category) %>% 
  mutate(t.across = (ha %*% ba) / (ba %*% ba)) %>% # each fish id has same t.across value
  mutate(d = ha - t.across * ba) %>% 
  summarize(dist = sqrt(sum(d^2))) %>% # dist is mismatch here
  rename(parent_mismatch = dist)

# individual-level mismatch
ind_parent_mismatch <- parent_mismatch_input_df2 %>% 
  select(parent_mismatch, FishID, category, FW_pop) %>% 
  rename(mismatch = parent_mismatch) 

# bring it home to make a full DF with parents and hybrids
parent_hybrid_mismatch.df <- hybrid_mismatch.df %>% 
  select(hybrid_mismatch, FishID, category, FW_pop) %>% 
  rename(mismatch = hybrid_mismatch) %>% 
  bind_rows(., ind_parent_mismatch)
  
parent_mismatch_df <- parent_mismatch_input_df2 %>% 
  ungroup() %>% 
  group_by(FW_pop) %>%
  summarize(mean = mean(parent_mismatch)) %>% # finding average mismatch values for F1 and F2 hybrids
  rename(parent_mismatch = mean)

#* dataset with parent mismatch, parent trait divergence and freshwater population name####
mismatch_alltraits_analysis_df <- hybrid_mismatch.df %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  left_join(., parent_mismatch_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence))




# Dominance Effect #### 
# this looks at the MEAN hybrid 
hybrid_scaled_traits_pop_cat_mean.df <- hybrid_scaled_traits_ind.df %>% 
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(mean_pop_cat_value = mean(scaled_hybrid_value)) %>% 
  group_by(FW_pop, trait, category) %>%
  summarise(mean_pop_cat_value = mean(mean_pop_cat_value))

hybrid_mean_mm.df1 <- hybrid_scaled_traits_pop_cat_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  filter(trait != "Gape.Wt.correct.log") %>% 
  mutate(ha = mean_pop_cat_value - P1) %>% # hybrid value minus anadromous parent value
  select(-mean_pop_cat_value, - P1) %>%  # remove hybrid value and anadromous parent
  mutate(ha = as.numeric(ha))

hybrid_mean_mm.df2 <- FW_scaled_parent_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ba = P2 - P1) %>% # ba = freshwater popuation parent value minus anadrmous parent value
  select(-P2, -P1) 

# dataframe of mismatch of the mean hybrid
mean_hybrid_mismatch.df <- hybrid_mean_mm.df1 %>% 
  left_join(p1_p2_mv_mm, by = c('FW_pop', 'trait')) %>% # joining 'ac' dataset by common columns 'FW_pop' and 'trait'
  group_by(FW_pop, category) %>% 
  mutate(t.across = (ha %*% ba) / (ba %*% ba)) %>% # each fish id has same t.across value
  mutate(d = ha - t.across * ba) %>% 
  summarize(mean_pheno_hybrid_mismatch = sqrt(sum(d^2))) %>% # dist is mismatch here
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence))



# EFFECT of trait variation ON mismatch

# segregation variance effect
family_scaled_traits_pop_cat_mean.df <- hybrid_scaled_traits_ind.df %>% 
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(mean_fam_value = mean(scaled_hybrid_value),
            family_size = n()) %>% 
  filter(family_size >= 5) %>% 
  rowwise()

family_mean_mm.df1 <- family_scaled_traits_pop_cat_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  mutate(ha = mean_fam_value - P1) %>% # hybrid value minus anadromous parent value
  select(-mean_fam_value, - P1) %>%  # remove hybrid value and anadromous parent
  mutate(ha = as.numeric(ha))

family_mean_mm.df2 <- FW_scaled_parent_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ba = P2 - P1) %>% # ba = freshwater popuation parent value minus anadrmous parent value
  select(-P2, -P1) 

# dataframe of mismatch of the mean hybrid
mean_family_mismatch.df <- family_mean_mm.df1 %>% 
  left_join(family_mean_mm.df2, by = c('FW_pop', 'trait')) %>% # joining 'ac' dataset by common columns 'FW_pop' and 'trait'
  group_by(FW_pop, Cross, category) %>% 
  mutate(t.across = (ha %*% ba) / (ba %*% ba)) %>% # each fish id has same t.across value
  mutate(d = ha - t.across * ba) %>% 
  summarize(mean_pheno_hybrid_mismatch = sqrt(sum(d^2))) %>% # dist is mismatch here
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence))

family_effect_of_var_on_hybrid_mismatch.df <- hybrid_mismatch.df %>% 
  left_join(., mean_family_mismatch.df, by = c("FW_pop", "Cross", "category")) %>% 
  mutate(var_caused_mismatch = hybrid_mismatch - mean_pheno_hybrid_mismatch) 


# write_files ####
write_csv(parent_dists_df, file  = 'data/derived_data/multivariate_parent_dist.csv')
write_csv(pairwise_mismatch.df, file = 'data/derived_data/pairwise_mismatch.csv')
write_csv(mismatch_alltraits_analysis_df, file = "data/derived_data/mismatch_alltraits_analysis_df.csv")


