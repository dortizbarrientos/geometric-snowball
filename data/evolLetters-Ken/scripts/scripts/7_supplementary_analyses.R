# Supplemental analyses

#///////////////////////////////#
#### PAIRWISE TRAIT MISMATCH ####
#///////////////////////////////#

# create smaller df of divergent traits
div_trait_data_smaller <- dominance_dataset_diff_sd_or_stats %>% 
  # make a new col that is T if either stats *or (|)* SD is true
  mutate(diff = ifelse(Different.or.not.stats == T | Different_or_not_sd ==T, T , F)) %>% 
  select(FW_pop, traits, diff) %>% 
  # rename traits as trait
  rename(trait = traits)

# dataset where we filter only to hybrids
filtered_traits_df <- standardized_trait_df %>% 
  # bring in div TF data
  left_join(., div_trait_data_smaller, by = c("FW_pop", "trait")) %>% 
  # filter only to hybrids (# if we want to filter only divergent traits, please add "diff == T")
  filter(category %in% c("F1", "F2", "P2"),
         # , diff == T,
         trait != "size_correction_Eye.Diam.log.res") %>% #%>% # we want to calculate mismatch for FW parents AND hybrids.
  mutate(scaled_val = as.numeric(scaled_val)) 

# generate a df with mean trait values for each pop/cat
filtered_traits_mean_hyb.df <- filtered_traits_df %>% 
  group_by(FW_pop, Cross, category, trait) %>% 
  summarise(scaled_val = mean(scaled_val),
            family_size = n()) %>% 
  group_by(FW_pop, category, trait) %>% 
  summarise(scaled_val = mean(scaled_val)) %>% 
  mutate(FishID = paste("mean", FW_pop, category, sep = "_"))

# now take traits and make them pairwise groups
## first step is to determine how many trait pairs can be made for each individual
pairwise_trait_df <- filtered_traits_df %>% 
  group_by(FW_pop, Cross, category, FishID) %>% 
  bind_rows(., filtered_traits_mean_hyb.df) %>% # this allows us to calculate the MEDIAN FISH's pairwise mismatch
  select(trait) %>% 
  unique() %>% 
  summarise(number_of_traits = as.character(length(trait))) %>% 
  filter(number_of_traits > 1) %>% 
  # not totally sure why but filtering is dropping FishID as a grouping var, so adding it back in
  group_by(FW_pop, Cross, category, FishID) %>% 
  # uncomment below to practice with just one fish
  # filter(FishID == "KT_Gac_2426") %>%
  # uncomment below to practice with just two fish
  # filter(FishID %in% c("KT_Gac_2426", "KT_Gac_108")) %>%
  do(make_pairwise_groups(n = .$number_of_traits)) 

# make the DF really long; one row per trait per pair per individual
pairwise_trait_df_long <- filtered_traits_df %>%
  bind_rows(., filtered_traits_mean_hyb.df) %>% # this allows us to calculate the MEDIAN FISH's pairwise mismatch
  # uncomment below to practice with just one fish
  group_by(FW_pop, Cross, category, FishID) %>% # some grouping here is redundant
  mutate(TraitNo = seq_along(trait)) %>% 
  group_by(FishID) %>%  
  select(FishID, trait, TraitNo, scaled_val) %>% 
  left_join(., pairwise_trait_df, by = c("FishID", "TraitNo")) %>% # join in the traits_by_group col
  separate_rows(Traits_by_group, sep = " ") %>% # 
  group_by(FishID, Traits_by_group) %>% 
  # Give each trait in a group a number, 1 or 2
  mutate(TraitNo_in_group = row_number()) 

# make a df of PARENT trait values that can be merged with hybrids
fw_parent_filtered_traits_df <- standardized_trait_df %>% 
  # bring in TF data
  left_join(., div_trait_data_smaller, by = c("FW_pop", "trait")) %>% 
  # filter also only to hybrids (and if we want to filter to divergent traits, below write:diff == T)
  filter(category %in% c("P2")) %>% 
  # select(FW_pop:value) %>%  
  group_by(FW_pop, Cross, trait) %>% 
  summarise(mean = mean(scaled_val)) %>% #average value by grouping by freshwater population, cross, trait
  group_by(FW_pop,trait) %>% 
  summarise(mean = mean(mean)) # average value by grouping by freshwater population and trait

lcr_trait_means <- standardized_trait_df %>% 
  filter(category == c("P1")) %>% # only include anadromous parent
  # select(FW_pop:value) %>% 
  group_by(trait, Cross) %>%  
  summarise(mean = mean(scaled_val)) %>% # calculate mean by grouping by trait and cross
  group_by(trait) %>% 
  summarise(mean = mean(mean)) # calculate mean by grouping by trait

# marine traits df
pairwise_trait_means_p1 <- pairwise_trait_df_long %>% 
  select(FishID, FW_pop, trait, TraitNo, TraitNo_in_group, Traits_by_group) %>% 
  # bring in the parent data
  left_join(., lcr_trait_means, by = c("trait")) %>% 
  select(-trait, -TraitNo) %>%
  ungroup() %>% 
  spread(key = TraitNo_in_group, value = mean) %>% # spread the mean values
  rename(p1_t1 = `1`, p1_t2 = `2`)  

# create a DF containing FW trait means that can be matched w/ hybrid DF
pairwise_trait_means_p2 <- pairwise_trait_df_long %>% 
  select(FishID, FW_pop, trait, TraitNo, TraitNo_in_group, Traits_by_group) %>% 
  # bring in the parent data
  left_join(., fw_parent_filtered_traits_df, by = c("FW_pop", "trait")) %>% #freshwater parent data
  select(-trait, -TraitNo) %>% 
  ungroup() %>% 
  spread(key = TraitNo_in_group, value = mean) %>% # spread mean values
  rename(p2_t1 = `1`, p2_t2 = `2`) 

# for now just filter to a single trait pair to keep the computing easy
# each trait pair should have two rows
pairwise_trait_df_simple <- pairwise_trait_df_long %>% 
  select(FishID, FW_pop, category, Traits_by_group, TraitNo_in_group, scaled_val, Cross) %>% 
  # spread into a common format that can be merged later
  spread(key = TraitNo_in_group, value = scaled_val) %>% 
  rename(h_t1 = `1`, h_t2 = `2`) %>% 
  # add the marine parent means
  # and the LCR parent means
  left_join(., pairwise_trait_means_p1, by = c("FishID", "FW_pop", "Traits_by_group")) %>% # join anadromous parent data
  # add the FW parent means
  left_join(., pairwise_trait_means_p2, by = c("FishID", "FW_pop", "Traits_by_group")) # join freshwater parent data

# DF with pairwise trait info for each fish
pairwise_trait_name_info <- pairwise_trait_df_long %>% 
  ungroup() %>% 
  select(FishID, trait, Traits_by_group) %>% 
  mutate(Traits_by_group = as.numeric(Traits_by_group)) %>% 
  arrange(FishID, Traits_by_group, trait) %>% 
  group_by(FishID, Traits_by_group) %>% 
  mutate(TraitOrder = row_number()) %>%
  spread(key = TraitOrder, value = trait) %>% 
  mutate(Group_w_TraitNames = paste(`1`, `2`, sep = "_")) %>% 
  select(FishID, Traits_by_group, Group_w_TraitNames)

# simple df of fish ID and category
fid_cat_df <- standardized_trait_df %>% 
  ungroup() %>% 
  select(FishID, category) %>% # only include category and FishID
  filter(category %in% c("P2", "F1", "F2")) # only include hybrids within category column

# * calculate 'mismatch' for each trait pair ####
pairwise_mismatch.df <- pairwise_trait_df_simple %>% 
  drop_na(h_t1:p2_t2) %>% 
  mutate(mismatch = dist2d(p1_t1 = p1_t1, p1_t2 = p1_t2, 
                           p2_t1 = p2_t1, p2_t2 = p2_t2, 
                           h_t1 = h_t1, h_t2 = h_t2)) %>% 
  mutate(parent_dist_focal_pair = sqrt((p1_t1 - p2_t1)^2 + (p1_t2 - p2_t2)^2)) %>% 
  ungroup() %>% 
  mutate(Traits_by_group = as.numeric(Traits_by_group)) %>% 
  left_join(., pairwise_trait_name_info, by = c("FishID", "Traits_by_group")) %>% # include 'pairwise_trait_name_info' data
  select(FishID, FW_pop, mismatch, parent_dist_focal_pair, Group_w_TraitNames, Cross) %>% 
  left_join(parent_dists_df, by = "FW_pop") %>% # join parent trait divergence data
  left_join(., fid_cat_df, by = "FishID") %>%  # join category, fish ID data
  distinct()  # removes duplicates

# pairwise_mismatch_summaries
pw_mismatch_summaries_all_inds <- pairwise_mismatch.df %>% 
  na.omit() %>% 
  # mutate(mismatch = log(mismatch)) %>% # toggle this for log
  group_by(FW_pop, Cross, category, Group_w_TraitNames) %>%
  summarise(family_size = n(),
            family_mean_mismatch = mean(mismatch),
            sd_mismatch = sd(mismatch)) %>% 
  group_by(FW_pop, category, Group_w_TraitNames) %>% 
  mutate(n_minus1_sd2 = (family_size - 1) * sd_mismatch^2,
         n_fams_sd = ifelse(is.na(sd_mismatch) == T, F, T),
         n_groups_sd = sum(n_fams_sd),
         fam_size_sd = ifelse(family_size > 1, family_size, NA)) %>% 
  summarise(mean_mismatch = weighted.mean(family_mean_mismatch, family_size),
            sample_size = sum(family_size),
            sd_mismatch = sqrt(sum(n_minus1_sd2, na.rm = T) / (sum(fam_size_sd, na.rm=T) - n_groups_sd))) %>% 
  unique()

# pool 'mismatch' across all parents
pw_mismatch_summaries_Ps <- pw_mismatch_summaries_all_inds %>% 
  filter(category == "P2") %>% 
  ungroup() %>% 
  rename(p_mismatch = mean_mismatch,
         p_n = sample_size, 
         p_sd = sd_mismatch) %>% 
  select(-category) %>% 
  group_by(Group_w_TraitNames) %>% 
  mutate(n_minus1_sd2 = (p_n - 1) * p_sd^2,
         n_groups = n()) %>% 
  summarise(pooled_p_mismatch = weighted.mean(p_mismatch, p_n),
            pooled_p_n = sum(p_n),
            pooled_p_sd = sqrt(sum(n_minus1_sd2) / (sum(p_n) - n_groups))) %>% 
  unique() 



# 
pw_mismatch_hyb_tests <- pw_mismatch_summaries_all_inds %>% 
  # filter(category %in% c("F1", "F2")) %>% 
  rename(h_mismatch = mean_mismatch,
         h_n = sample_size,
         h_sd = sd_mismatch) %>% 
  left_join(., pw_mismatch_summaries_Ps, by = c("Group_w_TraitNames")) %>% 
  rowwise() %>% 
  mutate(p = t.test2(m1 = h_mismatch, m2 = pooled_p_mismatch, 
                     s1 = h_sd ,s2 = pooled_p_sd,
                     n1 = h_n, n2 = pooled_p_n,
                     m0=0,equal.variance=FALSE)) %>% 
  mutate(diff = ifelse(p < (0.05 / nrow(.)), T, F),
         h_bigger = ifelse(h_mismatch > pooled_p_mismatch, T, F),
         diff_and_bigger = ifelse(diff == T & h_bigger == T, T, F))


# Snowball ####
origin_cat <- c(0,0)
F1_F2 <- c("F1", "F2")

# simple df containing the 'origin' to fit intercept through 0
origin_tib <- tibble(category = F1_F2, parent_trait_divergence = origin_cat, num_diff_and_bigger = origin_cat)

# how many trait pairs are more mismatched than expected?
pw_mismatch_diff_counts <- pw_mismatch_hyb_tests %>% 
  group_by(FW_pop, category) %>% 
  summarise(num_diff_and_bigger = sum(diff_and_bigger),
            num_bigger = sum(h_bigger)) %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)) %>% 
  bind_rows(origin_tib)


# Variance Effect ####
# subtract the mean effects!
# the mean hybrid is the one with 'na' in category
pairwise_mismatch.df_median <- pairwise_mismatch.df %>% 
  filter(is.na(category)) %>%
  select(-FW_pop, -category, -Cross, -parent_trait_divergence) %>%
  separate(FishID, into = c("type", "FW_pop", "category")) %>% 
  select(-type, -parent_dist_focal_pair) %>% 
  rename(median_hyb_mismatch = mismatch)

variance_effect_pairwise <- pairwise_mismatch.df %>% 
  filter(category %in% c("F1", "F2")) %>% 
  left_join(., pairwise_mismatch.df_median, by = c("FW_pop", "category", "Group_w_TraitNames")) %>% 
  mutate(variance_effect = mismatch - median_hyb_mismatch)


# Plate evolution effect ####
plate_means.df <- standardized_trait_df %>% 
  filter(trait == "N.plates",
         category %in% c("P1", "P2")) %>% 
  group_by(FW_pop, Cross, category) %>% 
  summarise(mean_pop_cat_value = mean(scaled_val),
            family_size = n()) %>% 
  group_by(FW_pop, category) %>%
  summarise(mean_pop_cat_value = weighted.mean(mean_pop_cat_value, family_size))

LCR_plate_scaled_mean <- plate_means.df %>% filter(FW_pop == "LCR")

parent_plate_midpt.df <- plate_means.df %>% 
  filter(category == "P2") %>% 
  select(-category) %>% 
  mutate(parent_plate_midpoint = (mean_pop_cat_value + LCR_plate_scaled_mean$mean_pop_cat_value) / 2) %>% 
  select(-mean_pop_cat_value)

# replace actual scaled plate mean with additive mean
plate_adjusted_hybrid_scaled_traits_pop_cat_mean.df <- hybrid_scaled_traits_ind.df %>% 
  group_by(FW_pop, Cross, trait, category) %>% 
  summarise(mean_pop_cat_value = mean(scaled_hybrid_value),
            family_size = n()) %>% 
  group_by(FW_pop, trait, category) %>%
  summarise(mean_pop_cat_value = weighted.mean(mean_pop_cat_value, family_size)) %>% 
  left_join(., parent_plate_midpt.df, by = "FW_pop") %>% 
  mutate(mean_pop_cat_value = ifelse(trait == "N.plates", parent_plate_midpoint, mean_pop_cat_value))


plate_adjusted_hybrid_mean_mm.df1 <- plate_adjusted_hybrid_scaled_traits_pop_cat_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ha = mean_pop_cat_value - P1) %>% # hybrid value minus anadromous parent value
  select(-mean_pop_cat_value, - P1) %>%  # remove hybrid value and anadromous parent
  mutate(ha = as.numeric(ha))

plate_adjusted_hybrid_mean_mm.df2 <- FW_scaled_parent_mean.df %>% 
  left_join(LCR_scaled_parent_means.df, by = c("trait")) %>% # include anadromous parent PC data by a common column called 'trait'
  filter(trait != "size_correction_Eye.Diam.log.res") %>% 
  mutate(ba = P2 - P1) %>% # ba = freshwater popuation parent value minus anadrmous parent value
  select(-P2, -P1)

# calculating multivariate mismatch for mean hybrids
plate_adjusted_mean_hybrid_mismatch.df <- plate_adjusted_hybrid_mean_mm.df1 %>% 
  left_join(plate_adjusted_hybrid_mean_mm.df2, by = c('FW_pop', 'trait')) %>% # joining 'ac' dataset by common columns 'FW_pop' and 'trait'
  group_by(FW_pop, category) %>% 
  mutate(t.across = (ha %*% ba) / (ba %*% ba)) %>% # each fish id has same t.across value
  mutate(d = ha - t.across * ba) %>% 
  summarize(plate_adjust_mean_pheno_hybrid_mismatch = sqrt(sum(d^2)))  # dist is mismatch here
# end plate adjusted 




multivar_mismatch_for_model <- mismatch_alltraits_analysis_df %>% 
  left_join(., p_mismatch_average.df, by = "FW_pop") %>% 
  mutate(excess_mismatch = hybrid_mismatch - average_p_mismatch,
         excess_mismatch_perpop = hybrid_mismatch - mean_pop_mm) %>% 
  filter(FishID %in% all_traits_measured_names$FishID) # only keep fish with all traits measured



# Does mismatch differ in hybrids among crosses? Or cross categories?
parent_hybrid_mismatch_cat.lmer <- lmer(hybrid_mismatch ~ FW_pop + category + (1|Cross), multivar_mismatch_for_model)
anova(parent_hybrid_mismatch_cat.lmer, ddf = "Kenward-Roger")

# This is the main model underlying Fig. 3A
mv_hyb_mm.lmer <- lmer(excess_mismatch ~ parent_trait_divergence + category + (1|Cross), data = multivar_mismatch_for_model)
anova_parent_hybrid_mismatch <- anova(mv_hyb_mm.lmer, ddf = "Kenward-Roger", type = 3)
main_multivar_mismatch.visreg <- visreg(mv_hyb_mm.lmer, "parent_trait_divergence", 
                                        by = "category", 
                                        overlay = TRUE,
                                        partial = TRUE)
summary(mv_hyb_mm.lmer)
emtrends(mv_hyb_mm.lmer, "category", var = "parent_trait_divergence")

r.squaredGLMM(mv_hyb_mm.lmer)
r2glmm::r2beta(mv_hyb_mm.lmer, partial = T)

# Seperate model with interaction term
mv_hyb_mm.lmer.int <- lmer(excess_mismatch ~ parent_trait_divergence * category + (1|Cross), data = multivar_mismatch_for_model)
emtrends(mv_hyb_mm.lmer.int, "category", var = "parent_trait_divergence")
anova_parent_hybrid_mismatch <- anova(mv_hyb_mm.lmer.int, ddf = "Kenward-Roger", type = 3)

# * * PAIRWISE MISMATCH ####
# F1 hybrids
pairwise_lms_multiple_regression_F1 <- pairwise_mismatch.df %>% 
  filter(category == "F1") %>% # only include F1
  ungroup() %>% # 
  na.omit() %>% # 
  # now generate means of each pop
  group_by(FW_pop, parent_trait_divergence, Group_w_TraitNames, Cross) %>% 
  summarise(mismatch = mean(mismatch)) %>% 
  group_by(FW_pop, parent_trait_divergence, Group_w_TraitNames, Cross) %>% 
  summarise(mismatch = mean(mismatch)) %>% 
  select(Group_w_TraitNames, mismatch, parent_trait_divergence) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_trait_divergence, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# F1 Pairwise analysis of slopes
slope_info_pairwise_multiple_F1 <- pairwise_lms_multiple_regression_F1 %>% 
  unnest(tidied) %>% 
  filter(term == "parent_trait_divergence") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(category = "F1")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_F1$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_F1$p.value)
sum(ifelse(slope_info_pairwise_multiple_F1$p.value < 0.05, T, F)) / length(slope_info_pairwise_multiple_F1$p.value)


# frac positive | sig
slope_info_pairwise_multiple_F1 %>% 
  filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>%
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# F2 hybrids
pairwise_lms_multiple_regression_F2 <- pairwise_mismatch.df %>% 
  filter(category == "F2") %>% # only include F2
  ungroup() %>% # 
  na.omit() %>% # 
  # now generate means of each pop
  group_by(FW_pop, parent_trait_divergence, Group_w_TraitNames, Cross) %>% 
  summarise(mismatch = mean(mismatch)) %>% 
  group_by(FW_pop, parent_trait_divergence, Group_w_TraitNames, Cross) %>% 
  summarise(mismatch = mean(mismatch)) %>% 
  select(Group_w_TraitNames, mismatch, parent_trait_divergence) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_trait_divergence, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# F2 Pairwise analysis of slopes
slope_info_pairwise_multiple_F2 <- pairwise_lms_multiple_regression_F2 %>% 
  unnest(tidied) %>% 
  filter(term == "parent_trait_divergence") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(category = "F2")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_F2$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_F2$p.value)
sum(ifelse(slope_info_pairwise_multiple_F2$p.value < 0.05, T, F)) / length(slope_info_pairwise_multiple_F2$p.value)


# frac positive | sig
slope_info_pairwise_multiple_F2 %>% 
  filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>%
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# put all the pairwise info into slopes
pairwise_mismatch_slopes.df <- bind_rows(slope_info_pairwise_multiple_F1, slope_info_pairwise_multiple_F2)

# * F1 PW Dom MM Effect ####
# Pairwise effects of dominance on mismatch
pairwise_lms_multiple_regression_median_hyb <- pairwise_mismatch.df %>% 
  filter(is.na(category)) %>% 
  select(-FW_pop, -category, -Cross) %>% 
  separate(FishID, into = c("type", "FW_pop", "category")) %>%
  filter(category == "F1") %>% 
  select(-type) %>% 
  ungroup() %>% # 
  na.omit() %>% # 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_dist_focal_pair, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# Median Hybrid Pairwise analysis of slopes
slope_info_pairwise_multiple_median <- pairwise_lms_multiple_regression_median_hyb %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(type = "median")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_median$p.value < 0.05, T, F))
sum(ifelse(slope_info_pairwise_multiple_median$estimate > 0, T, F)) 

# frac positive | sig
slope_info_pairwise_multiple_median %>% 
  filter(p.value < 0.05) %>% 
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# * F1 PW Var Mm Effect ####
pairwise_lms_multiple_regression_var_effect <- variance_effect_pairwise %>% 
  filter(category == "F1") %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(Group_w_TraitNames) %>% 
  select(Group_w_TraitNames, variance_effect, parent_dist_focal_pair, Cross) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(variance_effect ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment) # augment
  ) 

# F2 slope analysis
slope_info_pairwise_multiple_var_effect <- pairwise_lms_multiple_regression_var_effect %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value) %>% 
  mutate(category = "F1")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

sum(ifelse(slope_info_pairwise_multiple_var_effect$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_var_effect$p.value < 0.05)

slope_info_pairwise_multiple_var_effect %>% 
  filter(p.value < 0.05) %>% 
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))


# * F2 PW Dom MM Effect ####
pairwise_lms_multiple_regression_median_hyb <- pairwise_mismatch.df %>% 
  filter(is.na(category)) %>% 
  select(-FW_pop, -category, -Cross) %>% 
  separate(FishID, into = c("type", "FW_pop", "category")) %>%
  filter(category == "F2") %>% 
  select(-type) %>% 
  ungroup() %>% # 
  na.omit() %>% # 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_dist_focal_pair, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# Median Hybrid Pairwise analysis of slopes
slope_info_pairwise_multiple_median <- pairwise_lms_multiple_regression_median_hyb %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(type = "median")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_median$p.value < 0.05, T, F))
sum(ifelse(slope_info_pairwise_multiple_median$estimate > 0, T, F))

# frac positive | sig
slope_info_pairwise_multiple_median %>% 
  # filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# * F2 PW Var MM Effect ####

pairwise_lms_multiple_regression_var_effect <- variance_effect_pairwise %>% 
  filter(category == "F2") %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(Group_w_TraitNames) %>% 
  select(Group_w_TraitNames, variance_effect, parent_dist_focal_pair, Cross) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(variance_effect ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment) # augment
  ) 

# F2 slope analysis
slope_info_pairwise_multiple_var_effect <- pairwise_lms_multiple_regression_var_effect %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value) %>% 
  mutate(category = "F2")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

sum(ifelse(slope_info_pairwise_multiple_var_effect$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_var_effect$p.value < 0.05)

slope_info_pairwise_multiple_var_effect %>% 
  # filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))



multivar_mismatch_for_model <- mismatch_alltraits_analysis_df %>% 
  left_join(., p_mismatch_average.df, by = "FW_pop") %>% 
  mutate(excess_mismatch = hybrid_mismatch - average_p_mismatch,
         excess_mismatch_perpop = hybrid_mismatch - mean_pop_mm) %>% 
  filter(FishID %in% all_traits_measured_names$FishID) # only keep fish with all traits measured



# Does mismatch differ in hybrids among crosses? Or cross categories?
parent_hybrid_mismatch_cat.lmer <- lmer(hybrid_mismatch ~ FW_pop + category + (1|Cross), multivar_mismatch_for_model)
anova(parent_hybrid_mismatch_cat.lmer, ddf = "Kenward-Roger")

# This is the main model underlying Fig. 3A
mv_hyb_mm.lmer <- lmer(excess_mismatch ~ parent_trait_divergence + category + (1|Cross), data = multivar_mismatch_for_model)
anova_parent_hybrid_mismatch <- anova(mv_hyb_mm.lmer, ddf = "Kenward-Roger", type = 3)
main_multivar_mismatch.visreg <- visreg(mv_hyb_mm.lmer, "parent_trait_divergence", 
                                        by = "category", 
                                        overlay = TRUE,
                                        partial = TRUE)
summary(mv_hyb_mm.lmer)
emtrends(mv_hyb_mm.lmer, "category", var = "parent_trait_divergence")

r.squaredGLMM(mv_hyb_mm.lmer)
r2glmm::r2beta(mv_hyb_mm.lmer, partial = T)

# Seperate model with interaction term
mv_hyb_mm.lmer.int <- lmer(excess_mismatch ~ parent_trait_divergence * category + (1|Cross), data = multivar_mismatch_for_model)
emtrends(mv_hyb_mm.lmer.int, "category", var = "parent_trait_divergence")
anova_parent_hybrid_mismatch <- anova(mv_hyb_mm.lmer.int, ddf = "Kenward-Roger", type = 3)

# * * PAIRWISE MISMATCH ####
# F1 hybrids
pairwise_lms_multiple_regression_F1 <- pairwise_mismatch.df %>% 
  filter(category == "F1") %>% # only include F1
  ungroup() %>% # 
  na.omit() %>% # 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair, Cross) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(mismatch ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# F1 Pairwise analysis of slopes
slope_info_pairwise_multiple_F1 <- pairwise_lms_multiple_regression_F1 %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(category = "F1")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_F1$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_F1$p.value)
sum(ifelse(slope_info_pairwise_multiple_F1$p.value < 0.05, T, F)) / length(slope_info_pairwise_multiple_F1$p.value)


# frac positive | sig
slope_info_pairwise_multiple_F1 %>% 
  filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>%
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

## F2 hybrids
pairwise_lms_multiple_regression_F2 <- pairwise_mismatch.df %>% 
  filter(category == "F2") %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(Group_w_TraitNames) %>% 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair, Cross) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(mismatch ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment) # augment
  ) 

# F2 slope analysis
slope_info_pairwise_multiple_F2 <- pairwise_lms_multiple_regression_F2 %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value) %>% 
  mutate(category = "F2")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

sum(ifelse(slope_info_pairwise_multiple_F2$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_F2$p.value < 0.05)

slope_info_pairwise_multiple_F2 %>% 
  # filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  # group_by(positive_slope) %>%
  ungroup() %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# put all the pairwise info into slopes
pairwise_mismatch_slopes.df <- bind_rows(slope_info_pairwise_multiple_F1, slope_info_pairwise_multiple_F2)

# * F1 PW Dom MM Effect ####
# Pairwise effects of dominance on mismatch
pairwise_lms_multiple_regression_median_hyb <- pairwise_mismatch.df %>% 
  filter(is.na(category)) %>% 
  select(-FW_pop, -category, -Cross) %>% 
  separate(FishID, into = c("type", "FW_pop", "category")) %>%
  filter(category == "F1") %>% 
  select(-type) %>% 
  ungroup() %>% # 
  na.omit() %>% # 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_dist_focal_pair, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# Median Hybrid Pairwise analysis of slopes
slope_info_pairwise_multiple_median <- pairwise_lms_multiple_regression_median_hyb %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(type = "median")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_median$p.value < 0.05, T, F))
sum(ifelse(slope_info_pairwise_multiple_median$estimate > 0, T, F)) 

# frac positive | sig
slope_info_pairwise_multiple_median %>% 
  filter(p.value < 0.05) %>% 
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# * F1 PW Var Mm Effect ####
pairwise_lms_multiple_regression_var_effect <- variance_effect_pairwise %>% 
  filter(category == "F1") %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(Group_w_TraitNames) %>% 
  select(Group_w_TraitNames, variance_effect, parent_dist_focal_pair, Cross) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(variance_effect ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment) # augment
  ) 

# F2 slope analysis
slope_info_pairwise_multiple_var_effect <- pairwise_lms_multiple_regression_var_effect %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value) %>% 
  mutate(category = "F1")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

sum(ifelse(slope_info_pairwise_multiple_var_effect$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_var_effect$p.value < 0.05)

slope_info_pairwise_multiple_var_effect %>% 
  filter(p.value < 0.05) %>% 
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))


# * F2 PW Dom MM Effect ####
pairwise_lms_multiple_regression_median_hyb <- pairwise_mismatch.df %>% 
  filter(is.na(category)) %>% 
  select(-FW_pop, -category, -Cross) %>% 
  separate(FishID, into = c("type", "FW_pop", "category")) %>%
  filter(category == "F2") %>% 
  select(-type) %>% 
  ungroup() %>% # 
  na.omit() %>% # 
  select(Group_w_TraitNames, mismatch, parent_dist_focal_pair) %>%
  group_by(Group_w_TraitNames) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(mismatch ~ parent_dist_focal_pair, data = .x, na.action = na.omit)), # linear mixed model 
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment)# augment
  )

# Median Hybrid Pairwise analysis of slopes
slope_info_pairwise_multiple_median <- pairwise_lms_multiple_regression_median_hyb %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value)%>% 
  mutate(type = "median")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

# Fraction of slopes significant and positive
sum(ifelse(slope_info_pairwise_multiple_median$p.value < 0.05, T, F))
sum(ifelse(slope_info_pairwise_multiple_median$estimate > 0, T, F))

# frac positive | sig
slope_info_pairwise_multiple_median %>% 
  # filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

# * F2 PW Var MM Effect ####

pairwise_lms_multiple_regression_var_effect <- variance_effect_pairwise %>% 
  filter(category == "F2") %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(Group_w_TraitNames) %>% 
  select(Group_w_TraitNames, variance_effect, parent_dist_focal_pair, Cross) %>% 
  nest(-Group_w_TraitNames) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(variance_effect ~ parent_dist_focal_pair + (1|Cross), data = .x, na.action = na.omit)), # linear mixed model
    tidied = purrr::map(fit, tidy), # tidy
    glanced = purrr::map(fit, glance),# glanced
    augmented = purrr::map(fit, augment) # augment
  ) 

# F2 slope analysis
slope_info_pairwise_multiple_var_effect <- pairwise_lms_multiple_regression_var_effect %>% 
  unnest(tidied) %>% 
  filter(term == "parent_dist_focal_pair") %>%  # include parent trait divergence only
  select(Group_w_TraitNames, estimate, p.value) %>% 
  mutate(category = "F2")
# filter(p.value < 0.05) %>% # only include trait pairs for which p value is less than 0.05

sum(ifelse(slope_info_pairwise_multiple_var_effect$p.value < 0.05, T, F))
length(slope_info_pairwise_multiple_var_effect$p.value < 0.05)

slope_info_pairwise_multiple_var_effect %>% 
  # filter(p.value < 0.05) %>%
  mutate(positive_slope = ifelse(estimate > 0, T, F)) %>% 
  group_by(positive_slope) %>% 
  summarise(n=n(),
            mean_abs_slope = mean(abs(estimate)))

