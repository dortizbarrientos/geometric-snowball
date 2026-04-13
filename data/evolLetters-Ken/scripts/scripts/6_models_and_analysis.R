# Data analysis

library(broom) # elegent tools to deal with complex objects
library(broom.mixed) # elegent tools to deal with complex objects (mixed models)
library(emmeans) # extract estimated marginal means and CIs
library(lme4) # run mixed models
library(lmerTest) # test mixed models
library(MuMIn) # 
library(pspearman) # spearman tests
library(tidyverse) # data processing tools
library(visreg) # visualize regressions

# read data ####
family_parent_dists_df_w_info <- read_csv("data/family_parent_dists_df_w_info.csv")
mismatch_alltraits_analysis_df <- read_csv(file = 'data/derived_data/mismatch_alltraits_analysis_df.csv')

# 
parent_divergence.aov <- aov(parent_trait_divergence ~ FW_pop, family_parent_dists_df_w_info)
parent_divergence.aov.visreg <- visreg(parent_divergence.aov)
# Df for vector of ordered factor levels
parent_divergence.factor_levels_DF <- family_parent_dists_df_w_info %>% 
  group_by(FW_pop) %>% 
  summarise(mean_dist = mean(parent_trait_divergence)) %>% 
  arrange(mean_dist)

parent_divergence.aov.visregRes.DF <- parent_divergence.aov.visreg$res %>% 
  as_tibble() %>% 
  mutate(FW_pop = fct_reorder(FW_pop, visregRes, .fun='mean'))
parent_divergence.aov.visregFit.DF <- parent_divergence.aov.visreg$fit 

# Mismatch ####
# * MULTIVARIATE analysis ####
# Mismatch by cateogry 

parent_mismatch_unique.df <- mismatch_alltraits_analysis_df %>% 
  select(parent_trait_divergence, parent_mismatch) %>% 
  unique()
  
# * Parent mismatch variation among populations
ind_parent_by_pop.lmer <- lmer(mismatch ~ FW_pop + (1|Cross), data = ind_parent_mismatch)
anova(ind_parent_by_pop.lmer, ddf = "Kenward-Roger")

p_mismatch_average.df <- ind_parent_mismatch %>% 
  group_by(Cross, FW_pop) %>% 
  summarise(mean_cross_mm = mean(mismatch)) %>% 
  group_by(FW_pop) %>% 
  summarise(mean_pop_mm = mean(mean_cross_mm))

average_p_mismatch <- mean(p_mismatch_average.df$mean_pop_mm)



# * Main multivariate mismatch hypothesis ####
popmeans_multivar_mismatch_for_model <- mismatch_alltraits_analysis_df %>% 
  group_by(Cross, FW_pop, category) %>% 
  summarise(n = n(), 
            mean_fam_mismatch = mean(hybrid_mismatch),
            parent_trait_divergence = mean(parent_trait_divergence)) %>% 
  group_by(FW_pop, category) %>% 
  summarise(mean_mismatch = weighted.mean(mean_fam_mismatch, n),
            parent_trait_divergence = mean(parent_trait_divergence))

# This is the main model underlying Fig. 3A
mv_hyb_mm.lmer <- lmer(mean_mismatch ~ parent_trait_divergence + category + (1|FW_pop), data = popmeans_multivar_mismatch_for_model)
anova(mv_hyb_mm.lmer, ddf = "Kenward-Roger")
summary(mv_hyb_mm.lmer)
main_multivar_mismatch.visreg <- visreg(mv_hyb_mm.lmer, "parent_trait_divergence", 
                                        by = "category", 
                                        overlay = TRUE,
                                        partial = TRUE)

main_multivar_mismatch.visreg_res_DF <- main_multivar_mismatch.visreg$res
main_multivar_mismatch.visreg_fit_DF <- main_multivar_mismatch.visreg$fit


mv_hyb_mm.lm <- lm(mean_excess_mm ~ category + parent_mismatch + parent_trait_divergence, data = popmeans_multivar_mismatch_for_model)
car::Anova(mv_hyb_mm.lm, type = 3)
summary(mv_hyb_mm.lm)
visreg(mv_hyb_mm.lm)


mv_hyb_mm.lm <- lm(mean_excess_mm ~ parent_trait_divergence, data = popmeans_multivar_mismatch_for_model)
summary(mv_hyb_mm.lm)
###


# EFFECT OF DOMINANCE ON MISMATCH ####
dom_effect_mm.lmer <- lmer(mean_pheno_hybrid_mismatch ~ parent_trait_divergence * category + (1|FW_pop), data = mean_hybrid_mismatch.df)
anova(dom_effect_mm.lmer, type = 3, ddf = "Kenward-Roger")
summary(dom_effect_mm.lmer)
emtrends(dom_effect_mm.lmer, pairwise~category, var = "parent_trait_divergence", infer = T)

dom_effect.visreg <- visreg(dom_effect_mm.lmer, "parent_trait_divergence", by = 'category', overlay = T)


# * Dominance among populations within traits ####
trait_vec <- unique(raw_scaled_traits_for_dominance.df$trait)

magic_for()
for (i in 1:length(trait_vec)) {
  f1_trait_df_within_loop <- raw_scaled_traits_for_dominance.df_F1 %>% 
    filter(trait == trait_vec[i]) 
  
  model_in_loop <- lmer(z_dom_scaled_val ~ FW_pop + (1|Cross), data = f1_trait_df_within_loop)
  Tukey_in_loop <- emmeans(model_in_loop, list(pairwise ~ FW_pop), adjust = "tukey") 
  Tukey_in_loop_tibble <- Tukey_in_loop$`pairwise differences of FW_pop` %>% 
    as_tibble()
  frac_tukey_sig <- mean(Tukey_in_loop_tibble$p.value < 0.05)
  trait <- trait_vec[i]
  
  put(trait, frac_tukey_sig)
}
Tukey_Fracs.df <- magic_result_as_dataframe()
# mean fraction of pairwise differences that are significant
mean(Tukey_Fracs.df$frac_tukey_sig)

# DOMINANCE x DIVERGENCE ####
# is there an association between dominance + divergence?
# (one datum per pop)
dominance_divergence_regression <- scaled_traits_for_dominance.df %>%
  select(trait, FW_pop, category, z_dom_scaled_val, parent_trait_divergence) %>% 
  group_by(trait) %>% 
  nest() %>%  
  mutate(regression = map(data, ~lm(z_dom_scaled_val ~ parent_trait_divergence * category, data = .))) %>%
  mutate(tidy = map(regression, broom::tidy)) %>% 
  unnest(tidy) %>% 
  select(-data, -regression) %>% 
  filter(!term %in% c("(Intercept)", "categoryF2")) %>% 
  filter(term == "parent_trait.divergence") %>% 
  mutate(different = ifelse(p.value < 0.05, T, F))
  
# trait by trait
dominance_divergence_model_trait_RF.lmer <- lmer(z_dom_scaled_val ~ parent_trait_divergence * category + (1|trait) + (1|FW_pop), data = scaled_traits_for_dominance.df %>% mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)))
summary(dominance_divergence_model_trait_RF.lmer)
anova(dominance_divergence_model_trait_RF.lmer, type = 3)

# VARIANCE x DIVERGENCE ####
variance_model <- lm(mean_variance ~ parent_trait_divergence * category, data = mean_variances)
car::Anova(variance_model, type = 3)
variance_model.visreg <- visreg(variance_model, "parent_trait_divergence", "category", overlay = T)

# get the slopes of the two groups
emtrends(variance_model, ~category, var = "parent_trait_divergence", infer = T)

# trait by trait
variance_model_trait_RF.lmer <- lmer(var ~ parent_trait_divergence * category + (1|trait), data = individual_trait_variances.df)
anova(variance_model_trait_RF.lmer, type = 3)
# variance_model.visreg <- visreg(variance_model_trait_RF.lmer, "parent_trait_divergence", "category", overlay = T)

# get the slopes of the two groups
emtrends(variance_model, ~category, var = "parent_trait_divergence")

# * mismatch effect of variance ####
family_effect_of_var_on_hybrid_mismatch_mean.df <- family_effect_of_var_on_hybrid_mismatch.df %>% 
  group_by(FW_pop, category, Cross, parent_trait_divergence) %>% 
  summarise(var_caused_mismatch = mean(var_caused_mismatch)) %>% 
  group_by(FW_pop, category, parent_trait_divergence) %>% 
  summarise(var_caused_mismatch = mean(var_caused_mismatch))

var_effect_mismatch.lmer <- lmer(var_caused_mismatch ~ parent_trait_divergence * category + (1|FW_pop), family_effect_of_var_on_hybrid_mismatch_mean.df)
anova(var_effect_mismatch.lmer, type = 3, ddf = "Kenward-Roger")
summary(var_effect_mismatch.lmer)
emtrends(var_effect_mismatch.lmer, pairwise~category, var = "parent_trait_divergence", infer = T)
ind_var_effect.visreg <- visreg(var_effect_mismatch.lmer, "parent_trait_divergence", by = "category", overlay = T)

# write files ####
write_csv(pairwise_mismatch_slopes.df, "data/derived_data/pairwise_mismatch_slopes.df.csv")

          