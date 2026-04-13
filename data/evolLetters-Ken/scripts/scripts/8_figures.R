# Script

# LAOD PACKAGES ####
library(cowplot)
library(grid) # for figure annotations
library(ggpubr)
library(ggspatial)
library(lme4)
library(lmerTest)
library(rgeos)
library(rnaturalearth)
library(sf)
library(tidyverse) # to do relating to data cleaning
library(xtable)

# load datasets ####
combined_dataset_dominance <- read_csv('data/derived_data/combined_dominance_df.csv')
correlation_matrix_summaries <- read_csv('data/derived_data/correlation_matrix_summaries.csv') # correlation coefs and p vals
repeatability_df <- read_csv('data/derived_data/repeatability_df.csv')
mouth_correction_data_log <- read_csv('data/derived_data/size_mouth_corrected_pheno_data.csv')
trait_dominance_dataset <- read_csv('data/derived_data/trait_dominance_df.csv') # dominance of different traits; 
F1_standardized_dominance.df <- read.csv("data/derived_data/individual_f1_dominance_df.csv")
pairwise_mismatch.df <- read_csv('data/derived_data/pairwise_mismatch.csv')
pairwise_mismatch_slopes.df <- read_csv('data/derived_data/pairwise_mismatch_slopes.df.csv')
parent_dists.df <- read_csv('data/derived_data/multivariate_parent_dist.csv')
parent_means.df <- read_csv('data/derived_data/paren')
standardized_trait_df <- read_csv('data/derived_data/standardized_trait_df.csv')

# MAIN TEXT ####

# * Fig. 2A - Map ####
# Note: this map script only procudes a 'skeleton'. it was touched up quite a bit in Illustrator for the main fig.
world <- ne_countries(scale = "medium", returnclass = "sf")

world_points <- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

mappy<- par(usr=c(32.993358, -128.431168, 72.660698, -64.974139)) 

a1<- ggplot(data = world) + geom_point(x= -122.760440, y = 49.014735, size = 4, shape = 23, fill = "darkred") + 
  geom_sf() + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(axis.title.y = element_text(vjust = 0.2, size = 30),axis.title.x = element_text(vjust = 0.5, size = 30), plot.title = element_text(hjust = 0.5, face = "bold", size = 30)) + geom_sf(color = "black", fill = "antiquewhite")  + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + 
  annotation_scale(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in")) +
  coord_sf(xlim = c(-125.2, -122.40), ylim=c(48.5, 49.90), expand = FALSE) + 
  geom_text(data= world_points,aes(x=X, y=Y, label=name), color = "darkblue", fontface = "bold", check_overlap = TRUE) +
  annotate(geom = "text", x = -125.893867, y = 49.131450, label = "*",  fontface = "italic", fontface = "bold", color = "red", size = 6) +    
  annotate(geom = "text", x = -123.870369, y = 49.052725, label = "*",  fontface = "italic", fontface = "bold", color = "red", size = 6) +
  annotate(geom = "text", x = -123.178657, y = 49.194656, label = "*",  fontface = "italic", fontface = "bold", color = "red", size = 6) +
  annotate(geom = "text", x = -122.771928, y = 49.015127, label = "*",  fontface = "italic", fontface = "bold", color = "blue4", size = 6) +
  annotate(geom = "text", x =  -124.525167, y = 49.708722, label = "*", fontface = "italic", color = "blue4",fontface = "bold", size = 6, check_overlap = T) +
  annotate(geom = "text", x = -124.565571, y = 49.745935, label = "*", fontface = "italic", color = "blue4",fontface = "bold", size = 6) +
  annotate(geom = "text", x = -124.107367, y = 49.664196, label = "*", fontface = "italic",fontface = "bold", color = "blue4", size = 6) +
  annotate(geom = "text", x = -125.029256, y = 48.837784, label = "*", fontface = "italic",fontface = "bold", color = "blue4", size = 6) +
  annotate(geom = "text", x = -123.969894, y = 49.747864, label = "*", fontface = "italic", fontface = "bold",color = "blue4", size = 6) +
  annotate(geom = "text", x = -123.969556, y = 49.729643, label = "*", fontface = "italic",fontface = "bold", color = "blue4", size = 6) +
  annotate(geom = "text", x = -124.023049, y = 49.613181, label = "*", fontface = "italic",fontface = "bold", color = "blue4", size = 6) +
  annotate(geom = "text", x = -124.507107, y = 49.695493, label = "*", fontface = "italic", fontface = "bold",color = "blue4", size = 6, check_overlap = T) +
  annotate(geom = "text", x = -123.509106, y = 48.874348, label = "*", fontface = "italic", fontface = "bold",color = "blue4", size = 6)

# inset map
a2 <- ggplot(data = world) + geom_point(x= -122.760440, y = 49.014735, size = 4, shape = 23, fill = "darkred") + 
  geom_sf() + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  geom_sf(color = "black", fill = "antiquewhite") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-145.306168, -80), ylim=c(36.741868, 64.700696), expand = FALSE) + 
  geom_text(data= world_points,aes(x=X, y=Y, label=name), color = "black", fontface = "bold", check_overlap = TRUE) +
  geom_point(x = -123.509106, y = 48.874348, size = 6, shape = 21, fill = "darkred")

gg_inset_map1 <-  ggdraw() +
  draw_plot(a1) +
  draw_plot(a2, x = 0.71, y = 0.65, width = 0.3, height = 0.3)   

# the map was modified in adobe illustrator; so the code above just gives the 'skeleton' of this map

# * Fig. 3---Divergence-Mismatch Relationship ####
mean_p_mismatch <- mean(parent_mismatch_df$parent_mismatch)

multivariate_mismatch_divergence_plot <- ggplot(main_multivar_mismatch.visreg_res_DF, aes(x = parent_trait_divergence, y = visregRes, colour = category, fill = category)) + 
  xlab("phenotypic divergence between\nparent populations") + 
  ylab("mean mismatch\nin hybrids") +
  geom_smooth(data = main_multivar_mismatch.visreg_fit_DF %>% rename(visregRes = visregFit), aes(ymin = visregLwr, ymax = visregUpr, fill = category, colour = category), stat = "identity", alpha = 0.4) + 
  geom_point() + 
  scale_colour_manual(values = c("magenta", "darkorchid4")) +
  scale_fill_manual(values = c("magenta", "darkorchid4")) +
  theme_bw() + 
  geom_hline(yintercept = mean_p_mismatch, linetype = 2) + 
  theme(aspect.ratio = 0.55, 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") 

ggsave(multivariate_mismatch_divergence_plot, filename = '../PCD_polished_fig_panels/mismatch_divergence.pdf', units = "mm", height = 55, width = 95, useDingbats = F)

# * Fig 4---Causes of mismatch ####
# * * Dominance-effect ####
dominance_cause_mismatch.gg <- ggplot(data = dom_effect.visreg$res, aes(parent_trait_divergence, visregRes)) + 
  geom_point(aes(colour = category)) +
  geom_smooth(data = dom_effect.visreg$fit %>% rename(visregRes = visregFit), aes(ymin = visregLwr, ymax = visregUpr, fill = category, colour = category), stat = "identity", alpha = 0.25) + 
  scale_colour_manual(values = c("magenta", "darkorchid4")) +
  scale_fill_manual(values = c("magenta", "darkorchid4")) +
  xlab("phenotypic divergence between\nparent populations") +
  ylab("effect of dominance on\nmismatch in hybrids") + 
  theme_bw() + 
  theme(aspect.ratio = 0.67,
        legend.position = "none") 

# * * Variance-effect ####
variance_cause_mismatch.gg <- ggplot(data = ind_var_effect.visreg$res, aes(parent_trait_divergence, visregRes, colour = category)) +
  geom_smooth(data = ind_var_effect.visreg$fit %>% rename(visregRes = visregFit), aes(ymin = visregLwr, ymax = visregUpr, fill = category, colour = category), stat = "identity", alpha = 0.25) + 
  geom_point() +
  scale_colour_manual(values = c("magenta", "darkorchid4")) +
  scale_fill_manual(values = c("magenta", "darkorchid4")) +
  xlab("phenotypic divergence between\nparent populations") +
  ylab("mean effect of variation on\nmismatch in hybrids") + 
  theme_bw() +
  theme(aspect.ratio = 0.67,
        legend.position = "none")

height_effect_panels <- 2.6

ggsave(dominance_cause_mismatch.gg, filename = '../PCD_polished_fig_panels/effect_dom.pdf', height = height_effect_panels, width = 1.5*height_effect_panels)
ggsave(variance_cause_mismatch.gg, filename = '../PCD_polished_fig_panels/effect_var.pdf', height = height_effect_panels, width = 1.5*height_effect_panels)

#/////////////////////////////////////#
#### SUPPLEMENTARY ####
#/////////////////////////////////////#

# * Fig. S1---Repeatability ####

repeatability_names_vec <- sort(unique(repeatability_df$trait))
repeatability_codes_vec <- c("#AFR", "BD", "BW", "#DFR", "ED", "FDS", "GW", "#GR", "GRL", "HD", "#LAP", "PF", "PG", "PSL", "SDS", "SNT", "SL")
repeability_names_codes.tib <- tibble(trait = repeatability_names_vec, code = repeatability_codes_vec)

repeatability_df_codes <- repeatability_df %>% 
  left_join(., repeability_names_codes.tib, by = "trait")

repeatability_fig <- ggplot(repeatability_df_codes, aes(x = first_measurement, y = second_measurement)) + 
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_smooth(method = "lm", formula = 'y ~ x') + 
  facet_wrap(~code, scales = "free") + 
  xlab('first measurement') + 
  ylab('second measurement') + 
  theme(aspect.ratio = 1) + 
  theme_bw()

ggsave(repeatability_fig, filename = '../../../Apps/Overleaf/Thompson_KA_PhD_Thesis/figures/app_figs/D_1.pdf', width = 10, height = 10)

# * Fig. S2---Family FW-Marine Div ####
# Generate DF Of full names
pop_codes <- sort(parent_divergence.factor_levels_DF$FW_pop)
pop_names <- c("Bullock",
               "Cranby",
               "Klein",
               "Little Campbell River",
               "Little Quarry Benthic",
               "Little Quarry Limnetic",
               "North",
               "Paq",
               "Pachena",
               "Priest Benthic",
               "Priest Limnetic",
               "Paxton Benthic",
               "Paxton Limnetic")

pop_names_codes.tibble <- tibble(FW_pop = pop_codes, pop_names)


parent_divergence.aov.visregRes.DF_WNAMES <- parent_divergence.aov.visregRes.DF %>% 
  left_join(., pop_names_codes.tibble, by = "FW_pop")

parent_divergence.aov.visregFit.DF_WNAMES <- parent_divergence.aov.visregFit.DF %>% 
  left_join(., pop_names_codes.tibble, by = "FW_pop")

ptd_boxplot <- ggplot(parent_divergence.aov.visregRes.DF_WNAMES, aes(x = reorder(pop_names, visregRes, mean), y = visregRes)) + 
  geom_jitter(width = 0.25, alpha = 0.5, colour = "red") + 
  geom_point(parent_divergence.aov.visregFit.DF_WNAMES, mapping = aes(x = pop_names, y = visregFit), inherit.aes = T, colour = "black") + 
  geom_errorbar(parent_divergence.aov.visregFit.DF_WNAMES, mapping =  aes(x = pop_names, ymin = visregLwr, ymax = visregUpr), inherit.aes = F, colour = "black") +
  xlab('population') + 
  ylab('family mean distance to\nmean marine phenotype') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        aspect.ratio = 0.67) + 
  ylim(0, 11)

ggsave(ptd_boxplot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_2.pdf', height = 100, width = 120, units = "mm")

# * Fig. S3 - # div. traits ####
# * * A---Number ####
FW_pop <- "origin"
number_of_div_traits <- 0
parent_trait_divergence <- 0

num_divergent_ptd_stats.df <- t.test_dataset %>% 
  # select(FW_pop, trait, diff_stats) %>% 
  rowwise() %>% 
  mutate(nrow = nrow(.),
         sig_bonferroni = t_test_p < (0.05/nrow)) %>% 
  filter(sig_bonferroni == T) %>% 
  group_by(FW_pop) %>% 
  summarise(number_of_div_traits = n()) %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)) %>% 
  mutate(num_poss_mismatch = choose(number_of_div_traits, 2))

num_divergent_stats_only_fig <- ggplot(num_divergent_ptd_stats.df, aes(x = parent_trait_divergence, y = number_of_div_traits)) + 
  geom_point(size = 1) + 
  ylab('# divergent traits') + 
  xlab('phenotypic divergence between\nparent populations') + 
  geom_smooth(method = "lm") +
  theme_bw()

summary(lm(number_of_div_traits ~ parent_trait_divergence, data = num_divergent_ptd_stats.df))
anova(lm(number_of_div_traits ~ parent_trait_divergence, data = num_divergent_ptd_stats.df))

# * * B---Specific traits ####
diverg_ENCE_ENT_df <- t.test_dataset %>% 
  left_join(., parent_dists_df, by = "FW_pop") %>% 
  arrange(parent_trait_divergence) %>% 
  group_by(parent_trait_divergence) %>% 
  mutate(rank_div = cur_group_id()) %>% 
  left_join(., names_codes, by = "trait") %>% 
  ungroup() %>% 
  arrange(trait)




div_ence_ent_plot <- ggplot(diverg_ENCE_ENT_df %>% mutate(divergent = t_test_p < (0.05/nrow(.))), aes(x = rank_div, y = reorder(code, diff_stats, mean), fill = divergent)) + 
  geom_point(pch = 21, size = 2) + 
  xlab('rank order of phenotypic divergence') + 
  ylab('trait code') + 
  scale_y_discrete(limits=rev) + 
  scale_fill_manual(values = c("white", "black")) + 
  scale_x_continuous(breaks=1:12) +
  theme_bw()


div_traits_plot <- plot_grid(num_divergent_stats_only_fig, div_ence_ent_plot, labels = "AUTO", rel_widths = c(1, 1.5))

ggsave(div_traits_plot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_3.pdf', height = 4, width = 9)


# * Fig. S4: 2D mismatch ####
# * * A---PCH ####
input_data <- standardized_trait_df %>% 
  filter(abs(scaled_val) < 5) %>% 
  filter(trait %in% c("N.plates", "size_correction_Pelvic.Girdle.Lt.log.res")) %>% 
  select(FishID, Cross.Type, category, trait, scaled_val) %>% 
  pivot_wider(names_from = trait, values_from = scaled_val)  

# parent_mean_df
parent_means <- input_data %>% 
  filter(category %in% c("P1", "P2")) %>% 
  select(Cross.Type, N.plates, size_correction_Pelvic.Girdle.Lt.log.res) %>% 
  group_by(Cross.Type) %>% 
  summarise(N.plates = mean(N.plates, na.rm = T),
            RakerNo = mean(size_correction_Pelvic.Girdle.Lt.log.res, na.rm = T)) 

parent_means_PCH_lcr <- parent_means %>% 
  filter(Cross.Type %in% c("LCR-LCR", "PCH-PCH"))

plot_data_PCH_F2 <- input_data %>% 
  filter(Cross.Type == "LCR-PCH") %>% 
  select(N.plates, size_correction_Pelvic.Girdle.Lt.log.res) %>% 
  na.omit() %>% 
  rowwise() %>% 
  mutate(N.plates = N.plates + rnorm(1, 0, .2),
         size_correction_Pelvic.Girdle.Lt.log.res = size_correction_Pelvic.Girdle.Lt.log.res + rnorm(1, 0, .1))

# slope calc
m_p = (parent_means_PCH_lcr$N.plates[2] - parent_means_PCH_lcr$N.plates[1]) / (parent_means_PCH_lcr$RakerNo[2] - parent_means_PCH_lcr$RakerNo[1])
b_p = parent_means_PCH_lcr$N.plates[1] - m_p * parent_means_PCH_lcr$RakerNo[1]
m_h = - 1 / m_p

# df of line segments
segment_df_PCH <- plot_data_PCH_F2 %>% 
  rename(x = size_correction_Pelvic.Girdle.Lt.log.res,
         y = N.plates) %>% 
  # we want to make new line segments that START at the 2D phenotypes and end on the line connecting parents
  mutate(
    # slope of new line segments are perpendicular to old (negative reciprocal)
    m_h = m_h,
    # calculate y-intercept of new lines (y = mx + b; b = y - mx)
    b_h = y - m_h * x,
    # calculate the POINT OF INTERSECTION with the parent line
    # first calculate the x-coordinate
    xend = (b_p - b_h) / (m_h - m_p),
    yend = m_h * xend + b_h, 
    mismatch = sqrt((xend - x)^2 + (yend - y)^2))

mismatch_vis_plot_PCH <- ggplot(plot_data_PCH_F2, aes(x = size_correction_Pelvic.Girdle.Lt.log.res, y = N.plates)) + 
  geom_abline(slope = m_p, intercept = -b_p, size = 1.5, alpha = 0.5) + 
  # geom_abline(slope = -1/m, intercept = 0, size = 1.5, alpha = 0.75) + # visualizing perpindicular
  geom_segment(data = segment_df_PCH, aes(x, y, xend = xend, yend = yend), alpha = 0.75, linetype = 2) +   
  geom_point(colour = "magenta") + 
  geom_point(aes(x = parent_means_PCH_lcr$RakerNo[1], y = parent_means_PCH_lcr$N.plates[1]), colour = "red", size = 4) + # PCH parent
  geom_point(aes(x = parent_means_PCH_lcr$RakerNo[2], y = parent_means_PCH_lcr$N.plates[2]), colour = "blue", size = 4) +   # LCR parent
  xlim(c(1, -4.3)) +
  ylim(c(1.5, -2)) +
  xlab('length of pelvic girdle') +
  ylab('# lateral plates') + 
  theme_bw() + 
  ggtitle("Pachena Lake x Marine Hybrids") + 
  coord_fixed()
mismatch_vis_plot_PCH
mean(segment_df_PCH$mismatch)

# * * B---PXB ####

parent_means_PXB_lcr <- parent_means %>% 
  filter(Cross.Type %in% c("LCR-LCR", "PXB-PXB"))

plot_data_PXB_F2 <- input_data %>% 
  filter(Cross.Type == "LCR-PXB") %>% 
  select(N.plates, size_correction_Pelvic.Girdle.Lt.log.res) %>% 
  na.omit() %>% 
  rowwise() %>% 
  mutate(N.plates = N.plates + rnorm(1, 0, .2),
         size_correction_Pelvic.Girdle.Lt.log.res = size_correction_Pelvic.Girdle.Lt.log.res + rnorm(1, 0, .1))

# slope calc
m_p = (parent_means_PXB_lcr$N.plates[2] - parent_means_PXB_lcr$N.plates[1]) / (parent_means_PXB_lcr$RakerNo[2] - parent_means_PXB_lcr$RakerNo[1])
b_p = parent_means_PXB_lcr$N.plates[1] - m_p * parent_means_PXB_lcr$RakerNo[1]
m_h = - 1 / m_p

# df of line segments
segment_df_PXB <- plot_data_PXB_F2 %>% 
  rename(x = size_correction_Pelvic.Girdle.Lt.log.res,
         y = N.plates) %>% 
  # we want to make new line segments that START at the 2D phenotypes and end on the line connecting parents
  mutate(
    # slope of new line segments are perpendicular to old (negative reciprocal)
    m_h = m_h,
    # calculate y-intercept of new lines (y = mx + b; b = y - mx)
    b_h = y - m_h * x,
    # calculate the POINT OF INTERSECTION with the parent line
    # first calculate the x-coordinate
    xend = (b_p - b_h) / (m_h - m_p),
    yend = m_h * xend + b_h, 
    mismatch = sqrt((xend - x)^2 + (yend - y)^2))

mismatch_vis_plot_PXB <- ggplot(plot_data_PXB_F2, aes(x = size_correction_Pelvic.Girdle.Lt.log.res, y = N.plates)) + 
  geom_abline(slope = m_p, intercept = -b_p, size = 1.5, alpha = 0.5) + 
  # geom_abline(slope = -1/m, intercept = 0, size = 1.5, alpha = 0.75) + # visualizing perpindicular
  geom_segment(data = segment_df_PXB, aes(x, y, xend = xend, yend = yend), alpha = 0.75, linetype = 2) +   
  geom_point(colour = "magenta") + 
  geom_point(aes(x = parent_means_PXB_lcr$RakerNo[1], y = parent_means_PXB_lcr$N.plates[1]), colour = "red", size = 4) + # PXB parent
  geom_point(aes(x = parent_means_PXB_lcr$RakerNo[2], y = parent_means_PXB_lcr$N.plates[2]), colour = "blue", size = 4) +   # LCR parent
  xlim(c(1, -4.3)) +
  ylim(c(1.5, -2)) +
  xlab('length of pelvic girdle') +
  ylab('# lateral plates') + 
  coord_fixed() + 
  theme_bw() + 
  ggtitle("Paxton Lake Benthic x Marine Hybrids")
mismatch_vis_plot_PXB
mean(segment_df_PXB$mismatch)

mismatch_multiplot <- cowplot::plot_grid(mismatch_vis_plot_PCH, mismatch_vis_plot_PXB)
ggsave(mismatch_multiplot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_4.pdf', height = 3.5, width = 8)

# * Fig. S5---Dom ± SE across traits ####
dominance_across_traits <- ggplot(trait_mean_scaled_traits_for_dominance.df %>% filter(category == "F1"), aes(x = reorder(code, z_dom_scaled_val), y = z_dom_scaled_val)) + 
  geom_point() + 
  xlab(expression("trait")) + 
  ylab(expression("mean dominance coefficient ± SE")) +
  geom_errorbar(width = .3, aes(ymin = z_dom_scaled_val - se, 
                                ymax = z_dom_scaled_val + se), colour = "black") + 
  geom_hline(yintercept = 0.5, colour = "red", lty = 2) + 
  geom_hline(yintercept = 0, colour = "black", lty = 2) + 
  geom_hline(yintercept = 1, colour = "black", lty = 2) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5, axis.text.x.bottom = element_text(color = "black", angle = 40, size = 8, vjust = 1, hjust = 1)) + 
  annotate("label", label = "anadromous-like", size = 3.5, x = 13.6, y = 0.000, color = "black") + 
  annotate("label", label = "freshwater-like", size = 3.5, x = 14.0, y = 1.00, color = "black")
  # facet_wrap(~category, nrow = 2)

ggsave('../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_5.pdf', height = 8.5, width = 6.5)

# * Fig. S6---Dominance evolution ####
# * * PanelA---Plates ####
Cross_Type_FW_pop_Fish.df <- standardized_trait_df %>% 
  select(FishID, FW_pop, category, Cross) %>% 
  unique()

LCR_plate_sd.df <- standardized_trait_df %>% 
  filter(Cross.Type == "LCR-LCR",
         trait == "N.plates") %>% 
  summarise(lcr_plate_sd = sd(raw),
            lcr_plate_mean = mean(raw))

LCR_plates_mean <- LCR_plate_sd.df$lcr_plate_mean
LCR_plates_sd <- LCR_plate_sd.df$lcr_plate_sd

parent_dominance_means <- mouth_correction_data_log %>% 
  select(FishID, N.plates) %>% 
  left_join(., Cross_Type_FW_pop_Fish.df, by = "FishID") %>% 
  filter(category %in% c("P1", "P2")) %>% 
  na.omit() %>% 
  group_by(FW_pop, Cross) %>% 
  summarise(mean_cross_plates = mean(N.plates)) %>% 
  group_by(FW_pop) %>% 
  summarise(mean_fw_parent_plates = mean(mean_cross_plates)) %>% 
  filter(FW_pop != "NOR") %>% 
  arrange(mean_fw_parent_plates) %>% 
  mutate(row_number = row_number(),
         colour = ifelse(row_number %% 2 ==0, T, F),
         ptd_plate_sd = ((LCR_plates_mean - mean_fw_parent_plates) / LCR_plates_sd))

parental_plate_midpoints.df <- parent_dominance_means %>% 
  filter(FW_pop != "LCR") %>% 
  rowwise() %>% 
  mutate(midpoint = mean(c(mean_fw_parent_plates, LCR_plates_mean))) %>% 
  select(FW_pop, mean_fw_parent_plates, ptd_plate_sd, midpoint) %>% 
  mutate(colour = "unique")

Hybrid_plates_raw <- standardized_trait_df %>% 
  filter(trait == "N.plates",
         category %in% c("F1", "F2")) %>%
  select(FishID, category, FW_pop, raw) %>% 
  rename(N.plates = raw) %>% 
  left_join(., parent_dominance_means, by = "FW_pop") %>% 
  na.omit()

plates_evolution <- ggplot(Hybrid_plates_raw, aes(x = mean_fw_parent_plates, y = N.plates, colour = colour)) +
  geom_line(data = parental_plate_midpoints.df, aes(x = mean_fw_parent_plates, y = midpoint), size = 2) +
  geom_jitter(width = 0.1, size = 0.6) + 
  geom_smooth(method = "loess", colour = "red") + 
  scale_colour_manual(values = c("black", "grey59", "steelblue2")) +
  xlab('mean number of lateral armour plates\nin freshwater population ') + 
  ylab('number of lateral\narmour plates') +
  # ylab(expression(paste('F'[2],' hybrid lateral plate #'))) +
  theme_bw() + 
  theme(legend.position = "none",
        aspect.ratio = 1) + 
  facet_wrap(~category) + 
  scale_x_reverse()
plates_evolution

# * * PanelB---Dom Coefs ####
evolving_traits_df <- scaled_traits_for_dominance.df %>% 
  filter(trait %in% c("N.plates", "size_correction_Pelvic.Spine.Lt.log.res")) %>% 
  mutate(trait = ifelse(trait == "N.plates", "(i) number of lateral\narmour plates", "(ii) length of pelvic\nspine"),
         parent_trait_divergence = as.numeric(parent_trait_divergence))

dom_coefs_plates <- ggplot(evolving_traits_df %>% filter(category == "F1"), aes(x = parent_trait_divergence, y = z_dom_scaled_val, colour = category)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("magenta", "darkorchid4")) +
  geom_hline(yintercept = 0.5, colour = "red", lty = 2) + 
  geom_hline(yintercept = 0, colour = "black", lty = 2) + 
  geom_hline(yintercept = 1, colour = "black", lty = 2) + 
  ylab('dominance coefficient') + 
  xlab('phenotypic divergence between\nparent populations') + 
  facet_wrap(~trait) + 
  theme_bw() + 
  theme(legend.position = "none",
        aspect.ratio = 1)

dominance_evolution_plot <- cowplot::plot_grid(plates_evolution, dom_coefs_plates, labels = "AUTO", ncol = 1)
ggsave(dominance_evolution_plot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/figures/supp/PCD_S_6.pdf', units = "mm", useDingbats=FALSE, width = 200, height = 200)


# * Fig. S7 - Variance patterns ####
var_plot <- ggplot(variance_model.visreg$res, aes(x = parent_trait_divergence, y = visregRes, colour = category)) + 
  geom_jitter(alpha = 0.95, size = 0.75, width = 0.035) + 
  geom_smooth(data = variance_model.visreg$fit %>% rename(visregRes = visregFit), aes(ymin = visregLwr, ymax = visregUpr, fill = category, colour = category), stat = "identity") + 
  ylab("mean variance across traits") + 
  xlab("phenotypic distance from anadromous") +
  scale_colour_manual(values = c("magenta", "darkorchid4")) +
  scale_fill_manual(values = c("magenta", "darkorchid4")) +
  theme_bw() + 
  theme(aspect.ratio = 0.67,
        legend.position = "none")

ggsave(var_plot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_7.pdf', units = "mm", height = 82, width = 85)


# * Fig. S8: SNOWBALL ####
snowball_fig <- ggplot(pw_mismatch_diff_counts %>% filter(category != "P2"), aes(x = parent_trait_divergence, y = num_diff_and_bigger)) + 
  geom_point() + 
  xlim(0, 10.5) + 
  geom_smooth(method = "lm", formula = 'y ~ 0 + x', se = F) +
  geom_smooth(method = "lm", formula = 'y  ~ 0 + x + I(x^2)', se = F, colour = "red") +
  facet_wrap(~category) + 
  xlab('phenotypic distance from anadromous') + 
  ylab('# mismatched trait pairs') + 
  theme_KT_PCD

# STATS for the Fig
# F1s
mod1 <- lm(num_diff_and_bigger ~ 0 + parent_trait_divergence, data = pw_mismatch_diff_counts %>% filter(category == "F1"))
mod2 <- lm(num_diff_and_bigger ~ 0 + parent_trait_divergence + I(parent_trait_divergence^2), data = pw_mismatch_diff_counts %>% filter(category == "F1"))
AIC(mod1)
AIC(mod2)

# F2s

mod1 <- lm(num_diff_and_bigger ~ 0 + parent_trait_divergence, data = pw_mismatch_diff_counts %>% filter(category == "F2"))
mod2 <- lm(num_diff_and_bigger ~ 0 + parent_trait_divergence + I(parent_trait_divergence^2), data = pw_mismatch_diff_counts %>% filter(category == "F2"))
AIC(mod1)
AIC(mod2)

ggsave(snowball_fig, filename = '../../../Apps/Overleaf/Thompson_KA_PhD_Thesis/figures/app_figs/D_8.pdf', height = 4, width = 7)

# * Fig. S9 - Pairwise Mismatch Patterns ####
pairwise_mismatch_slopes.df_F1 <- pairwise_mismatch_slopes.df %>% filter(category == "F1")
pairwise_mismatch_slopes.df_F2 <- pairwise_mismatch_slopes.df %>% filter(category == "F2") %>% mutate(sig = ifelse(p.value < 0.05, T, F))

pw_mismatch_ptd_estimate_plot <- ggplot(pairwise_mismatch_slopes.df %>% mutate(sig = ifelse(p.value < 0.05, T, F)), aes(x = estimate, fill = category, alpha = sig)) +
  geom_histogram() + 
  scale_alpha_discrete(range = c(.4, 1)) + 
  xlab('effect of divergence on mismatch (pairwise)') + 
  ylab("frequency") + 
  scale_fill_manual(values = c("magenta", "darkorchid4")) +
  theme_bw()  + 
  geom_vline(xintercept = 0) + 
  scale_x_continuous(breaks = seq(-0.3, 0.4, by = 0.15)) + 
  scale_y_continuous(breaks = seq(0, 12, by = 3)) + 
  theme(aspect.ratio = 0.5, 
        # strip.background = element_blank(),
        # strip.text.x = element_blank(),
        legend.position = "none") + 
  facet_wrap(~category)

ggsave(pw_mismatch_ptd_estimate_plot, filename = '../../../Apps/Overleaf/PCD_Phenotypic_Consequences_of_Divergence/Figures/supp/PCD_S_9.pdf', units = "mm", height = 50, width = 180)
