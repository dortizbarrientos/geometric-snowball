# 08_repeatability
library(broom)
library(purrr)
library(tidyverse)

linear_msrs_1 <- read_csv('data/raw_data/raw_msr_set_1.csv')
linear_msrs_2 <- read_csv('data/raw_data/repeatability/ken_repeatability_fish_2.csv')
scope_measures_original_additional <- read_csv('data/raw_data/raw_msr_set_2.csv')
ken_fins <- read_csv('data/raw_data/repeatability/ken_repeatability_fins.csv')
ken_scope <- read_csv('data/raw_data/repeatability/ken_repeatability_scope.csv')


# linear measurement datasets

scope_measures_original_additional_long <- scope_measures_original_additional %>% 
  select(TubeLabel, N.plates:GillRakerL) %>% 
  rename(FishID = TubeLabel,
         Dorsal.Fin.Rays = Dorsal.Ray,
         Anal.Fin.Rays = Anal.Ray,
         Gill.Raker.Ct.OuterArch.Lower = GillRakerNo,
         Gill.raker.lt_um_40x_mag = GillRakerL) %>% 
  pivot_longer(cols = N.plates:Gill.raker.lt_um_40x_mag, names_to = "trait", values_to = "first_measurement") %>% 
  na.omit()

scope_msrs_long_complete <- linear_msrs_1 %>% 
  filter(FishID %in% linear_msrs_2$fish) %>% 
  select(FishID, N.plates, Dorsal.Fin.Rays:Gill.raker.lt_um_40x_mag) %>% 
  pivot_longer(cols = N.plates:Gill.raker.lt_um_40x_mag, names_to = "trait", values_to = "first_measurement") %>% 
  bind_rows(., scope_measures_original_additional_long)


linear_msrs_1_long <- linear_msrs_1 %>% 
  filter(FishID %in% linear_msrs_2$fish) %>% 
  mutate(Eye.Diam = eye.measurements.2) %>% 
  select(FishID, Standard.Lt:Pelvic.Spine.Lt) %>% 
  pivot_longer(cols = Standard.Lt:Pelvic.Spine.Lt, names_to = "trait", values_to = "first_measurement")

fins_processed <- ken_fins %>% 
  pivot_wider(names_from = rep, values_from = length_mm) %>% 
  rename(first_measurement = `1`, 
         second_measurement = `2`) %>% 
  select(-msr_order) %>% 
  mutate(trait = "pectoral_fin_length")

scope_original_and_repeat <- ken_scope %>% 
  pivot_longer(N.plates:Gill.raker.lt_um_40x_mag, names_to = "trait", values_to = "second_measurement") %>% 
  na.omit() %>% 
  rename(FishID = fish) %>% 
  left_join(., scope_msrs_long_complete, by = c("FishID", "trait"))


repeatability_df <- linear_msrs_2 %>% 
  pivot_longer(cols = Standard.Lt:Pelvic.Spine.Lt, names_to = "trait", values_to = "second_measurement") %>% 
  rename(FishID = fish) %>% 
  left_join(., linear_msrs_1_long, by = c("FishID", "trait")) %>% 
  na.omit() %>% 
  bind_rows(., fins_processed, scope_original_and_repeat)

# do the repeatability regressions for all traits! evaluate the stats
fitted_models <- repeatability_df %>% 
  group_by(trait) %>%
  do(repeat_lm = cor.test(.$second_measurement, .$first_measurement, method = "pearson")) %>% 
  ungroup() %>% 
  mutate(repeat_coef = map(repeat_lm, tidy)) %>% 
  unnest(repeat_coef) 

write_csv(repeatability_df, 'data/derived_data/repeatability_df.csv')

