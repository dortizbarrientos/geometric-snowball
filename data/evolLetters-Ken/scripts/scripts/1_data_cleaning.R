# PCD: Phenotypic consequences of divergence
# Script authors: Ken A. Thompson & Avneet K. Chhina
# This script takes raw phenotype datasets and cleans/standardizes them for use in subsequent analyses
# Script 1 of 9

# PACKAGES ####
library(tidyverse) # to do everything relating to data cleaning

# function to calculate parent euclidean distance
parent_euclidean <- function(df, which_LCR) {
  # vector of p1 vals
  x1 = as_vector(df[2:ncol(df)]) 
  rbind(x1, which_LCR) %>% 
    dist(.) # euclidean distance of parents
}

# Critical values
replace_0_val <- 0.1
fish_to_exclude <- c("KT_Gac_212", "KT_Gac_818", "KT_Gac_896",  "KT_Gac_1806", "KT_Gac_1807", "KT_Gac_1186", "KT_Gac_201", "KT_Gac_218", "KT_Gac_3387", "KT_Gac_2426") # these fish had broken second dorsal spines and/or failed to inflate their swim bladder

# create vectors for naming crosses
f1_CrossName_vec <- c("LCR-LQL", "LCR-PLB", "LCR-PLL", 
                      "LCR-KLN", "LCR-CRN", "LCR-LQB",
                      "LCR-BUL", "LCR-NOR", "LCR-PAQ",
                      "LCR-PCH", "BUL-LCR", "LCR-PXL", 
                      "LCR-PXB")

p2_CrossName_vec <- c("LQL-LQL", "LQB-LQB", "PCH-PCH",
                      "PXL-PXL", "PLB-PLB", "PXB-PXB",
                      "KLN-KLN", "PAQ-PAQ", "LQL-LQL", 
                      "BUL-BUL", "PLL-PLL", "NOR-NOR",
                      "CRN-CRN", "PLL-PLL_F2")

f2_CrossName_vec <- c("LCR-LQB_F2", "LCR-CRN_F2", "LCR-KLN_F2",
                      "LCR-PCH_F2", "LCR-PAQ_F2", "LCR-NOR_F2",
                      "LCR-PLB_F2", "LCR-LQL_F2", "LCR-BUL_F2",
                      "LCR-PLL_F2", "LCR-PXB_F2", "LCR-PXL_F2")

# READ DATA ####
F1_P1_P2_dataset <- read.csv("data/raw_data/raw_msr_set_1.csv") # mostly contains measurements for P1 (anadromous),P2 (freshwater), and F1 hybrids
F2_P1_P2_dataset <- read.csv("data/raw_data/raw_msr_set_2.csv") # mostly contains measurements for P1, P2, F2 
linear_measurement_dataset <- read.csv("data/raw_data/raw_additional_linear_msr.csv") # contains linear measuremnts of fish
additional_measurements_dataset <- read.csv("data/raw_data/raw_msr_set_3.csv") # four set of raw fish measurements

# DATA CLEANING ####
# DF with linear measurement trait names
linear_msr_info_df <- linear_measurement_dataset %>% 
  fill(FishID) %>% #fills the empty cells with corresponding Fish IDs (fish names).
  filter(FishID == "KT_Gac_2096") %>% # picking one fish with complete measurements to retain trait names
  select(trait, number) #only trait and number is selected in the dataset

# DF containing info for all fish
linear_measurements_editing <- F2_P1_P2_dataset %>% 
  rename(FishID = TubeLabel) %>% # for a later join
  mutate(FishID = as.character(FishID)) %>% # convert this column in a character variable
  select(FishID, Cross.Type, Cross)  

# DF containing info for fish with mouth open/closed for later correction
mouth_open_close_naming <- linear_measurement_dataset %>% 
  mutate(FishID = as.character(FishID)) %>% # this was done because later we will left_join with a character variable.
  select(Photo.Mouth.Open, FishID) # adding Photo.Mouth.Open column that we will left join later.

# complete linear measurement df which has many NAs for meta-data
edited_linear_measurements_dataset <- linear_measurement_dataset %>% 
  fill(FishID) %>% # fill empty cells with corresponding fish ID
  fill(Cross.Type) %>% # fill empty cells wth corresponding Cross.Type
  mutate(FishID = as.character(FishID)) %>% # converting column into character column
  select(Cross.Type, FishID, number, trait, measurement,Photo.Mouth.Open) %>% 
  left_join(linear_msr_info_df, by = "number") %>% # left joining by "number" becaue it is a common ID between table and dataset
  mutate(trait.y = fct_explicit_na(trait.y)) %>% # this had some implicit NAs that needed to to explicitly placed
  group_by(FishID,trait.y, Cross.Type) %>% 
  summarize(measurement = mean(measurement)) %>%  # duplicated data is averaged in measurement column 
  spread(key = trait.y, value = measurement) %>% # now each column represents a trait
  left_join(linear_measurements_editing, by = "FishID") %>% # left joining by common ID "FishID"
  left_join(mouth_open_close_naming, by = "FishID") # left joining by common ID "FishID"

# this gives the list of all FishIDs that are in edited_linear_measurements_dataset only
# (Because some linear measurements were accidentally made twice)
fish_Ids <- levels(factor(edited_linear_measurements_dataset$FishID))  # vector of fish Ids

# Create dataset that only has FishIds and has all the trait values for those particular fish that are in fishids 
all_fish_incomplete_1 <- F2_P1_P2_dataset %>% 
  filter(TubeLabel%in% fish_Ids) %>% # filters the fish ids from the column with all fish (fish from above 'edited linear measurements dataset') )
  select(Cross.Type,TubeLabel, N.plates, Dorsal.Ray:GillRakerL, fin.measurement) %>% 
  rename(FishID = TubeLabel) %>% # same column names are needed to merge all datasets
  mutate(FishID = as.character(FishID)) %>% # to left join with a character variable
  left_join(., edited_linear_measurements_dataset, by = "FishID") %>% # joining with another dataset by common identity
  select(-Cross.Type.x, -Cross.Type.y)  %>% 
  mutate(Photo.Mouth.Open = as_factor(Photo.Mouth.Open))

# deleting rows with fishids and fixing F1 dataset
F1_P2_dataset <- F1_P1_P2_dataset %>% 
  select(Cross.Type,FishID, Standard.Lt:Head.Lt,Body.Dpt:Second.Dor.Spine.Lt,Gape.Wt:Pelvic.Spine.Lt,N.plates:Gill.raker.lt_um_40x_mag, eye.measurements.2,pectoral.fin.dissection.measurement,Cross,Photo.Mouth.Open) %>% 
  rename(Dorsal.Ray = Dorsal.Fin.Rays, # renaming to make sure that all names in all datasets are same for later bind_rows
         Anal.Ray = Anal.Fin.Rays, 
         GillRakerNo = Gill.Raker.Ct.OuterArch.Lower,
         GillRakerL = Gill.raker.lt_um_40x_mag,
         Eye.Diam = eye.measurements.2,
         fin.measurement = pectoral.fin.dissection.measurement) %>% 
  mutate(Photo.Mouth.Open = as_factor(Photo.Mouth.Open))

# complete dataset has all the fish with all the traits (except additional_new_measurements)
complete_dataset <- F2_P1_P2_dataset %>% 
  select(Cross.Type, TubeLabel, Standard.Lt:Second.Dor.Spine.Lt, 
                Gape.Wt:Pelvic.Spine.Lt, N.plates, Dorsal.Ray:GillRakerL,
                fin.measurement, Cross, Photo.Mouth.Open) %>% 
  filter(!TubeLabel %in% fish_Ids) %>% # to avoid duplicating fish data, I am getting rid of some fish 
  rename(FishID = TubeLabel) %>% 
  rename(GillRakerNo = GillRakerNo) %>% 
  bind_rows(all_fish_incomplete_1, F1_P2_dataset)

#adding in additional_measurements_dataset
additional_new_measurements <- additional_measurements_dataset %>% 
  filter(Got.pics.or.not %in% c("yes", "done", "yep")) %>% # these are fish for which measurements were taken
  select(-fin_tube_KT, -Got.pics.or.not, -RackNo, -RackPos,-Photo_date,-category) %>% 
  mutate(Photo.Mouth.Open = as_factor(Photo.Mouth.Open))

new_fish_Ids <- levels(factor(additional_new_measurements$FishID)) # fish ids that are in additonal_new_measurements data

# df containing all fish with all trait measurements
new_complete_dataset <- complete_dataset %>% 
  filter(!FishID %in% new_fish_Ids) %>% # to avoid duplicating fish data, I am getting rid of some fish in here
  bind_rows(additional_new_measurements) # joining datasets

# assigning categories like P1,P2,F1,F2 to crosses and populations. P1 = anadromous, P2 = freshwater, F1 and F2 hybrids
editing_complete_dataset_1 <- new_complete_dataset %>%
  mutate(category = ifelse(Cross.Type %in% c("LCR-LCR"), "P1", # P1 are anadromous fish 
                           ifelse(Cross.Type %in% f1_CrossName_vec, "F1", 
                                  ifelse(Cross.Type %in% p2_CrossName_vec, "P2", 
                                         ifelse(Cross.Type %in% f2_CrossName_vec, "F2", NA)))))

# * Fixing unique errors ####
# KT_Gac_1809 was measured twice, so we took the mean of both measurements
fixing_dataset_mean_fish_1809 <- editing_complete_dataset_1 %>% 
  filter(FishID == "KT_Gac_1809") %>% # only have 'KT_Gac_1809' fish
  summarise_at(., vars(Standard.Lt:fin.measurement), mean) %>%  # take average of multiple measurements of KT_Gac_1809
  mutate(category = "F2") %>% # this is filling in the fish information
  mutate(Photo.Mouth.Open = "0") %>% 
  mutate(Cross = "LCR-BUL_F2.1") %>% 
  mutate(Cross.Type = "LCR-BUL_F2") %>% 
  mutate(FishID = "KT_Gac_1809")

# *convert all zero linear measurements to NA #### 
editing_complete_dataset <- editing_complete_dataset_1 %>%
  filter(!FishID %in% fish_to_exclude) %>% 
  # mutate(First.Dor.Spine.Lt = ifelse(FishID == "BUL-BUL4-4" & First.Dor.Spine.Lt == "0.016", "0", First.Dor.Spine.Lt))
  mutate(First.Dor.Spine.Lt = as.numeric(First.Dor.Spine.Lt)) %>% #june 5 edit
  mutate(Pelvic.Girdle.Lt = ifelse(Pelvic.Girdle.Lt == 0, NA, Pelvic.Girdle.Lt)) %>% #  replace 0 with NA; we later replace with pre-determined val
  mutate(Pelvic.Spine.Lt = ifelse(Pelvic.Spine.Lt == 0, NA, Pelvic.Spine.Lt)) %>% #  replace 0 with NA; we later replace with pre-determined val
  mutate(First.Dor.Spine.Lt = ifelse(First.Dor.Spine.Lt == 0, NA, First.Dor.Spine.Lt)) %>% #  replace 0 with NA; we later replace with pre-determined val
  mutate(Second.Dor.Spine.Lt = ifelse(Second.Dor.Spine.Lt == 0, NA, Second.Dor.Spine.Lt)) %>%  #  replace 0 with NA; we later replace with pre-determined val
  mutate(Dorsal.Ray = ifelse(FishID == "KT_Gac_3372" & Dorsal.Ray == 19, 10, Dorsal.Ray)) %>% # correcting measurement error
  mutate(Dorsal.Ray = as.integer(Dorsal.Ray)) %>% 
  mutate(Head.Lt = ifelse(FishID == "KT_Gac_3256" & Head.Lt == 6.791, 9.36, Head.Lt)) %>% # correcting measurement error 
  mutate(Head.Lt = ifelse(FishID == "PAQ-PAQ1-10" & Head.Lt == 7.047, 11.21, Head.Lt)) %>% # correcting measurement error
  mutate(First.Dor.Spine.Lt = ifelse(FishID == "BUL-BUL4-4" & First.Dor.Spine.Lt == 0.016, replace_0_val, First.Dor.Spine.Lt)) %>% # correcting measurement error
  filter(!Photo.Mouth.Open == "NA") %>% 
  filter(FishID != "BLANK") %>% 
  filter(FishID != "blank") %>% 
  filter(FishID != "Blank") %>%
  filter(FishID != "KT_Gac_201") %>% # there was a measurement error with second dorsal spine, but i fixed it and added in "df_201"
  filter(FishID != "KT_Gac_1809") %>% # duplicates in this dataset, but i fixed it and added in "df_201"
  filter(FishID != "KT_Gac_218") %>% # there was a measurement error with second dorsal spine, but i fixed it and added in "df_201"
  unique() %>% 
  filter(category %in% c("P1","P2","F1","F2")) %>% 
  # next few lines convert NA values to zero for traits that are 'measured but missing'
  mutate(Pelvic.Girdle.Lt = ifelse(is.na(Body.Wt) == F, replace_na(Pelvic.Girdle.Lt, replace_0_val), Pelvic.Girdle.Lt)) %>% 
  mutate(Pelvic.Spine.Lt = ifelse(is.na(Body.Wt) == F, replace_na(Pelvic.Spine.Lt, replace_0_val), Pelvic.Spine.Lt)) %>% 
  mutate(First.Dor.Spine.Lt = ifelse(is.na(Standard.Lt) == F, replace_na(First.Dor.Spine.Lt, replace_0_val), First.Dor.Spine.Lt)) %>%  
  mutate(Second.Dor.Spine.Lt = ifelse(is.na(Standard.Lt) == F, replace_na(Second.Dor.Spine.Lt, replace_0_val), Second.Dor.Spine.Lt)) 

# SIZE CORRECTION ####
# * log-log SL regressions
log_transformed_size_corrected_data <- editing_complete_dataset %>%
  mutate(size_correction_snout.Lt.log.res = resid(lm(log(Snout.Lt) ~ log(Standard.Lt),na.action = na.exclude)), 
         size_correction_Head.Lt.log.res = resid(lm(log(Head.Lt) ~ log(Standard.Lt), na.action = na.exclude)), 
         size_correction_Eye.Diam.log.res = resid(lm(log(Eye.Diam) ~ log(Standard.Lt),na.action = na.exclude)), 
         size_correction_Body.Dpt.log.res = resid(lm(log(Body.Dpt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_First.Dor.Spine.Lt.log.res = resid(lm(log((First.Dor.Spine.Lt)) ~ log((Standard.Lt)), na.action = na.exclude)),
         size_correction_Second.Dor.Spine.Lt.log.res = resid(lm(log(Second.Dor.Spine.Lt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_Gape.Wt.log.res = resid(lm(log(Gape.Wt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_Body.Wt.log.res = resid(lm(log(Body.Wt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_Pelvic.Girdle.Lt.log.res = resid(lm(log(Pelvic.Girdle.Lt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_Pelvic.Spine.Lt.log.res = resid(lm(log(Pelvic.Spine.Lt) ~ log(Standard.Lt), na.action = na.exclude)),
         size_correction_Gill_Raker_Lt.log.res = resid(lm(log(GillRakerL) ~ log(Standard.Lt), na.action = na.exclude)), 
         size_correction_pectoral.fin.lt.res = resid(lm(log(fin.measurement) ~ log(Standard.Lt), na.action = na.exclude)),
         gape_wt_mouth_coeff_1_log = (lm(size_correction_Gape.Wt.log.res ~ Photo.Mouth.Open))$coefficients[2], 
         gape_wt_mouth_coeff_2_log = (lm(size_correction_Gape.Wt.log.res ~ Photo.Mouth.Open))$coefficients[3],
         gape_wt_mouth_coeff_0_log = (lm(size_correction_Gape.Wt.log.res ~ Photo.Mouth.Open))$coefficients[1],
         snout_lt_coeff_2_log = (lm(size_correction_snout.Lt.log.res ~ Photo.Mouth.Open))$coefficients[3],
         snout_lt_coeff_1_log = (lm(size_correction_snout.Lt.log.res ~ Photo.Mouth.Open))$coefficients[2],
         snout_lt_coeff_0_log = (lm(size_correction_snout.Lt.log.res ~ Photo.Mouth.Open))$coefficients[1], 
         head_lt_coeff_2_log = (lm(size_correction_Head.Lt.log.res ~ Photo.Mouth.Open))$coefficients[3],
         head_lt_coeff_1_log = (lm(size_correction_Head.Lt.log.res ~ Photo.Mouth.Open))$coefficients[2],
         head_lt_coeff_0_log = (lm(size_correction_Head.Lt.log.res ~ Photo.Mouth.Open))$coefficients[1])



# In this data everything has been size corrected####
mouth_correction_data_log<- log_transformed_size_corrected_data %>% 
  mutate(Gape.Wt.correct.log=ifelse(Photo.Mouth.Open == "1", size_correction_Gape.Wt.log.res - gape_wt_mouth_coeff_1_log, 
                                    ifelse(Photo.Mouth.Open == "2",size_correction_Gape.Wt.log.res - gape_wt_mouth_coeff_2_log,
                                           ifelse(Photo.Mouth.Open == "0",size_correction_Gape.Wt.log.res - gape_wt_mouth_coeff_0_log, NA)))) %>%
  mutate(snout.lt.correct.log=ifelse(Photo.Mouth.Open == "1", size_correction_snout.Lt.log.res - snout_lt_coeff_1_log, 
                                     ifelse(Photo.Mouth.Open == "2",size_correction_snout.Lt.log.res- snout_lt_coeff_2_log,
                                            ifelse(Photo.Mouth.Open == "0",size_correction_snout.Lt.log.res - snout_lt_coeff_0_log, NA)))) %>% 
  mutate(head.lt.correct.log=ifelse(Photo.Mouth.Open == "1", size_correction_Head.Lt.log.res - head_lt_coeff_1_log, 
                                    ifelse(Photo.Mouth.Open == "2",size_correction_Head.Lt.log.res - head_lt_coeff_2_log,
                                           ifelse(Photo.Mouth.Open == "0",size_correction_Head.Lt.log.res - head_lt_coeff_0_log, NA)))) %>% 
  distinct() %>%
  dplyr::select(-"gape_wt_mouth_coeff_1_log",-"gape_wt_mouth_coeff_2_log",-"gape_wt_mouth_coeff_0_log",-"snout_lt_coeff_2_log",-"snout_lt_coeff_1_log",-"snout_lt_coeff_0_log",-"head_lt_coeff_2_log", -"head_lt_coeff_1_log", -"head_lt_coeff_0_log")

# fish that should be retained based on number of traits measured
fish_with_all_msrments <- mouth_correction_data_log %>% 
  select(FishID, Snout.Lt:fin.measurement) %>% 
  pivot_longer(Snout.Lt:fin.measurement, names_to = "trait", values_to = "value") %>% 
  group_by(FishID) %>% 
  na.omit() %>% 
  summarise(n_traits_msrd = n()) %>% 
  filter(n_traits_msrd == 16)
    
trait_names <- mouth_correction_data_log

# create dataset for later joining of FW_pop
Cross_FW_pop_df <- mouth_correction_data_log %>% 
  mutate(FW_pop = Cross.Type) %>% 
  separate(FW_pop, into = c("first", "second"), sep = "-") %>% # Separate population name into two components. For example, separate LCR-LQB into "LCR" in 'first' column and "LQB" in second column
  mutate(second = ifelse(second == "LCR", NA, second)) %>% # if second column has 'LCR', replace it with NA, otherwise keep it the same
  mutate(first = ifelse(first == "LCR", NA, first)) %>% # if first column has 'LCR', replace it with NA otherwise keep it the same
  rename(FW_pop = second) %>% 
  ungroup() %>% 
  select(Cross.Type, FW_pop) %>% 
  unique() %>% 
  separate(FW_pop, into = c("FW_pop", "delete_me", sep = "_")) %>% 
  select(Cross.Type, FW_pop) %>% 
  mutate(FW_pop = ifelse(Cross.Type == "LCR-LCR", "LCR", FW_pop))

category_df <- mouth_correction_data_log %>% 
  select(FishID, category) %>% 
  unique()

# Changing trait names in the column trait in order to make them more easy to read and identify. 
trait_names <- c(size_correction_Eye.Diam.log.res = "eye diameter", 
                 size_correction_Body.Dpt.log.res = "body depth",
                 size_correction_First.Dor.Spine.Lt.log.res = "first dorsal spine lt",
                 size_correction_Second.Dor.Spine.Lt.log.res = "second dorsal spine lt",
                 size_correction_Body.Wt.log.res = "body width",
                 size_correction_Pelvic.Girdle.Lt.log.res = "pelvic girdle lt",
                 size_correction_Pelvic.Spine.Lt.log.res = "pelvic spine length",
                 size_correction_Gill_Raker_Lt.log.res = "gill raker lt",
                 Gape.Wt.correct.log = "gape width",
                 snout.lt.correct.log = "snout length",
                 head.lt.correct.log = "head length",
                 N.plates = "no. of plates",
                 Dorsal.Ray = "dorsal ray",
                 Anal.Ray = "anal ray",
                 Standard.Lt = "standard length",
                 size_correction_pectoral.fin.lt.res = "pectoral fin length",
                 GillRakerNo = "gill raker no.")

# * standardize traits ####
standardized_trait_df_for_pairwise <- mouth_correction_data_log %>% 
  left_join(Cross_FW_pop_df, by = "Cross.Type") %>% 
  select(FishID, FW_pop, Cross.Type, Cross, names(trait_names)) %>% 
  mutate(Parent_Hybrid = ifelse(Cross.Type %in% c(p2_CrossName_vec, "LCR-LCR"), "Parent", "Hybrid")) %>% # if parent; call parent; else hybrid
  pivot_longer(cols = names(trait_names), names_to = "trait", values_to = "raw") %>% 
  # na.omit() %>% 
  group_by(trait) %>% 
  mutate(scaled_val = scale(raw)) %>%  # standardize traits to mean 0 sd 1
  left_join(category_df, by = "FishID") %>% # bring in 'category' variable
  filter(trait != "Standard.Lt") %>% 
  filter(trait != "size_correction_Eye.Diam.log.res")  # drop Eye because it's not repeatable

# * standardize traits ####
standardized_trait_df <- mouth_correction_data_log %>% 
  left_join(Cross_FW_pop_df, by = "Cross.Type") %>% 
  select(FishID, FW_pop, Cross.Type, Cross, names(trait_names)) %>% 
  mutate(Parent_Hybrid = ifelse(Cross.Type %in% c(p2_CrossName_vec, "LCR-LCR"), "Parent", "Hybrid")) %>% # if parent; call parent; else hybrid
  pivot_longer(cols = names(trait_names), names_to = "trait", values_to = "raw") %>% 
  na.omit() %>% 
  group_by(trait) %>% 
  mutate(scaled_val = scale(raw)) %>%  # standardize traits to mean 0 sd 1
  left_join(category_df, by = "FishID") %>% # bring in 'category' variable
  filter(trait != "Standard.Lt") %>% 
  filter(trait != "size_correction_Eye.Diam.log.res") %>%  # drop Eye because it's not repeatable
  filter(!FishID %in% fish_to_exclude)

# * parent means df ####
# Means of all 'pure' parents
stickle_scaled_parent_means <- standardized_trait_df %>%
  filter(Parent_Hybrid == "Parent") %>% # only include parents
  # giving each family (cross) and population an equal weighting when calculating means 
  group_by(FW_pop, trait, Cross) %>% 
  summarise(mean_scaled_val = mean(scaled_val, na.rm = T)) %>% #average calculated by grouping by freshwater population, PC, cross, and Parent ID
  group_by(trait, FW_pop) %>% 
  summarise(mean_scaled_val = mean(mean_scaled_val, na.rm = T)) %>% # average calculated by grouping by freshwater population, PC, and Parent ID  
  pivot_wider(names_from = trait, values_from = mean_scaled_val)

# LCR (anadromous) parent means 
LCR_scaled_parent_means.df <- stickle_scaled_parent_means %>% 
  filter(FW_pop == "LCR") %>% # choose only Parent ID == A (i.e., anadromous parent)
  gather(Anal.Ray:snout.lt.correct.log, key = 'trait', value = 'value') %>% # gather all 16 PCs in the column called 'trait' and all values in the column called 'value'
  ungroup() %>% 
  select(-FW_pop) %>%  # remove FW_pop column
  rename(P1 = value)

# which fish have all traits measured? (used later)
all_traits_measured_names <- standardized_trait_df %>% 
  select(FishID, trait, raw) %>% 
  pivot_wider(names_from = trait, values_from = raw) %>% 
  na.omit()

# Subset data to only hybrids
hybrid_scaled_traits_ind.df <- standardized_trait_df %>% 
  filter(Parent_Hybrid == "Hybrid") %>% # only include hybrids
  select(-raw, -Parent_Hybrid) %>% # remove the columns "ParentID" and "Parent_Hybrid"
  rename(scaled_hybrid_value = scaled_val) 

# scales parent values (used later on!)
fw_parent_dataset_for_p_mismatch <- standardized_trait_df %>% 
  filter(trait != "Standard.Lt") %>% 
  filter(category == "P2") %>% # only include hybrids
  select(-raw, -Parent_Hybrid) %>% # remove the columns "ParentID" and "Parent_Hybrid"
  rename(scaled_hybrid_value = scaled_val) 

# df with information for each fish
FishID_info <- standardized_trait_df %>% 
  ungroup() %>% 
  select(FishID:Cross, category) %>% 
  unique()

# Generate parent mean trait vectors ####
stickle_parent_pheno_data <- standardized_trait_df %>% 
  filter(category %in% c("P1", "P2")) %>% # retain only parents; P1 = LCR (anadromous); P2 = Freshwater parent
  group_by(FW_pop, Cross, trait) %>% # group by population   # removed cross
  summarise(mean_value = mean(scaled_val)) %>% # mean value of each trait
  group_by(FW_pop, trait) %>% 
  summarise(mean_value = mean(mean_value)) %>%   # mean value of each trait
  # now need it wide
  pivot_wider(names_from = trait, values_from = mean_value)

# trait means for FW pops
stickle_FW_parent_pheno_data <- stickle_parent_pheno_data %>% 
  pivot_longer(cols = Anal.Ray:snout.lt.correct.log, names_to = "trait", values_to = "value") %>% 
  rename(FW_scaled_mean = value)

FW_scaled_parent_mean.df <- stickle_parent_pheno_data %>% 
  pivot_longer(cols = Anal.Ray:snout.lt.correct.log, names_to = "trait", values_to = "value") %>% 
  rename(P2 = value)

# trait means for LCR
LCR_dat <- stickle_parent_pheno_data %>% 
  filter(FW_pop == "LCR") %>% 
  pivot_longer(cols = Anal.Ray:snout.lt.correct.log, names_to = "trait", values_to = "value") # LCR 15 traits; long

LCR_trait_mean_vec <- LCR_dat$value # get a vector of the mean LCR traits 

# DIST BTT PARENTS ####
parent_dists_df <- stickle_parent_pheno_data %>% 
  rowwise() %>% # by rows
  filter(FW_pop != "LCR") %>% # exclude LCR
  do({
    result = as_tibble(.) # making it a dataframe
    result$p_dist = parent_euclidean(result, which_LCR = LCR_trait_mean_vec)# p_dist is parent euclidean distance
    result
  }) %>% 
  select(FW_pop, p_dist) %>% # only include freshwater population name and parent eucliden distance
  rename(parent_trait_divergence = p_dist)

# Individual-level parent dist
stickle_parent_INDIVIDUAL_pheno_data <- standardized_trait_df %>% 
  filter(category %in% c("P1", "P2")) %>% # retain only parents; P1 = LCR (anadromous); P2 = Freshwater parent
  select(FishID, trait, scaled_val) %>%
  na.omit() %>% 
  # now need it wide
  pivot_wider(names_from = trait, values_from = scaled_val)

stickle_parent_FAMILY_pheno_data <- standardized_trait_df %>% 
  filter(category %in% c("P1", "P2")) %>% # retain only parents; P1 = LCR (anadromous); P2 = Freshwater parent
  group_by(FW_pop, Cross, trait) %>% 
  summarise(scaled_val = mean(scaled_val)) %>% 
  ungroup() %>% 
  select(Cross, trait, scaled_val) %>%
  na.omit() %>% 
  # now need it wide
  pivot_wider(names_from = trait, values_from = scaled_val)


# Distance of each individual FW parent individual to the mean LCR phenotype
ind_parent_dists_df <- stickle_parent_INDIVIDUAL_pheno_data %>% 
  lect(FishID, p_dist) %>% # only include freshwater population name and parent eucliden distance
  rename(parent_trait_divergence = p_dist)  

# Distance of each individual FW parent individual to the mean LCR phenotype
fam_parent_dists_df <- stickle_parent_FAMILY_pheno_data %>%
  rowwise() %>% # by rows
  # should purrr this
  do({
    result = as_tibble(.) # making it a dataframe
    result$p_dist = parent_euclidean(result, which_LCR = LCR_trait_mean_vec)# p_dist is parent euclidean distance
    result
  }) %>% 
  select(Cross, p_dist) %>% # only include freshwater population name and parent eucliden distance
  rename(parent_trait_divergence = p_dist)

family_parent_dists_df_w_info <- fam_parent_dists_df %>% 
  left_join(., FishID_info, by = "Cross") %>% 
  select(-FishID) %>% 
  unique() %>% 
  mutate(parent_trait_divergence = as.numeric(parent_trait_divergence)) 

# write final file of phenotype data
write_csv(family_parent_dists_df_w_info, "data/derived_data/family_parent_dists_df_w_info.csv")
write_csv(parent_dists_df, "data/derived_data/parent_dists_df.csv")
write_csv(LCR_scaled_parent_means.df, "data/derived_data/LCR_scaled_parent_means.df.csv")
write_csv(standardized_trait_df, "data/derived_data/standardized_trait_df.csv")
write_csv(mouth_correction_data_log,"data/derived_data/size_mouth_corrected_pheno_data.csv")
write_csv(category_df, 'data/derived_data/category_df.csv')
write_csv(hybrid_scaled_traits_ind.df, 'data/derived_data/hybrid_scaled_traits_ind.df.csv')

