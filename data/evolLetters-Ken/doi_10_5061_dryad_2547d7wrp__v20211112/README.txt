README for "Adaptive divergence and the evolution of hybrid trait mismatch in threespine stickleback"
Written by: Ken A. Thompson (ken.thompson@zoology.ubc.ca; kthomp1063@gmail.com)

2021-06-26

The main project directory contains two folders: 'data' and 'scripts'. It also contains an 'R Project' file that interfaces with RStudio but is not neccessary.  

# SCRIPTS FODLER #
The 'scripts' folder contains the analysis scripts used to generate the results of the study. The scripts should be fully reproducible and sufficiently annotated for the reviewers to reproduce our findings. Each script works on its own, but ideally the scripts could be run in order.

1_data_cleaning.R---takes raw data and generates data for downstream analysis
2_trait_mismatch.R---calculates mismatch for each individual and each pair of traits
3_patterns_dominance_variation.R---calculates dominance and trait variation, and their effects on mismatch
4_reapeatability.R---conducts the repeatability analysis.
5_F2_correlation_matrix.R---generates the correlation matrices using F2 hybrid data.
6_models_and_analysis.R---runs statistical models presented in the MS.
7_supplementary_analyses.R---analyses concerning statements of 'see online R script' etc. not in paper, but here just in case.
8_figures.R---generates figures that are in the paper.

# DATA FOLDER #
The 'data' folder contains two folders: 'raw_data' and 'derived_data'.

The 'raw_data' folder contains all of the datashets into which raw measurements were input. The formats vary a bit among dataframes but the raw data are preserved as they are. These are raw dataframes, warts and all! They are preserved in their entirety in the event that we need to go back to double check something that seemed unimportant at the time. The 'repeatability' subfolder contains measurements that I (Ken) made on all fish. Dataset descriptions are below. A guide to our measurements is included in the manuscript as a figure, as are all methods. A complete dataset of proceessed data is available in the 'derived_data' folder (described below). Users interested in a dataset of size-corrected and organized measurements should look at the 'derived_data'---the raw data will probably be of little use to you!

# raw_data/raw_additional_linear_msr.csv
Category: Type of fish (P1 = LCR, P2 = FW parent; F1/F2: Hybrids)
Cross.Type: Unique cross type
FishID: Unique fish label; we fill 'down' in R; matches paper label in tube
number: measurement ID for a given fish
angle: additional output from ImageJ
measurement: length in 'mm'
Photo.Mouth.Open: score of whether fish's mouth is closed (0), partially open (0.5), or open (1)
trait: trait name
Cross: unique family

# raw_data/raw_msr_set_1
tempfishno: a unique number temporarily assigned to each fish during early stages of data processing
finfisihno: essentially 'row number'
Cross: unique family
Cross.Type: Unique cross type
Dam.Pop: source of cytoplasm
Sire.Pop: source of y (in males) 
Parent.or.Hybrid: whether the fish is pure or a hybrid
FishNo: fish number within family
FishID: Unique fish label; matches paper label in tube
TailBoxNo: Cryobox where tail tissue was kept
TailTubeNo: Position of tail tissue in cryobox
SampleDate: data that the measurements were made 
FishTubeRack: Unique rack (i.e, 25-tube stand) where tube is kept
FishTubeRack_Pos: Position in rack (1 to 25)
PhotoNote: Note re: photo
PhotoDate: Date photo was taken
Standard.Lt: Standard length (mm)
Snout.Lt: Snouth length (mm)
Head.Lt: Head length (mm)
Eye.Diam: Eye diameter (mm); note this was not used in analysis due to low repeatability
Body.Dpt: Bodth depth (mm)
First.Dor.Spine.Lt: Length of first dorsal spine (mm)
Second.Dor.Spine.Lt: Length of second dorsal spine (mm)
Pectoral.Fin.Lt: Length of pectoral fin (mm)
Gape.Wt: Width of gape (mm)
Body.Wt: Wdith of body (mm)
Pelvic.Girdle.Lt: Length of pelvic girdle (mm)
Pelvic.Spine.Lt: Length of pelvic spine (mm)
Photo.Mouth.Open: score of whether fish's mouth is closed (0), partially open (0.5), or open (1)
N.plates: Count of lateral plate number
Dorsal.Fin.Rays: Count of rays on the dorsal fin
Anal.Fin.Rays: Count of rays on the anal fin
Gill.Raker.Ct.OuterArch.Lower: Count of gill rakers
Gill.raker.lt_um_40x_mag: Length of longest gill raker (ocular units)
Raker.Note: comment on rakers
random.no.: number used at a step in data processing
random.no: number used at a step in data processing
eye.measurements.2: repeated measurement of the eye (we tried at one point to see if we could improve accuracy, and failed)
rough stuff(don't need): Additional output from imageJ
rough stuff(don't need)_1: Additional output from imageJ
pectoral.fin.dissection.measurement: Unique number for preserved pectoral fin
Pectoral.fin.number: unique fin number in rack
Pectoral.fin.rack.no: unique fin rack
right/left: was measurement made on left or right pectoral fin

# raw_data/raw_msr_set_2
TubeLabel: Unique fish ID; paper label
Cross: unique family
Cross.Type: Unique cross type
RackNo: Unique rack (i.e, 25-tube stand) where tube is kept
RackPos: Position in rack (1 to 25)
ID.No: Temporary processing variable
Fish_Unique_ID: Temporary processing variable
Note: a note
Standard.Lt: Standard length (mm)
Snout.Lt: Snouth length (mm)
Head.Lt: Head length (mm)
Eye.Diam: Eye diameter (mm); note this was not used in analysis due to low repeatability
Body.Dpt: Bodth depth (mm)
First.Dor.Spine.Lt: Length of first dorsal spine (mm)
Second.Dor.Spine.Lt: Length of second dorsal spine (mm)
Pectoral.Fin.Lt: placeholder (replaced with 'fin.measurement' in the script)
Gape.Wt: Width of gape (mm)
Body.Wt: Wdith of body (mm)
Pelvic.Girdle.Lt: Length of pelvic girdle (mm)
Pelvic.Spine.Lt: Length of pelvic spine (mm)
Photo.Mouth.Open: score of whether fish's mouth is closed (0), partially open (0.5), or open (1)
Photo_date: Date photo was taken
N.plates: Count of lateral plate number
Dorsal.Ray: Count of rays on the dorsal fin
Anal.Ray: Count of rays on the anal fin
GillRakerNo: Count of gill rakers
GillRakerL: Length of longest gill raker (ocular units)
Raker.Note: comment on rakers
fin.no: Unique number for preserved pectoral fin
fin.rack.no: unique fin rack
random.no: random number used during data processing
random.no_1: random number used during data processing
fin.measurement: Length of pectoral fin (mm)
X34: silly column of all NAs added by R
X35: silly column of all NAs added by R
X36: silly column of all NAs added by R
X37: silly column of all NAs added by R
X38: silly column of all NAs added by R
X39: silly column of all NAs added by R
X40: silly column of all NAs added by R
X41: silly column of all NAs added by R

# raw_data/raw_msr_set_3
Cross.Type: Unique cross type
FishID: Unique fish ID; paper label
Standard.Lt: Standard length (mm)
Snout.Lt: Snouth length (mm)
Head.Lt: Head length (mm)
Eye.Diam: Eye diameter (mm); note this was not used in analysis due to low repeatability
Body.Dpt: Bodth depth (mm)
First.Dor.Spine.Lt: Length of first dorsal spine (mm)
Second.Dor.Spine.Lt: Length of second dorsal spine (mm)
Gape.Wt: Width of gape (mm)
Body.Wt: Wdith of body (mm)
Pelvic.Girdle.Lt: Length of pelvic girdle (mm)
Pelvic.Spine.Lt: Length of pelvic spine (mm)
N.plates: Count of lateral plate number
Dorsal.Ray: Count of rays on the dorsal fin
Anal.Ray: Count of rays on the anal fin
GillRakerNo: Count of gill rakers
GillRakerL: Length of longest gill raker (ocular units)
fin.measurement: Length of pectoral fin (mm)
Cross: unique family
Photo.Mouth.Open: score of whether fish's mouth is closed (0), partially open (0.5), or open (1)
category: Type of fish (P1 = LCR, P2 = FW parent; F1/F2: Hybrids)
fin_tube_KT: unique fin tube
RackNo: Unique rack (i.e, 25-tube stand) where tube is kept
RackPos: Position in rack (1 to 25)
Photo_date: Date photo was taken
Got pics or not: whether a fish was photographed

# raw_data/repeatability/ken_repeatability_fins.csv
FishID: Unique fish ID; paper label
msr_order: order that the photos were taken
length_mm: length in mm
rep: measurement 1 or measurement 2

# raw_data/repeatability/ken_repeatability_fish_2.csv
fish: Unique fish ID; paper label 
Standard.Lt: Standard length (mm)
Snout.Lt: Snouth length (mm)
Head.Lt: Head length (mm)
Eye.Diam: Eye diameter (mm); note this was not used in analysis due to low repeatability
Body.Dpt: Bodth depth (mm)
First.Dor.Spine.Lt: Length of first dorsal spine (mm)
Second.Dor.Spine.Lt: Length of second dorsal spine (mm)
Gape.Wt: Width of gape (mm)
Body.Wt: Wdith of body (mm)
Pelvic.Girdle.Lt: Length of pelvic girdle (mm)
Pelvic.Spine.Lt: Length of pelvic spine (mm)

# raw_data/repeatability/ken_repeatability_scope.csv
fish: Unique fish ID; paper label
N.plates: Count of lateral plate number
Dorsal.Fin.Rays: Count of rays on the dorsal fin
Anal.Fin.Rays: Count of rays on the anal fin
Gill.Raker.Ct.OuterArch.Lower: Count of gill rakers
Gill.raker.lt_um_40x_mag: Length of longest gill raker (ocular units)
