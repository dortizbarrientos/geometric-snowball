#### Calculating genetic distances between pops on the phylogeny ####

# Load the ape package
library(ape)

# Read the Nexus file
tree <- read.nexus("~Phylogeny.nex")
# Calculate the distance matrix of the tree
distance_matrix <- cophenetic(tree)
heatmap(as.matrix(distance_matrix))
# Calculate distance matrix for each parapatric pair. e.g.,
distance_matrix["D32", "H12"]

#### Environmental associations ####

# Load the data table
CrossesDH <- read.delim("~CrossesDH-ReciprocalsAveraged.txt", header=T)
# Subset the data into DD and HH
DDdata <- CrossesDH[CrossesDH$CrossType.1 == 'DD',]
HHdata <- CrossesDH[CrossesDH$CrossType.1 == 'HH',]

# Dunes
# SeedSet

# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIseedSet~Clade+EnvDist.MDS.5PCs, data = DDdata))

summary(lm(RIseedSet~Clade+EnvDist.MDS.5PCs, data = DDdata))

model2 <- lm(RIseedSet~Clade+EnvDist.MDS.5PCs, data = DDdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(DDdata$RIseedSet ~ DDdata$EnvDist.MDS.5PCs)

# Viability
# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIviability~Clade+EnvDist.MDS.5PCs, data = DDdata))

summary(lm(RIviability~Clade+EnvDist.MDS.5PCs, data = DDdata))

model2 <- lm(RIviability~Clade+EnvDist.MDS.5PCs, data = DDdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(DDdata$RIviability ~ DDdata$EnvDist.MDS.5PCs)


# Headlands
# SeedSet
# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIseedSet~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata))

summary(lm(RIseedSet~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata))

model2 <- lm(RIseedSet~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(HHdata$RIseedSet ~ HHdata$EnvDist.MDS.5PCs)

# Viability
# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIviability~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata))

summary(lm(RIviability~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata))

model2 <- lm(RIviability~Genetic.distance+EnvDist.MDS.5PCs, data = HHdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(HHdata$RIviability ~ HHdata$EnvDist.MDS.5PCs)



#### Phenotypic associations ####

# Headlands
# SeedSet
# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIseedSet~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata))

summary(lm(RIseedSet~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata))

model2 <- lm(RIseedSet~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(HHdata$RIseedSet ~ HHdata$PhenoDist.MDS.7PCs)

# Viability
# Association between ecological divergence and reproductive isolation, whilst correcting for divergence time
anova(lm(RIviability~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata))

summary(lm(RIviability~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata))

model2 <- lm(RIviability~Genetic.distance+PhenoDist.MDS.7PCs, data = HHdata)
anova(model2)
beta <- model2$coefficients 
abline(beta[1], beta[2], col=1) 
for (i in 2:3) {  abline(beta[1]+beta[i],beta[2], col=i-1)   }	# add lines for other groups 

plot(HHdata$RIviability ~ HHdata$PhenoDist.MDS.7PCs)