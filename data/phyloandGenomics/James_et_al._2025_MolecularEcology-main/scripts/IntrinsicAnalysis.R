# Read input file
CrossesDH <- read.delim("~/CrossesDH.txt")
View(CrossesDH)

library(lmerTest)

#### Linear models to detect the effect of ecotype and clade on RI, with cross type as a random factor ####

# Seed set
# Run the model, with CrossType as a random factor
SeedSet <- lmer(RIseedSet.SobelChen ~ Clade * Ecotype + (1 | CrossType), data = CrossesDH)

#Do the anova of the model to get the P-values
anova(SeedSet)

#As the interaction is not significant, remove it from the model
SeedSet <- lmer(RIseedSet.SobelChen ~ Clade + Ecotype + (1 | CrossType), data = CrossesDH)
anova(SeedSet)
summary_model <- summary(SeedSet)

# Extract variance components
total_variance <- attr(VarCorr(SeedSet), "sc")^2
fixed_variance <- sum(diag(vcov(SeedSet)))

# Calculate approximate R-squared
r_squared <- fixed_variance / (fixed_variance + total_variance)

# Print approximate R-squared
print(r_squared)


# Viability

# Run the model, with CrossType as a random factor
Viability <- lmer(RIviability.SobelChen ~ Clade + Ecotype + (1 | CrossType), data = CrossesDH)

#Do the anova of the model to get the P-values
anova(Viability)
summary_model <- summary(Viability)

# Extract variance components
total_variance <- attr(VarCorr(Viability), "sc")^2
fixed_variance <- sum(diag(vcov(Viability)))

# Calculate approximate R-squared
r_squared <- fixed_variance / (fixed_variance + total_variance)

# Print approximate R-squared
print(r_squared)


#### Models with the ecotype crosses ####

# Seed set
# Run the model
SeedSet <- lm(RIseedSet.SobelChen ~ Clade + EcotypeCross, data = CrossesDH)
anova(SeedSet)
summary(SeedSet)

# Viability
# Run the model
Viability <- lm(RIviability.SobelChen ~ Clade + EcotypeCross, data = CrossesDH)
anova(Viability)
summary(Viability)


#### Figure 3A ####

clade_order <- c('Within', 'Between')

# Function to plot with SE bars
# Fix for the summarySE function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # We don't need rlang explicitly loaded
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column using the measurevar name
  # This is the part that needs to change
  names(datac)[names(datac) == "mean"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Change the "Same" values to "aSame" so they plot in the correct order
CrossesDH$Ecotype[CrossesDH$Ecotype == "Same"] <- "aSame"

# Using the function on the dataset for plotting
Figure3A <- summarySE(CrossesDH, measurevar="RIseedSet.SobelChen", groupvars=c("Clade", "Ecotype"))

clade_order <- c('Within', 'Between')

RISEEDSET <- ggplot(Figure3A, aes(fill=Ecotype, y=RIseedSet.SobelChen, x=Clade)) + 
  geom_bar(position="dodge", stat = "summary", fun = mean, colour= "black", linewidth=0.2) + 
  scale_x_discrete(limits = clade_order) +
  labs(y = "Reproductive isolation (F1 seed set)") + 
  coord_cartesian(ylim = c(-1, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7),
        axis.ticks.length=unit(.17, "cm"), axis.ticks = element_line(colour = "black", linewidth=0.4),
        panel.background = element_blank()) + 
  theme(legend.position = "none") + # Remove the legend
  theme(text=element_text(family="Helvetica Neue Light", size=17)) +  
  theme(axis.text.x = element_text(colour="black")) + 
  theme(axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept=0, color = "black", linewidth = 0.2)+ 
  theme(legend.position = "none") + # Ensure legend is removed
  theme(legend.text = element_text(size=12)) + 
  theme(legend.title = element_text(size=14)) +
  scale_fill_manual(values = c("#ADDFAD", "#43ac43"), name = "Ecotype comparison") +
  geom_errorbar(aes(ymin=RIseedSet.SobelChen-se, ymax=RIseedSet.SobelChen+se), width=.1, linewidth=0.4,
                position=position_dodge(0.9))  +
  theme(legend.spacing.y = unit(0.1, 'cm'))  +
  guides(fill = guide_legend(byrow = TRUE)) +
  geom_jitter(data = CrossesDH, aes(x = Clade, y = RIseedSet.SobelChen), 
              color = "#6c6c6c", size = 2, shape = 16, alpha = 0.7, 
              position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.9))
RISEEDSET

ppi <- 300
png("RIseedset.png", width=6*ppi, height=4.2*ppi, res=ppi)
plot(RISEEDSET)
dev.off()

#### Figure 3B ####
# Function to plot with SE bars
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column using the measurevar name
  # Fix: Use base R naming instead of tidyverse renaming
  names(datac)[names(datac) == "mean"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Change the "Same" values to "aSame" so they plot in the correct order
CrossesDH$Ecotype[CrossesDH$Ecotype == "Same"] <- "aSame"

# Using the function on the dataset for plotting
# Fix: Change RIviability to RIviability.SobelChen to match your new analysis
Figure3B <- summarySE(CrossesDH, measurevar="RIviability.SobelChen", groupvars=c("Clade", "Ecotype"))

RIVIABILITY <- ggplot(Figure3B, aes(fill=Ecotype, y=RIviability.SobelChen, x=Clade)) + 
  geom_bar(position="dodge", stat = "summary", fun = mean, colour= "black", linewidth=0.2) + 
  scale_x_discrete(limits = clade_order) +
  labs(y = "Reproductive isolation (F1 viability)") + 
  coord_cartesian(ylim = c(-1, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7),
        axis.ticks.length=unit(.17, "cm"), axis.ticks = element_line(colour = "black", linewidth=0.4),
        panel.background = element_blank()) + 
  theme(text=element_text(family="Helvetica Neue Light", size=17)) +  
  theme(axis.text.x = element_text(colour="black")) + 
  theme(axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept=0, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#ADDFAD", "#43ac43"), name = "Ecotype comparison") +
  geom_errorbar(aes(ymin=RIviability.SobelChen-se, ymax=RIviability.SobelChen+se), width=.1, linewidth=0.4,
                position=position_dodge(0.9))  +
  geom_jitter(data = CrossesDH, aes(x = Clade, y = RIviability.SobelChen), 
              color = "#6c6c6c", size = 2, shape = 16, alpha = 0.7, 
              position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.9)) + 
  theme(legend.position = "none") # Remove the legend

RIVIABILITY

ppi <- 300
png("RIviability.png", width=6*ppi, height=4.2*ppi, res=ppi)
plot(RIVIABILITY)
dev.off()

#### Figure 3C ####
CrossesDH$EcotypeCross[CrossesDH$EcotypeCross == "DD"] <- "a"
CrossesDH$EcotypeCross[CrossesDH$EcotypeCross == "HH"] <- "b"
CrossesDH$EcotypeCross[CrossesDH$EcotypeCross == "DH"] <- "c"
CrossesDH$EcotypeCross[CrossesDH$EcotypeCross == "HD"] <- "d"

Figure3C <- summarySE(CrossesDH, measurevar="RIseedSet.SobelChen", groupvars=c("Clade", "EcotypeCross"))

RISEEDSET <- ggplot(Figure3C, aes(fill=EcotypeCross, y=RIseedSet.SobelChen, x=Clade)) + 
  geom_bar(position="dodge", stat = "summary", fun = mean, colour= "black", linewidth=0.2) + 
  scale_x_discrete(limits = clade_order) +
  labs(y = "Reproductive isolation (F1 seed set)") + 
  coord_cartesian(ylim = c(-1, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7),
        axis.ticks.length=unit(.17, "cm"), axis.ticks = element_line(colour = "black", linewidth=0.4),
        panel.background = element_blank()) + 
  theme(legend.position = "none") +
  theme(text=element_text(family="Helvetica Neue Light", size=17)) +  
  theme(axis.text.x = element_text(colour="black")) + 
  theme(axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept=0, color = "black", linewidth = 0.2)+ 
  # Fix for legend position deprecation warning
  theme(legend.position.inside = c(0.13, 0.75))+
  theme(legend.text = element_text(size=12)) + 
  theme(legend.title = element_text(size=14)) +
  scale_fill_manual(values = c("#0072B1", "#D45E00", "#0072B180", "#D45E0080"), name = "Cross type") +
  geom_errorbar(aes(ymin=RIseedSet.SobelChen-se, ymax=RIseedSet.SobelChen+se), width=.1, linewidth=0.4,
                position=position_dodge(0.9))  +
  theme(legend.spacing.y = unit(0.1, 'cm'))  +
  guides(fill = guide_legend(byrow = TRUE)) +
  geom_jitter(data = CrossesDH, aes(x = Clade, y = RIseedSet.SobelChen), 
              color = "#6c6c6c", size = 2, shape = 16, alpha = 0.7, 
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9))

RISEEDSET

ppi <- 300
png("RIseedsetEcotypes.png", width=6*ppi, height=4.2*ppi, res=ppi)
plot(RISEEDSET)
dev.off()


#### Figure 3D ####

Figure3D <- summarySE(CrossesDH, measurevar="RIviability.SobelChen", groupvars=c("Clade", "EcotypeCross"))

RIVIABILITY <- ggplot(Figure3D, aes(fill=EcotypeCross, y=RIviability.SobelChen, x=Clade)) + 
  geom_bar(position="dodge", stat = "summary", fun = mean, colour= "black", linewidth=0.2) + 
  scale_x_discrete(limits = clade_order) +
  labs(y = "Reproductive isolation (F1 viability)") + 
  coord_cartesian(ylim = c(-1, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7),
        axis.ticks.length=unit(.17, "cm"), axis.ticks = element_line(colour = "black", linewidth=0.4),
        panel.background = element_blank()) + 
  theme(text=element_text(family="Helvetica Neue Light", size=17)) +  
  theme(axis.text.x = element_text(colour="black")) + 
  theme(axis.text.y = element_text(colour="black")) +
  geom_hline(yintercept=0, color = "black", linewidth = 0.2)+ 
  # Remove the legend completely
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#0072B1", "#D45E00", "#0072B180", "#D45E0080"), name = "Cross type") +
  geom_errorbar(aes(ymin=RIviability.SobelChen-se, ymax=RIviability.SobelChen+se), width=.1, linewidth=0.4,
                position=position_dodge(0.9)) +
  # These legend formatting options won't have any effect since legend is hidden
  # but keeping them doesn't hurt anything
  theme(legend.spacing.y = unit(0.1, 'cm'))  +
  guides(fill = guide_legend(byrow = TRUE)) +
  geom_jitter(data = CrossesDH, aes(x = Clade, y = RIviability.SobelChen), 
              color = "#6c6c6c", size = 2, shape = 16, alpha = 0.7, 
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9))

RIVIABILITY

ppi <- 300
png("RIviabilityEcotypes.png", width=6*ppi, height=4.2*ppi, res=ppi)
plot(RIVIABILITY)
dev.off()


#### Stats for differences from zero ####
CrossesDH <- read.delim("~/CrossesDH.txt")

# Subset the data into DD and HH
DDdata <- CrossesDH[CrossesDH$EcotypeCross == 'DD',]
HHdata <- CrossesDH[CrossesDH$EcotypeCross == 'HH',]
DDdataWithin <- DDdata[DDdata$Clade == 'Within',]
DDdataBetween <- DDdata[DDdata$Clade == 'Between',]
HHdataWithin <- HHdata[HHdata$Clade == 'Within',]
HHdataBetween <- HHdata[HHdata$Clade == 'Between',]

# T-tests to test whether each category (within and between clades for each ecotype) is different from 0
t.test(DDdataWithin$RIseedSet.SobelChen, mu = 0, alternative = "greater")
t.test(DDdataBetween$RIseedSet.SobelChen, mu = 0, alternative = "greater")
t.test(DDdataWithin$RIviability.SobelChen, mu = 0, alternative = "greater")
t.test(DDdataBetween$RIviability.SobelChen, mu = 0, alternative = "greater")
t.test(HHdataWithin$RIseedSet.SobelChen, mu = 0, alternative = "greater")
t.test(HHdataBetween$RIseedSet.SobelChen, mu = 0, alternative = "greater")
t.test(HHdataWithin$RIviability.SobelChen, mu = 0, alternative = "greater")
t.test(HHdataBetween$RIviability.SobelChen, mu = 0, alternative = "greater")
