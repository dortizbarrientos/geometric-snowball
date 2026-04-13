# Load libraries
library(lme4)
library(ggplot2)
library(extrafont)

# Read data file
data <- read.csv("~/Transplant.csv")

#### Seedling Establishment ####

# Extract and prepare Dune environment data
# Filter for records where both ENV and Type are "DUNE"
DUNE <- droplevels(data[data$ENV %in% "DUNE" & data$Type %in% "DUNE",])

# Extract and prepare Headland environment data  
# Filter for records where both ENV and Type are "HEAD"
HEAD <- droplevels(data[data$ENV %in% "HEAD" & data$Type %in% "HEAD",])

# Check current factor levels for TypeS in DUNE data
levels(DUNE$TypeS)

# Rename population codes for consistency
# Change "D1" to "AD1" in Dune data
DUNE$TypeS <- gsub("D1", "AD1", DUNE$TypeS)
# Change "H1" to "AH1" in Headland data  
HEAD$TypeS <- gsub("H1", "AH1", HEAD$TypeS)
# Verify the factor level changes
levels(as.factor(DUNE$TypeS))

# Convert variables to factors
DUNE$TypeS <- as.factor(DUNE$TypeS)
DUNE$block <- as.factor(DUNE$block)

# Convert variables to factors
HEAD$TypeS <- as.factor(HEAD$TypeS)
HEAD$block <- as.factor(HEAD$block)

# Fit generalized linear mixed-effects model for seedling establishment for Dune data
# L10: Binary response variable (seedling establishment success)
# TypeS: Fixed effect (population)
# block: Random effect to account for experimental design
m1 <- glmer(L10 ~ TypeS + (1|block), family="binomial", data=DUNE)
summary(m1)

# Plot the model effects
p1 <- plot_model(m1, type = "eff", terms = "TypeS", colors=c("#FFFFFF", "black", "black"), title = "Dune environment")

# Customize the plot for publication
p1 <- p1 + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        axis.ticks.length = unit(.17, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.4),
        panel.background = element_blank()) + 
  theme(text = element_text(family = "Helvetica Neue Light", size = 17)) +  
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.position = "none") +
  labs(y = "Probability of seedling establishment (%)", x = "Population", title = "Dune environment (LH)") + theme(plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(ylim = c(0.10, 0.65)) +
  
  # Add custom points on top of the existing plot
  geom_point(data = p1$data, 
             aes(x = x, y = predicted), 
             size = 3, 
             color = c("#0072B1", "black", "black")[as.numeric(p1$data$x)]) +
  
  # Add custom error bars on top of the existing plot
  geom_errorbar(data = p1$data, 
                aes(x = x, ymin = conf.low, ymax = conf.high), 
                width = 0.1, 
                color = c("#0072B1", "black", "black")[as.numeric(p1$data$x)])

p1

ppi <- 300
png("Dune-SeedlingEstablishment.png", width=5.5*ppi, height=4.5*ppi, res=ppi)
plot(p1)
dev.off()

# Fit generalized linear mixed-effects model for Headland data
m2 <- glmer(L10 ~ TypeS + (1|block), family="binomial", data=HEAD)
summary(m2)

# Plot the model effects
p2 <- plot_model(m2, type = "eff", terms = "TypeS", colors=c("#FFFFFF", "black", "black"), title = "Headland environment")

# Customize the plot for publication
p2 <- p2 + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        axis.ticks.length = unit(.17, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.4),
        panel.background = element_blank()) + 
  theme(text = element_text(family = "Helvetica Neue Light", size = 17)) +  
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.position = "none") +
  labs(y = "Probability of seedling establishment (%)", x = "Population", title = "Headland environment (LH)") + theme(plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(ylim = c(0.10, 0.65)) +
  
  # Add custom points on top of the existing plot
  geom_point(data = p2$data, 
             aes(x = x, y = predicted), 
             size = 3, 
             color = c("#D45E00", "black", "black")[as.numeric(p2$data$x)]) +
  
  # Add custom error bars on top of the existing plot
  geom_errorbar(data = p2$data, 
                aes(x = x, ymin = conf.low, ymax = conf.high), 
                width = 0.1, 
                color = c("#D45E00", "black", "black")[as.numeric(p2$data$x)])

p2

ppi <- 300
png("Headland-SeedlingEstablishment.png", width=5.5*ppi, height=4.5*ppi, res=ppi)
plot(p2)
dev.off()

#### Survival curves ####

# Load additional libraries for survival analysis
library(survival)
library(coxme)
library(multcomp)

# Fit Cox proportional hazards mixed-effects model for Dune data
# Surv(totaldays, CD320): Survival object with time and event indicator
# TypeS: Population type as fixed effect
# block: Random effect for experimental blocks
model1 <- coxme(Surv(totaldays,CD320) ~ TypeS + (1|block), data = DUNE)

# Extract population names for reference
popnames <- as.character(DUNE$TypeS)
popnames <- levels(DUNE$TypeS)

# Fit survival curves for each population for the Dune data
fit1 <- survfit(Surv(totaldays,CD320) ~ TypeS, data = DUNE)

# Perform pairwise comparisons between population types using Tukey's method
# Save results to text file
sink(file = "ghlt_DUNE.txt")
summary(glht(model1, linfct = mcp(TypeS = "Tukey")))
sink()

# Fit Cox proportional hazards mixed-effects model for Headland data
model2 <- coxme(Surv(totaldays,CD320) ~  TypeS + (1|block), data = HEAD)

# Extract population names for reference
popnames <- as.character(HEAD$TypeS)
popnames <- levels(HEAD$TypeS)

# Fit survival curves for each population for the Headland data
fit2 <- survfit(Surv(totaldays, CD320) ~ TypeS, data = HEAD)

# Perform pairwise comparisons between population types using Tukey's method
# Save results to text file
sink(file = "ghlt_HEAD.txt")
summary(glht(model2, linfct = mcp(TypeS = "Tukey")))
sink()


## Dune graph ##

# Set up the PNG output
png("Dune-Survivors.png", height = 6, width = 5.5, units = "in", res = 300)

# Set the font to Helvetica Neue Light using par()
par(family = "Helvetica Neue Light")

# Plot with the specified font
plot(fit1, col = c("#0072B1", "black", "black"), 
     xlab = "Time (days)", ylab = "Estimated S(t)", 
     main = "Dune Environment", lwd = 1.5, conf.int = TRUE, ylim = c(0, 0.9))

# Reset the font to the default after the plot if needed
par(family = "default")

# Close the PNG device
dev.off()

## Headland graph ##

# Set up the PNG output
png("Heaadland-Survivors.png", height = 6, width = 5.5, units = "in", res = 300)

# Set the font to Helvetica Neue Light using par()
par(family = "Helvetica Neue Light")

# Plot with the specified font
plot(fit2, col = c("#D45E00", "black", "black"), 
     xlab = "Time (days)", ylab = "Estimated S(t)", 
     main = "Headland Environment", lwd = 1.5, conf.int = TRUE, ylim = c(0, 0.9)); 

# Reset the font to the default after the plot if needed
par(family = "default")

# Close the PNG device
dev.off()

# Plots
plot(fit1, col = c("#0072B1","black","black"), xlab = "Time (days)", ylab = "Estimated S(t)", main = "Dune Environment", lwd = 2, conf.int=TRUE)
plot(fit2, col = c("chartreuse3","black","black"), xlab = "Time (days)", ylab = "Estimated S(t)", main = "Headland Environment", lwd = 2, conf.int=TRUE)


