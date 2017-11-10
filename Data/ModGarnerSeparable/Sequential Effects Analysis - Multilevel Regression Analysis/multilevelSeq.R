# Analysis for Garner Separable paper
rm(list = ls()) # Clear the global workspace
cat("\014")     # Clear the screen

# Change to current directory
directory = file.path("C:", "Users", "littled", "Dropbox", "Work", "2017 Garner", "MODIFIED GARNER SEPARABLE", "Sequential Effects Analysis - Multilevel Regression Analysis")
setwd(directory)

#install.packages("lme4")
#install.packages("lmerTest")
# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE) # To use rstan follow instructions here: https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows
# install.packages("brms") 

library(lmerTest)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(brms)
library(shinystan)

# 
adapt_delta_level = .9995

# Read similarity rating data
datafile <- "mldata_sep_rt_bri.dat" # ML data are repetitions vs all items except repetitions
bridata <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rep", "rt"))
bridata$logrt = log(bridata$rt)

datafile <- "mldata_sep_rt_idv_bri.dat" # ML data are repetitions vs all items except repetitions
bridataFI <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rt"))
bridataFI$logrt = log(bridataFI$rt)

datafile <- "mldata_sep_rt_sat.dat" # ML data are repetitions vs all items except repetitions
satdata <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rep", "rt"))
satdata$logrt = log(satdata$rt)

datafile <- "mldata_sep_rt_idv_sat.dat" # ML data are repetitions vs all items except repetitions
satdataFI <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rt"))
satdataFI$logrt = log(satdataFI$rt)

# Multilevel regression model
# Fixed Factors: con, tsk, itm, rep, Random Factors: sub
# lmer - use Sattherwaite correction - to get p-values: m1 <- lmer(rt~con + tsk + itm + rep + (1|sub), data)

#########################################################################################################################################################
# Useful plot check
# bcn <- data[data$con == condition & data$tsk == 1 & data$itm == 1, ] # near items from the line control task
# bcn['label'] = "rep"
# bcn$label[bcn$rep == 0] = "nonrep"
# ggplot(bcn, aes(rt, fill = label)) + geom_density(alpha = .2)
#########################################################################################################################################################

# Bayesian Regression
# Set sample prior on to compute Bayes Factors via Savage Dickey ratio using hypothesis
condition = 1  # line condition
data = bridata

# line, Control
# briControlRepetitionEffect <- lmer(rt ~ rep + (1|sub), bcn)
bcn <- data[data$con == condition & data$tsk == 1 & data$itm == 1, ] # datacode: bcn = brightnes, control, near item (b vs s - saturation, c vs l vs f - correlated vs filtering, n vs f vs a - far vs adj)
briControlRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
           data = bcn, 
           family = student(),
           prior = c(set_prior("normal(0, 1)", class = "b")), 
           warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
           control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(briControlRepetitionEffect)
hypothesis(briControlRepetitionEffect, "rep = 0")
hypothesis(briControlRepetitionEffect, "rep < 0")

# The coefficient for rep is estimated to be -.11 which implies that there is about 
# an 11% decrease in RT when there is a repetition

#########################################################################################################################################################
bcf <- data[data$con == condition & data$tsk == 1 & data$itm == 2, ]
# briControlPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briControlPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
               data = bcf, 
               family = student(),
               prior = c(set_prior("normal(0, 1)", class = "b")), 
               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
               control = list(adapt_delta = adapt_delta_level))
summary(briControlPushingEffect)
hypothesis(briControlPushingEffect, "rep = 0")
hypothesis(briControlPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .9999995
bca <- data[data$con == condition & data$tsk == 1 & data$itm == 3, ]
# briControlPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briControlPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                               data = bca, 
                               family = student(),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = fdt_adapt_delta_level))
summary(briControlPullingEffect)
hypothesis(briControlPullingEffect, "rep = 0")
hypothesis(briControlPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# line, Correlated
# briCorrelatedRepetitionEffect <- lmer(rt ~ rep + (1|sub), bln)
bln <- data[data$con == condition & data$tsk == 2 & data$itm == 1, ]
briCorrelatedRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = bln, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(briCorrelatedRepetitionEffect)
hypothesis(briCorrelatedRepetitionEffect, "rep = 0")
hypothesis(briCorrelatedRepetitionEffect, "rep < 0")
#########################################################################################################################################################

fdt_adapt_delta_level = .9999995
blf <- data[data$con == condition & data$tsk == 2 & data$itm == 2, ]
# briCorrelatedPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briCorrelatedPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                               data = blf, 
                               family = student(),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = fdt_adapt_delta_level))
summary(briCorrelatedPushingEffect)
hypothesis(briCorrelatedPushingEffect, "rep = 0")
hypothesis(briCorrelatedPushingEffect, "rep > 0")
#########################################################################################################################################################
bla <- data[data$con == condition & data$tsk == 2 & data$itm == 3, ]
# briCorrelatedPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briCorrelatedPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                               data = bla, 
                               family = student(),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = adapt_delta_level))
summary(briCorrelatedPullingEffect)
hypothesis(briCorrelatedPullingEffect, "rep = 0")
hypothesis(briCorrelatedPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# line, Filtering
# briFilteringRepetitionEffect <- lmer(rt ~ rep + (1|sub), bln)
bfn <- data[data$con == condition & data$tsk == 3 & data$itm == 1, ]
briFilteringRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                     data = bfn, 
                                     family = student(),
                                     prior = c(set_prior("normal(0, 1)", class = "b")), 
                                     warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                     control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(briFilteringRepetitionEffect)
hypothesis(briFilteringRepetitionEffect, "rep = 0")
hypothesis(briFilteringRepetitionEffect, "rep < 0")
#########################################################################################################################################################
bff <- data[data$con == condition & data$tsk == 3 & data$itm == 2, ]
# briFilteringPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briFilteringPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = bff, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(briFilteringPushingEffect)
hypothesis(briFilteringPushingEffect, "rep = 0")
hypothesis(briFilteringPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .99995
bfa <- data[data$con == condition & data$tsk == 3 & data$itm == 3, ]
# briFilteringPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briFilteringPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = bfa, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = fdt_adapt_delta_level))
summary(briFilteringPullingEffect)
hypothesis(briFilteringPullingEffect, "rep = 0")
hypothesis(briFilteringPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# Filtration vs Filtration IDV
dataidv = bridataFI
fdt_adapt_delta_level = .99995

# Rep vs Rep
idvN <- dataidv[dataidv$con == condition & dataidv$itm == 1, ]
briIrrelevantDimensionRep <- brm(formula = logrt ~ tsk + (1|sub), 
                                 data = idvN, 
                                 family = student(),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = fdt_adapt_delta_level))
summary(briIrrelevantDimensionRep)
hypothesis(briIrrelevantDimensionRep, "tsk = 0")
#########################################################################################################################################################
# Far vs Far
fdt_adapt_delta_level = .99995
idvF <- dataidv[dataidv$con == condition & dataidv$itm == 2, ]
briIrrelevantDimensionPush <- brm(formula = logrt ~ tsk + (1|sub), 
                                 data = idvF, 
                                 family = student(),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = fdt_adapt_delta_level))
summary(briIrrelevantDimensionPush)
hypothesis(briIrrelevantDimensionPush, "tsk = 0")
#########################################################################################################################################################
# Near vs Near
idvA <- dataidv[dataidv$con == condition & dataidv$itm == 3, ]
briIrrelevantDimensionPull <- brm(formula = logrt ~ tsk + (1|sub), 
                                  data = idvA, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(briIrrelevantDimensionPull)
hypothesis(briIrrelevantDimensionPull, "tsk = 0")
#########################################################################################################################################################
#########################################################################################################################################################

# Bayesian Regression
# Set sample prior on to compute Bayes Factors via Savage Dickey ratio using hypothesis
condition = 2  # Saturation condition
data = satdata

# Saturation, Control
# SatControlRepetitionEffect <- lmer(rt ~ rep + (1|sub), bcn)
scn <- data[data$con == condition & data$tsk == 1 & data$itm == 1, ]
SatControlRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = scn, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatControlRepetitionEffect)
hypothesis(SatControlRepetitionEffect, "rep = 0")
hypothesis(SatControlRepetitionEffect, "rep < 0")
#########################################################################################################################################################
scf <- data[data$con == condition & data$tsk == 1 & data$itm == 2, ]
# SatControlPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatControlPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                               data = scf, 
                               family = student(),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = adapt_delta_level))
summary(SatControlPushingEffect)
hypothesis(SatControlPushingEffect, "rep = 0")
hypothesis(SatControlPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .9999995
sca <- data[data$con == condition & data$tsk == 1 & data$itm == 3, ]
# SatControlPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatControlPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                               data = sca, 
                               family = student(),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = fdt_adapt_delta_level))
summary(SatControlPullingEffect)
hypothesis(SatControlPullingEffect, "rep = 0")
hypothesis(SatControlPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# Saturation, Correlated
# SatCorrelatedRepetitionEffect <- lmer(rt ~ rep + (1|sub), bln)
sln <- data[data$con == condition & data$tsk == 2 & data$itm == 1, ]
SatCorrelatedRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                     data = sln, 
                                     family = student(),
                                     prior = c(set_prior("normal(0, 1)", class = "b")), 
                                     warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                     control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatCorrelatedRepetitionEffect)
hypothesis(SatCorrelatedRepetitionEffect, "rep = 0")
hypothesis(SatCorrelatedRepetitionEffect, "rep < 0")
#########################################################################################################################################################
slf <- data[data$con == condition & data$tsk == 2 & data$itm == 2, ]
# SatCorrelatedPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatCorrelatedPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = slf, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatCorrelatedPushingEffect)
hypothesis(SatCorrelatedPushingEffect, "rep = 0")
hypothesis(SatCorrelatedPushingEffect, "rep > 0")
#########################################################################################################################################################
sla <- data[data$con == condition & data$tsk == 2 & data$itm == 3, ]
# SatCorrelatedPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatCorrelatedPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                  data = sla, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatCorrelatedPullingEffect)
hypothesis(SatCorrelatedPullingEffect, "rep = 0")
hypothesis(SatCorrelatedPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# Saturation, Filtering
# SatFilteringRepetitionEffect <- lmer(rt ~ rep + (1|sub), bln)
sfn <- data[data$con == condition & data$tsk == 3 & data$itm == 1, ]
SatFilteringRepetitionEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                    data = sfn, 
                                    family = student(),
                                    prior = c(set_prior("normal(0, 1)", class = "b")), 
                                    warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                    control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatFilteringRepetitionEffect)
hypothesis(SatFilteringRepetitionEffect, "rep = 0")
hypothesis(SatFilteringRepetitionEffect, "rep < 0")
#########################################################################################################################################################
sff <- data[data$con == condition & data$tsk == 3 & data$itm == 2, ]
# SatFilteringPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatFilteringPushingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                 data = sff, 
                                 family = student(),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
summary(SatFilteringPushingEffect)
hypothesis(SatFilteringPushingEffect, "rep = 0")
hypothesis(SatFilteringPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .99995
sfa <- data[data$con == condition & data$tsk == 3 & data$itm == 3, ]
# SatFilteringPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatFilteringPullingEffect <- brm(formula = logrt ~ rep + (1|sub), 
                                 data = sfa, 
                                 family = student(),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = fdt_adapt_delta_level))
summary(SatFilteringPullingEffect)
hypothesis(SatFilteringPullingEffect, "rep = 0")
hypothesis(SatFilteringPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# Filtration vs Filtration IDV
dataidv = satdataFI

# Rep vs Rep
sidvN <- dataidv[dataidv$con == condition & dataidv$itm == 1, ]
SatIrrelevantDimensionRep <- brm(formula = logrt ~ tsk + (1|sub), 
                                 data = sidvN, 
                                 family = student(),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
print(summary(SatIrrelevantDimensionRep), digits=4)
hypothesis(SatIrrelevantDimensionRep, "tsk = 0")
#########################################################################################################################################################
# Far vs Far
fdt_adapt_delta_level = .9999995
sidvF <- dataidv[dataidv$con == condition & dataidv$itm == 2, ]
SatIrrelevantDimensionPush <- brm(formula = logrt ~ tsk + (1|sub), 
                                  data = sidvF, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = fdt_adapt_delta_level))
summary(SatIrrelevantDimensionPush)
hypothesis(SatIrrelevantDimensionPush, "tsk = 0")
#########################################################################################################################################################
# Near vs Near
sidvA <- dataidv[dataidv$con == condition & dataidv$itm == 3, ]
SatIrrelevantDimensionPull <- brm(formula = logrt ~ tsk + (1|sub), 
                                  data = sidvA, 
                                  family = student(),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatIrrelevantDimensionPull)
hypothesis(SatIrrelevantDimensionPull, "tsk = 0")
#########################################################################################################################################################
#########################################################################################################################################################
allObjects = ls()
save(allObjects, file = "sep_RT.Rdata")
