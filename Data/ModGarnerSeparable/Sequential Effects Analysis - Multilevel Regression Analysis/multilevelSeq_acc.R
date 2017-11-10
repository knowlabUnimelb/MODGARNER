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
adapt_delta_level = .99995

# Read similarity rating data
datafile <- "mldata_sep_rt_bri_acc.dat" # ML data are repetitions vs all items except repetitions
bridataall <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rep", "acc"))

# Summarize data as counts
bridata <- aggregate(acc~sub+tsk+itm+rep, data=bridataall, FUN=sum, na.rm=T)
counts <-  aggregate(acc~sub+tsk+itm+rep, data=bridataall, FUN=length)
bridata$N <- counts$acc
bridata$err <- bridata$N - bridata$acc

datafile <- "mldata_sep_rt_idv_bri_acc.dat" # ML data are repetitions vs all items except repetitions
bridataallFI <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "acc"))

# Summarize data as counts
bridataFI <- aggregate(acc~sub+tsk+itm, data=bridataallFI, FUN=sum, na.rm=T)
counts <-  aggregate(acc~sub+tsk+itm, data=bridataallFI, FUN=length)
bridataFI$N <- counts$acc
bridataFI$err <- bridataFI$N - bridataFI$acc

datafile <- "mldata_sep_rt_sat_acc.dat" # ML data are repetitions vs all items except repetitions
satdataall <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "rep", "acc"))

# Summarize data as counts
satdata <- aggregate(acc~sub+tsk+itm+rep, data=satdataall, FUN=sum, na.rm=T)
counts <-  aggregate(acc~sub+tsk+itm+rep, data=satdataall, FUN=length)
satdata$N <- counts$acc
satdata$err <- satdata$N - satdata$acc

datafile <- "mldata_sep_rt_idv_sat_acc.dat" # ML data are repetitions vs all items except repetitions
satdataallFI <- read.table(datafile, sep = "\t", header = FALSE, row.names = NULL, col.names = c("con", "sub", "tsk", "itm", "acc"))

# Summarize data as counts
satdataFI <- aggregate(acc~sub+tsk+itm, data=satdataallFI, FUN=sum, na.rm=T)
counts <-  aggregate(acc~sub+tsk+itm, data=satdataallFI, FUN=length)
satdataFI$N <- counts$acc
satdataFI$err <- satdataFI$N - satdataFI$acc

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
bcn <- data[ data$tsk == 1 & data$itm == 1, ] # datacode: bcn = line, control, near item (b vs s - saturation, c vs l vs f - correlated vs filtering, n vs f vs a - far vs adj)
briControlRepetitionEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
           data = bcn, 
           family = binomial("logit"),
           prior = c(set_prior("normal(0, 1)", class = "b")), 
           warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
           control = list(adapt_delta = adapt_delta_level))
summary(briControlRepetitionEffect)
hypothesis(briControlRepetitionEffect, "rep = 0")
hypothesis(briControlRepetitionEffect, "rep < 0")
#########################################################################################################################################################
bcf <- data[data$tsk == 1 & data$itm == 2, ]
# briControlPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briControlPushingEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
               data = bcf, 
               family = binomial("logit"),
               prior = c(set_prior("normal(0, 1)", class = "b")), 
               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
               control = list(adapt_delta = adapt_delta_level))
summary(briControlPushingEffect)
hypothesis(briControlPushingEffect, "rep = 0")
hypothesis(briControlPushingEffect, "rep > 0")
#########################################################################################################################################################

fdt_adapt_delta_level = .999995

bca <- data[data$tsk == 1 & data$itm == 3, ]
# briControlPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briControlPullingEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                               data = bca, 
                               family = binomial("logit"),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = fdt_adapt_delta_level))
summary(briControlPullingEffect)
hypothesis(briControlPullingEffect, "rep = 0")
hypothesis(briControlPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# brightness, Correlated
# briCorrelatedRepetitionEffect <- lmer(rt ~ rep + (1|sub), bln)
bln <- data[data$tsk == 2 & data$itm == 1, ]
briCorrelatedRepetitionEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                                  data = bln, 
                                  family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(briCorrelatedRepetitionEffect)
hypothesis(briCorrelatedRepetitionEffect, "rep = 0")
hypothesis(briCorrelatedRepetitionEffect, "rep < 0")
#########################################################################################################################################################

blf <- data[data$tsk == 2 & data$itm == 2, ]
# briCorrelatedPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briCorrelatedPushingEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                               data = blf, 
                               family = binomial("logit"),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = adapt_delta_level))
summary(briCorrelatedPushingEffect)
hypothesis(briCorrelatedPushingEffect, "rep = 0")
hypothesis(briCorrelatedPushingEffect, "rep > 0")
#########################################################################################################################################################
bla <- data[data$tsk == 2 & data$itm == 3, ]
# briCorrelatedPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briCorrelatedPullingEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                               data = bla, 
                               family = binomial("logit"),
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
bfn <- data[data$tsk == 3 & data$itm == 1, ]
briFilteringRepetitionEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                                     data = bfn, 
                                    family = binomial("logit"),
                                     prior = c(set_prior("normal(0, 1)", class = "b")), 
                                     warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                     control = list(adapt_delta = adapt_delta_level))
summary(briFilteringRepetitionEffect)
hypothesis(briFilteringRepetitionEffect, "rep = 0")
hypothesis(briFilteringRepetitionEffect, "rep < 0")
#########################################################################################################################################################

bff <- data[data$tsk == 3 & data$itm == 2, ]
# briFilteringPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
briFilteringPushingEffect <- brm(formula = err | trials(N) ~ rep + (1|sub), 
                                  data = bff, 
                                 family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(briFilteringPushingEffect)
hypothesis(briFilteringPushingEffect, "rep = 0")
hypothesis(briFilteringPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .999995
bfa <- data[data$tsk == 3 & data$itm == 3, ]
# briFilteringPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
briFilteringPullingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                  data = bfa, 
                                 family = binomial("logit"),
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

# Rep vs Rep
idvN <- dataidv[dataidv$itm == 1, ]
briIrrelevantDimensionRep <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                 data = idvN, 
                                 family = binomial("logit"),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
summary(briIrrelevantDimensionRep)
hypothesis(briIrrelevantDimensionRep, "tsk = 0")
#########################################################################################################################################################
# Far vs Far
idvF <- dataidv[dataidv$itm == 2, ]
briIrrelevantDimensionPush <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                 data = idvF, 
                                 family = binomial("logit"),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
summary(briIrrelevantDimensionPush)
hypothesis(briIrrelevantDimensionPush, "tsk = 0")
#########################################################################################################################################################
# Near vs Near
idvA <- dataidv[dataidv$itm == 3, ]
briIrrelevantDimensionPull <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                  data = idvA, 
                                  family = binomial("logit"),
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
fdt_adapt_delta_level = .999995
# Saturation, Control
# SatControlRepetitionEffect <- lmer(rt ~ rep + (1|sub), bcn)
scn <- data[data$tsk == 1 & data$itm == 1, ]
SatControlRepetitionEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                  data = scn, 
                                  family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = fdt_adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatControlRepetitionEffect)
hypothesis(SatControlRepetitionEffect, "rep = 0")
hypothesis(SatControlRepetitionEffect, "rep < 0")
#########################################################################################################################################################

scf <- data[data$tsk == 1 & data$itm == 2, ]
# SatControlPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatControlPushingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                               data = scf, 
                               family = binomial("logit"),
                               prior = c(set_prior("normal(0, 1)", class = "b")), 
                               warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                               control = list(adapt_delta = adapt_delta_level))
summary(SatControlPushingEffect)
hypothesis(SatControlPushingEffect, "rep = 0")
hypothesis(SatControlPushingEffect, "rep > 0")
#########################################################################################################################################################
fdt_adapt_delta_level = .999995
sca <- data[data$tsk == 1 & data$itm == 3, ]
# SatControlPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatControlPullingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                               data = sca, 
                               family = binomial("logit"),
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
sln <- data[data$tsk == 2 & data$itm == 1, ]
SatCorrelatedRepetitionEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                     data = sln, 
                                     family = binomial("logit"),
                                     prior = c(set_prior("normal(0, 1)", class = "b")), 
                                     warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                     control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatCorrelatedRepetitionEffect)
hypothesis(SatCorrelatedRepetitionEffect, "rep = 0")
hypothesis(SatCorrelatedRepetitionEffect, "rep < 0")
#########################################################################################################################################################

slf <- data[data$tsk == 2 & data$itm == 2, ]
# SatCorrelatedPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatCorrelatedPushingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                  data = slf, 
                                  family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatCorrelatedPushingEffect)
hypothesis(SatCorrelatedPushingEffect, "rep = 0")
hypothesis(SatCorrelatedPushingEffect, "rep > 0")
#########################################################################################################################################################
sla <- data[data$tsk == 2 & data$itm == 3, ]
# SatCorrelatedPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatCorrelatedPullingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                  data = sla, 
                                  family = binomial("logit"),
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
sfn <- data[data$tsk == 3 & data$itm == 1, ]
SatFilteringRepetitionEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                    data = sfn, 
                                    family = binomial("logit"),
                                    prior = c(set_prior("normal(0, 1)", class = "b")), 
                                    warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                    control = list(adapt_delta = adapt_delta_level))
# Options to think about: use student() instead of gaussian() for robust regression
summary(SatFilteringRepetitionEffect)
hypothesis(SatFilteringRepetitionEffect, "rep = 0")
hypothesis(SatFilteringRepetitionEffect, "rep < 0")
#########################################################################################################################################################
sff <- data[data$tsk == 3 & data$itm == 2, ]
# SatFilteringPushingEffect <- lmer(rt ~ rep + (1|sub), bcf)
SatFilteringPushingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                 data = sff, 
                                 family = binomial("logit"),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
summary(SatFilteringPushingEffect)
hypothesis(SatFilteringPushingEffect, "rep = 0")
hypothesis(SatFilteringPushingEffect, "rep > 0")
#########################################################################################################################################################
sfa <- data[data$tsk == 3 & data$itm == 3, ]
# SatFilteringPullingEffect <- lmer(rt ~ rep + (1|sub), bca)
SatFilteringPullingEffect <- brm(formula = err|trials(N) ~ rep + (1|sub), 
                                 data = sfa, 
                                 family = binomial("logit"),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
summary(SatFilteringPullingEffect)
hypothesis(SatFilteringPullingEffect, "rep = 0")
hypothesis(SatFilteringPullingEffect, "rep > 0")
#########################################################################################################################################################
#########################################################################################################################################################
# Filtration vs Filtration IDV
dataidv = satdataFI

# Rep vs Rep
sidvN <- dataidv[dataidv$itm == 1, ]
SatIrrelevantDimensionRep <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                 data = sidvN, 
                                 family = binomial("logit"),
                                 prior = c(set_prior("normal(0, 1)", class = "b")), 
                                 warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                 control = list(adapt_delta = adapt_delta_level))
print(summary(SatIrrelevantDimensionRep), digits=4)
hypothesis(SatIrrelevantDimensionRep, "tsk = 0")
#########################################################################################################################################################
# Far vs Far
sidvF <- dataidv[dataidv$itm == 2, ]
SatIrrelevantDimensionPush <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                  data = sidvF, 
                                  family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatIrrelevantDimensionPush)
hypothesis(SatIrrelevantDimensionPush, "tsk = 0")
#########################################################################################################################################################
# Near vs Near
sidvA <- dataidv[dataidv$itm == 3, ]
SatIrrelevantDimensionPull <- brm(formula = err|trials(N) ~ tsk + (1|sub), 
                                  data = sidvA, 
                                  family = binomial("logit"),
                                  prior = c(set_prior("normal(0, 1)", class = "b")), 
                                  warmup = 1000, iter = 2000, chains = 4, sample_prior = TRUE, 
                                  control = list(adapt_delta = adapt_delta_level))
summary(SatIrrelevantDimensionPull)
hypothesis(SatIrrelevantDimensionPull, "tsk = 0")
#########################################################################################################################################################
#########################################################################################################################################################
allObjects = ls()
save(allObjects, file = "sep_ACC.Rdata")
