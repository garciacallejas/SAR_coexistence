# calculate structural metrics and plots for the observed values

# initial info ------------------------------------------------------------

library(tidyverse)
source("./R/structural_metrics.R")
source("./R/species_roles.R")

# initial data ------------------------------------------------------------

richness <- read.csv2("./data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)
abundances <- read.csv2(file = "./data/01_05_abundances.csv",header = T,stringsAsFactors = F)

hetlambdas <- read.csv2(file = "./data/01_lambda_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
hetalphas <- read.csv2(file = "./data/01_alpha_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)

homlambdas <- read.csv2(file = "./data/01_05_lambda_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
homalphas <- read.csv2(file = "./data/01_05_alpha_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)

sp.rates <- read.csv2(file = "./data/01_05_plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)

# obtain structural metrics -----------------------------------------------

# homogeneous space

hom_metrics <- structural_metrics(lambdas = homlambdas,
                                  alphas = homalphas,
                                  sp.rates = sp.rates)
hom_metrics$fit <- "homogeneous_space"

# heterogeneous space

het_metrics <- structural_metrics(lambdas = hetlambdas,
                                  alphas = hetalphas,
                                  sp.rates = sp.rates)
het_metrics$fit <- "heterogeneous_space"

# single dataframe

observed_metrics <- rbind(hom_metrics,het_metrics)

# species roles in the heterogeneous curve --------------------------------

het_sproles <- species_roles(structural_metrics = het_metrics,
                             lambdas = hetlambdas,
                             sp.rates = sp.rates,
                             abundances = abundances)

# write to disk -----------------------------------------------------------

write.csv2(observed_metrics,"./results/02_observed_structural_metrics.csv",row.names = FALSE)
write.csv2(het_sproles,"./results/04_observed_species_roles.csv",row.names = FALSE)
