# calculate structural metrics and plots for a number of random draws

# NOTE: this script should be repeated for each year you want to analyze
# by setting the variable 'base.year'
# and for each parameterization by setting "fit.type"
# (yes, a loop would do all years at once, but since calculations
# are intensive, I thought it was ok like that)

# NOTE 2: These computations are computationally expensive, so a parallel version
# using foreach is presented. Be aware of it if running the script...

# NOTE 4: as it takes a long time, I use a function "structural_metrics_lite"
# that returns only the minimum necessary to build the coexistence-are curves
# in script S3_03

# load packages -----------------------------------------------------------

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# initial info ------------------------------------------------------------
# load auxiliary functions

source("./R/structural_metrics_lite.R")
source("./R/species_roles.R")
source("./R/toolbox_coexistence.R")

# initial data ------------------------------------------------------------

fit.type <- "heterogeneous_time"
# fit.type <- "heterogeneous_both"

richness <- read.csv2("data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)
abundances <- read.csv2(file = "data/01_05_abundances.csv",header = T,stringsAsFactors = F)
lambdas <- read.csv2(file = paste("data/posterior_draws/S3_lambda_posterior_draws_",fit.type,".csv",sep=""),header = TRUE,stringsAsFactors = FALSE)
alphas <- read.csv2(file = paste("data/posterior_draws/S3_alpha_posterior_draws_",fit.type,".csv",sep=""),header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2(file = "data/01_05_plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
valid.sp <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

# which year to take as a baseline
base.year <- 2019
lambdas <- subset(lambdas, year == base.year)
alphas <- subset(alphas, year == base.year)
abundances <- subset(abundances, year == base.year)
richness <- subset(richness, year == base.year)

# resulting structures
posterior_draws_metrics <- NULL
posterior_draws_sproles <- NULL

# parameters --------------------------------------------------------------

replicates <-  1:96#max(lambdas$sample)

# ID to loop over ---------------------------------------------------------

id <- expand.grid(base.year,replicates)
id.char <- paste(id[,1],"_",id[,2],sep="")

# function for combining foreach output -----------------------------------
# from https://stackoverflow.com/questions/28348749/outptut-two-objects-using-foreach

# in my code, I need to return two dataframes, 
# that will be combined row-wise

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

# obtain null realizations ------------------------------------------------
cat(format(Sys.time(),usetz = TRUE)," - STARTING:",base.year,"-",fit.type,"with",workers,"cores and",length(replicates),"replicates")

results <- foreach(i.id = 1:length(id.char), 
                   .combine=comb.fun, 
                   .packages = 'tidyverse') %dopar% {
                     
                     # recover year, replicate from the ID
                     
                     i.year <- substr(id.char[i.id],1,4)
                     i.rep <- substr(id.char[i.id],6,nchar(id.char[i.id]))
                     
                     # -------------------------------------------------------------------------
                     
                     my.lambdas <- subset(lambdas, sample == i.rep)
                     my.lambdas <- my.lambdas[,c("year","plot","sp","lambda")]
                     
                     my.alphas <- subset(alphas, sample == i.rep)
                     my.alphas <- my.alphas[,c("year","plot","focal","neighbour","magnitude")]
                     
                     # calculate structural metrics --------------------------------------------
                     
                     metrics <- structural_metrics_lite(lambdas = my.lambdas, 
                                                        alphas = my.alphas, 
                                                        sp.rates = sp.rates)
                     
                     # calculate species roles -------------------------------------------------
                     
                     sproles <- species_roles(structural_metrics = metrics,
                                              lambdas = my.lambdas,
                                              sp.rates = sp.rates,
                                              abundances = abundances)
                     
                     # after calculating species roles, 
                     # keep only feasible communities
                     metrics_feas <- subset(metrics, feasibility == 1)
                     
                     # add id columns
                     sproles$replicate <- i.rep
                     
                     metrics_feas$replicate <- i.rep
                     
                     list(metrics_feas,sproles)
                     }

cat("\n",format(Sys.time(),usetz = TRUE)," - FINISHED:",base.year,"-",fit.type,"with",workers,"cores")

posterior_draws_metrics <- results[[1]]
posterior_draws_metrics$fit.type <- fit.type

posterior_draws_sproles <- results[[2]]
posterior_draws_sproles$fit.type <- fit.type

write.csv2(posterior_draws_metrics,file = paste("./results/posterior_draws/S3_posterior_draws_structural_metrics_",base.year,"_",fit.type,".csv",sep=""),row.names = FALSE)
write.csv2(posterior_draws_sproles,file = paste("./results/posterior_draws/S3_posterior_draws_species_roles_",base.year,"_",fit.type,".csv",sep=""),row.names = FALSE)

stopCluster(cl)
