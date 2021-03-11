# calculate structural metrics and plots for a number of null models
# NOTE: this script should be repeated for each year you want to analyze
# by setting the variable 'base.year'
# (yes, a loop would do all years at once, but since calculations
# are intensive, I thought it was ok like that)

# NOTE 2: These computations are computationally expensive, so a parallel version
# using foreach is presented. Be aware of it if running the script...

# load packages -----------------------------------------------------------

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 10
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# initial info ------------------------------------------------------------
# load auxiliary functions

source("./R/reshuffle_spatial_network.R")

source("./R/null_structural_metrics.R")
source("./R/species_roles.R")

# initial data ------------------------------------------------------------
# take as a baseline the homogeneous parameters
# as with these we have only the variability associated
# to the interaction links, while space is constant

richness <- read.csv2("./data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)
abundances <- read.csv2(file = "./data/01_05_abundances.csv",header = T,stringsAsFactors = F)
lambdas <- read.csv2(file = "./data/01_05_lambda_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
alphas <- read.csv2(file = "./data/01_05_alpha_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2(file = "./data/01_05_plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
valid.sp <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

# which year to take as a baseline
base.year <- 2015
lambdas <- subset(lambdas, year == base.year)
alphas <- subset(alphas, year == base.year)
abundances <- subset(abundances, year == base.year)
richness <- subset(richness, year == base.year)

# presence per plot
presence <- abundances %>% 
  filter(species %in% valid.sp) %>% 
  group_by(plot,species) %>% 
  summarise(ind = sum(individuals))
presence$presence <- presence$ind > 0
presence_df <- presence[,c("plot","species","presence")]
names(presence_df)[1] <- "site"

# for the function input
alphas <- alphas[,2:5]
names(alphas)[1] <- "site"

# resulting structures
null_metrics <- NULL
null_sproles <- NULL


# parameters --------------------------------------------------------------

num.replicates <- 60 # for an exact number of CPUs in cluster
models <- c("KD","RA") # reshuffle all, keep diagonal

# ID to loop over ---------------------------------------------------------

id <- expand.grid(models,base.year,1:num.replicates)
id.char <- paste(id[,1],"_",id[,2],"_",id[,3],sep="")

# function for combining foreach output -----------------------------------
# from https://stackoverflow.com/questions/28348749/outptut-two-objects-using-foreach

# in my code, I need to return two dataframes, 
# that will be combined row-wise

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

# obtain null realizations ------------------------------------------------

results <- foreach(i.id = 1:length(id.char), 
                   .combine=comb.fun, 
                   .packages = 'tidyverse') %dopar% {
                     
                     # recover year, model, and replicate from the ID
                     
                     i.year <- substr(id.char[i.id],4,7)
                     i.model <- substr(id.char[i.id],1,2)
                     i.rep <- substr(id.char[i.id],9,nchar(id.char[i.id]))
                     
                     keep.diag <- ifelse(i.model == "RA",FALSE,TRUE)
                     # generate null models ----------------------------------------------------
                     
                     null_edge_list <- reshuffle_spatial_network(edge_list = alphas,
                                                                 scale = "network",
                                                                 keep.diagonal = keep.diag)  
                     # recover original names
                     names(null_edge_list)[1] <- "plot"
                     null_edge_list$year <- base.year
                     null_edge_list <- null_edge_list[,c("year","plot","focal","neighbour","magnitude")]
                     null_edge_list_valid <- null_edge_list[complete.cases(null_edge_list),]
                     
                     # calculate structural metrics --------------------------------------------
                     # use the null version of the function
                     
                     # careful to include only valid links in alphas
                     metrics <- null_structural_metrics(lambdas = lambdas, 
                                                        alphas = null_edge_list_valid, 
                                                        sp.rates = sp.rates)
                     
                     # calculate species roles -------------------------------------------------
                     
                     sproles <- species_roles(structural_metrics = metrics,
                                              lambdas = lambdas,
                                              sp.rates = sp.rates,
                                              abundances = abundances)
                     
                     # after calculating species roles, 
                     # keep only feasible communities
                     metrics_feas <- subset(metrics, feasibility == 1)
                     
                     # add id columns
                     sproles$replicate <- i.rep
                     sproles$year <- i.year
                     sproles$null_model <- i.model
                     
                     metrics_feas$replicate <- i.rep
                     metrics_feas$year <- i.year
                     metrics_feas$null_model <- i.model
                     
                     list(metrics_feas,sproles)
                   }


null_metrics <- results[[1]]
null_sproles <- results[[2]]

write.csv2(null_metrics,file = paste("./results/null_replicates/null_structural_metrics_",base.year,".csv",sep=""),row.names = FALSE)
write.csv2(null_sproles,file = paste("./results/null_replicates/null_species_roles_",base.year,".csv",sep=""),row.names = FALSE)

stopCluster(cl)
