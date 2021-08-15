# calculate structural metrics and plots 
# for modified germination and seed survival rates

# NOTE: this script should be repeated for each year you want to analyze
# by setting the variable 'base.year'
# and for each parameterization by setting "fit.type"
# (yes, a loop would do all years at once, but since calculations
# are intensive, I thought it was ok like that)

# NOTE 2: These computations are computationally expensive, so a parallel version
# using foreach is presented. Be aware of it if running the script...

# NOTE 3: as it takes a long time, I use a function "structural_metrics_lite"
# that returns only the minimum necessary to build the coexistence-are curves

# load packages -----------------------------------------------------------

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 6
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# initial info ------------------------------------------------------------
# load auxiliary functions

source("./R/structural_metrics_lite.R")
source("./R/species_roles.R")
source("./R/toolbox_coexistence.R")

# initial data ------------------------------------------------------------
# set parameterization

# fit.type <- "heterogeneous_time"
fit.type <- "heterogeneous_both"

# which year to take as a baseline
base.year <- 2018

# replicates
replicates <- 1:96

# -------------------------------------------------------------------------
richness <- read.csv2("data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)
abundances <- read.csv2(file = "data/01_05_abundances.csv",header = T,stringsAsFactors = F)

if(fit.type == "heterogeneous_both"){
  lambdas <- read.csv2(file = "./data/01_lambda_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
  alphas <- read.csv2(file = "./data/01_alpha_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
}else{
  lambdas <- read.csv2(file = "./data/01_05_lambda_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
  alphas <- read.csv2(file = "./data/01_05_alpha_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
}

sim.rates <- read.csv2(file = "data/rate_sensitivity_analysis/S4_01_modified_vital_rates.csv",
                       header = TRUE,stringsAsFactors = FALSE)
analyses <- unique(sim.rates$analysis)
valid.sp <- unique(sim.rates$species.code)

lambdas <- subset(lambdas, year == base.year)
alphas <- subset(alphas, year == base.year)
abundances <- subset(abundances, year == base.year)
richness <- subset(richness, year == base.year)

# ID to loop over ---------------------------------------------------------

# switch the number of leading zeros depending on the desired number of replicates
id <- expand.grid(base.year,formatC(replicates, width = 2, format = "d", flag = "0"),analyses)
id.char <- paste(id[,1],"_",id[,2],"_",id[,3],sep="")

# function for combining foreach output -----------------------------------
# from https://stackoverflow.com/questions/28348749/outptut-two-objects-using-foreach

# in my code, I need to return two dataframes, 
# that will be combined row-wise

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

# obtain null realizations ------------------------------------------------
cat(format(Sys.time(),usetz = TRUE)," - STARTING:",base.year,"-",fit.type,"with",workers,"cores")

modified_rates_metrics <- foreach(i.id = 1:length(id.char), 
                                  .combine=rbind, 
                                  .packages = 'tidyverse') %dopar% {
                                    
                                    # recover year, replicate, analysis from the ID
                                    
                                    i.year <- substr(id.char[i.id],1,4)
                                    i.rep <- as.numeric(substr(id.char[i.id],6,7))
                                    i.analysis <- substr(id.char[i.id],9,nchar(id.char[i.id]))
                                    
                                    # -------------------------------------------------------------------------
                                    
                                    my.lambdas <- lambdas[,c("year","plot","sp","lambda")]
                                    my.alphas <- alphas[,c("year","plot","focal","neighbour","magnitude")]
                                    
                                    my.rates <- subset(sim.rates, analysis == i.analysis &
                                                         replicate == i.rep)
                                    
                                    # if this analysis-replicate combination is present, go on
                                    if(nrow(my.rates)>0){
                                      
                                      # transform rates to wide format
                                      my.rates.w <- my.rates %>% select(species.code,rate,value) %>%
                                        pivot_wider(names_from = rate,values_from = value)
                                      
                                      # calculate structural metrics --------------------------------------------
                                      
                                      metrics <- structural_metrics_lite(lambdas = my.lambdas, 
                                                                         alphas = my.alphas, 
                                                                         sp.rates = my.rates.w)
                                      
                                      # calculate species roles -------------------------------------------------
                                      # this is not necessary for this analysis
                                      
                                      # sproles <- species_roles(structural_metrics = metrics,
                                      #                          lambdas = my.lambdas,
                                      #                          sp.rates = my.rates.w,
                                      #                          abundances = abundances)
                                      
                                      # keep only feasible communities
                                      metrics_feas <- subset(metrics, feasibility == 1)
                                      
                                      # add id columns
                                      metrics_feas$replicate <- i.rep
                                      metrics_feas$analysis <- i.analysis
                                      
                                      metrics_feas
                                      
                                    }else{
                                      NULL
                                    }# if-else rep-analysis combination is present
                                    
                                  }# foreach

cat("\n",format(Sys.time(),usetz = TRUE)," - FINISHED:",base.year,"-",fit.type,"with",workers,"cores")

modified_rates_metrics$fit.type <- fit.type

write.csv2(modified_rates_metrics,file = paste("./results/rate_sensitivity_analysis/S4_modified_rates_structural_metrics_",base.year,"_",fit.type,".csv",sep=""),row.names = FALSE)

stopCluster(cl)
beepr::beep(3)
