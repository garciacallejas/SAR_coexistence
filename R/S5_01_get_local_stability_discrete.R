
library(tidyverse)

# initial data ------------------------------------------------------------
# set parameterization

# fit.type <- "homogeneous"
fit.type <- "heterogeneous"

abundances <- read.csv2(file = "data/01_05_abundances.csv",header = T,stringsAsFactors = F)
sp.rates <- read.csv2(file = "data/01_05_plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)

if(fit.type == "heterogeneous"){
  lambdas <- read.csv2(file = "./data/01_lambda_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
  alphas <- read.csv2(file = "./data/01_alpha_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
}else{
  lambdas <- read.csv2(file = "./data/01_05_lambda_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
  alphas <- read.csv2(file = "./data/01_05_alpha_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------
# auxiliary functions

# transform lambdas from annual plant models to LV growth rates (r)
growth_rate_equivalence <- function(lambda,g,s,sign = "negative"){
  if(sign == "negative"){
    r <- log((1-(1-g)*s)/g) - lambda
  }else{
    r <- log((1-(1-g)*s)/g) + lambda
  }
  r
}

# Function to evaluate the steady-state of a linear LV
model_stode = function(t,y,gr,parms=NULL,A) {
  # competition is encoded by alphas > 0, be aware
  dy = y*(gr-A%*%y)
  return(list(dy,1))
}

# Function to evaluate the Jacobian matrix of a linear LV	
model_J = function(t,y,r,parms=NULL,A) {
  # competition is encoded by alphas > 0, be aware
  dy = y*(r-A%*%y)
  return(as.list(dy))
}

# -------------------------------------------------------------------------
years <- unique(lambdas$year)
n.plots <- length(unique(lambdas$plot))
all.sp <- unique(lambdas$sp)

# -------------------------------------------------------------------------
# alphas as read here have signs switched,
# so competition are alphas > 0
# this is important later for calculating the
# growth rate equivalence between lambdas and growth rates of linear LV

# lambdas are stored in "real" dimensions, i.e. number of seeds
# but the linear LV works with model raw parameters
# so here, undo the transformation by taking log again

lambdas$lambda.log <- log(lambdas$lambda)

# -------------------------------------------------------------------------

local.stab.list <- list()

for(i.year in 1:length(years)){
  
  lambda.year <- subset(lambdas,year == years[i.year])
  alpha.year <- subset(alphas,year == years[i.year])
  n.plots <- max(c(unique(lambda.year$plot),unique(alpha.year$plot)))
  
  for(i.plot in 1:n.plots){
    
    lambda.plot <- subset(lambda.year,plot == i.plot)
    alpha.plot <- subset(alpha.year,plot == i.plot)
    
    # present species
    lambda.plot <- subset(lambda.plot, !is.na(lambda))
    my.sp <- sort(unique(lambda.plot$sp))
    alpha.plot <- subset(alpha.plot, focal %in% my.sp & neighbour %in% my.sp)              
    
    # there are some NAs left, on species that are suppossed to be present
    # so, set these interactions to zero
    alpha.plot$magnitude[which(is.na(alpha.plot$magnitude))] <- 0
    
    # obtain the interaction matrix
    alpha.matrix <- tidyr::spread(alpha.plot,key = neighbour, value = magnitude)
    alpha.matrix <- as.matrix(alpha.matrix[,my.sp])
    rownames(alpha.matrix) <- my.sp
    colnames(alpha.matrix) <- my.sp
    
    # double check
    alpha.matrix[is.na(alpha.matrix)] <- 0
    
    # calculate feasibility of all combinations
    if(length(my.sp)>1){
      for(nsp in 2:length(my.sp)){
        # draw combinations of nsp size
        combos <- t(combn(my.sp,nsp))
        
        # store the species present in every combination
        # for, later, the accumulation curves
        allspmatrix <- matrix(nrow = nrow(combos),ncol = length(all.sp))
        colnames(allspmatrix) <- all.sp
        pw_st <- as.data.frame(allspmatrix)
        pw_st$year <- years[i.year]
        pw_st$plot <- i.plot
        pw_st$max_eigenvalue <- NA_real_
        
        # for each combination, obtain maximum eigenvalue of the jacobian
        # at steady state
        for(s in 1:nrow(pw_st)){
          comb.species <- combos[s,]

          my.lambdas <- lambda.plot$lambda.log[lambda.plot$sp %in% comb.species]
          
          my.g <- sp.rates$germination.rate[sp.rates$species.code %in% comb.species]
          my.s <- sp.rates$seed.survival[sp.rates$species.code %in% comb.species]
          
          # combo matrix
          # as.matrix() for the case of single species combos
          A <- as.matrix(alpha.matrix[comb.species,comb.species])
          
          # just in case
          A[is.na(A)] <- 0
          
          # check conditions for calculating feasibility
          valid.combo <- FALSE
          
          # first, for pairs, species need to interact
          if(nsp == 2){
            offdiag <- A[row(A)!=col(A)]
            if(sum(offdiag != 0)>0){
              valid.combo <- TRUE
            }
            # second, for >2 sp, matrix must be invertible (det != 0) 
          }else{
            if(det(A) != 0){
              valid.combo <- TRUE
            }
          }
          
          if(valid.combo){
            
            sp.r <- growth_rate_equivalence(lambda = my.lambdas,
                                            g = my.g,
                                            s = my.s,sign = "positive")
            
            # local stability
            # initial abundances
            comb.abund <- subset(abundances, year == years[i.year] & 
                                   plot == i.plot &
                                   species %in% comb.species)
            
            init.abund <- comb.abund %>% 
              group_by(species) %>% 
              summarise(abund = sum(individuals))
            
            init.abund <- init.abund$abund
            
            # Steady state solution
            # the commented code is for continuous systems
            # eq.abund = rootSolve::stode(y = init.abund,
            #                             func=model_stode,
            #                             parms=NULL,
            #                             A=A,
            #                             gr=sp.r, 
            #                             positive = TRUE)[[1]]
            
            # but this is a discrete time one
            eq.abund = deSolve::ode(y = init.abund,
                                        func=model_stode,
                                        parms=NULL,
                                        A=A,
                                        gr=sp.r,
                                    times = 0:1e3,
                                    method = "iteration") 
            
            # Jacobian
            J = rootSolve::jacobian.full(y=eq.abund[nrow(eq.abund),2:(ncol(eq.abund)-1)],
                                         func=model_J,A=A,r=sp.r)
            
            # Local stability analysis
            if(any(is.na(J)) | any(is.infinite(J))){
              pw_st$max_eigenvalue[s] <- NA_real_
            }else{
              # in discrete-time systems, no eigenvalue can be < -1 or > 1
              # for the system to be locally stable
              # this is opposed to continuous-time systems, where the condition
              # is <= 0
              pw_st$max_eigenvalue[s] <- max(abs(as.double(eigen(J)$values)))	
            }

            # which sp
            pw_st[s,comb.species] <- TRUE
            pw_st[s,all.sp[which(!all.sp %in% comb.species)]] <- FALSE
            
          }# if valid combination
        }# for each sp combination
        
        # inefficient, but enough
        local.stab.list[[length(local.stab.list)+1]] <- pw_st
        
        # store individual iterations?
        # write.csv2(pw_st,paste("results/local_stability_analysis/S5_local_stability_",
        #                        fit.type,"_",years[i.year],"_",i.plot,".csv",sep=""),
        #            row.names = FALSE)
        
      }# for combinations of nsp
    }# if >1 my.sp
    cat(date()," - year:",years[i.year],", plot:",i.plot," completed\n",sep="")
  }# for i.plot
}# for i.year

local.stab.df <- bind_rows(local.stab.list)
validrows <- rowSums(local.stab.df[,all.sp])
local.stab.df.clean <- local.stab.df[which(!is.na(validrows)),]

write.csv2(local.stab.df.clean,paste("results/local_stability_analysis/S5_local_stability_",fit.type,sep=""),row.names = FALSE)

