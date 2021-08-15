#' obtain feasibility from model parameters
#' 
#' @param lambdas data frame, year, plot, sp, lambda 
#' @param alphas data frame, year, plot, focal, neighbour, magnitude
#' @param sp.rates data frame, species.code, germination.rate, seed.survival 
#'
#' @return data frame with feasibility of every combination
#' @export
#'
#' @examples
structural_metrics_lite <- function(lambdas, alphas, sp.rates){

  ##################
  # aux functions
  # source("./R/toolbox_coexistence.R")

  # transform lambdas from the negbin to r in a linear LV 
  growth_rate_equivalence <- function(lambda,g,s,sign = "negative"){
    if(sign == "negative"){
      r <- log((1-(1-g)*s)/g) - lambda
    }else{
      r <- log((1-(1-g)*s)/g) + lambda
    }
    r
  }
  ###################
  
  years <- unique(lambdas$year)
  n.plots <- length(unique(lambdas$plot))
  all.sp <- unique(lambdas$sp)
  
  ####################### IMPORTANT NOTE
  # alphas as read here have signs switched,
  # so competition are alphas > 0
  # this is important later for calculating the
  # growth rate equivalence between lambdas and growth rates of linear LV

  # lambdas are stored in "real" dimensions, i.e. number of seeds
  # but the linear LV works with model raw parameters
  # so here, undo the transformation by taking log again
  #######################-
  
  lambdas$lambda.log <- log(lambdas$lambda)
  
  feasres <- list()
  
  for(i.year in 1:length(years)){
    
    lambda.year <- subset(lambdas,year == years[i.year])
    alpha.year <- subset(alphas,year == years[i.year])
    n.plots <- max(c(unique(lambda.year$plot),unique(alpha.year$plot)))
    
    for(i.plot in 1:n.plots){
      
      lambda.plot <- subset(lambda.year,plot == i.plot)
      alpha.plot <- subset(alpha.year,plot == i.plot)
      
      # present species
      lambda.plot <- subset(lambda.plot, !is.na(lambda))
      my.sp <- sort(unique(intersect(lambda.plot$sp,alpha.plot$focal)))
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
          # pw_st$fit.type <- param.types[i.type]
          pw_st$feasibility <- 0
          pw_st$feasibility.domain <- 0
          pw_st$com.pair.differential <- 0
          
          # for each combination, obtain feasibility
          for(s in 1:nrow(pw_st)){
            species <- combos[s,]
            
            # my.lambdas <- lambda.plot$lambda[lambda.plot$sp %in% species]
            my.lambdas <- lambda.plot$lambda.log[lambda.plot$sp %in% species]
            
            my.g <- sp.rates$germination.rate[sp.rates$species.code %in% species]
            my.s <- sp.rates$seed.survival[sp.rates$species.code %in% species]
            
            # combo matrix
            # as.matrix() for the case of single species combos
            A <- as.matrix(alpha.matrix[species,species])
            
            # # just in case
            # A[is.na(A)] <- 0
            
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
              # feasibility
              feas <- try(test_feasibility(alpha = A,
                                           r = sp.r),
                          silent = TRUE)
              
              # if feasibility was computed, append it
              # domain <- NA_real_
              # cpdiff <- NA_real_
              # if(!inherits(feas,"try-error")){				
              #   fmetrics <- try(compute_overlap(A,1e3))
              #   if(class(fmetrics) != "try-error"){
              #     domain <- fmetrics$Omega
              #     cpdiff <- fmetrics$Omega - fmetrics$Omega_all
              #   }
              # }
              
              pw_st$feasibility[s] <- ifelse(inherits(feas,"try-error"),NA_real_,feas)
              # pw_st$feasibility.domain[s] <- domain
              # pw_st$com.pair.differential[s] <- cpdiff
              
              # which sp
              pw_st[s,species] <- TRUE
              pw_st[s,all.sp[which(!all.sp %in% species)]] <- FALSE
              
            }# if valid combination
          }# for each sp combination
          
          # terribly inefficient, but that's life
          feasres[[length(feasres) + 1]] <- pw_st
          #<- dplyr::bind_rows(feasres,pw_st)
          
        }# for combinations of nsp
      }# if >1 my.sp
    }# for i.plot
  }# for i.year

  feasres.all <- bind_rows(feasres)
  validrows <- rowSums(feasres.all[,all.sp])
  feasres.clean <- feasres.all[which(!is.na(validrows)),]

  feasres.clean

}