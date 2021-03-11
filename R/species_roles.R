#' obtain species coexistence roles
#'
#' given a set of communities, which sp enter by pairwise coex, indirect effects,
#' or are transient/dominant. See `structural_metrics` for calculating the communities.
#'
#' @param structural_metrics output from `structural_metrics` function.
#' @param lambdas data frame, year, plot, species, lambda.
#' @param sp.rates species germination and seed survival rates.
#' @param abundances data frame, year, plot, subplot, species, individuals. (Carcoles/data/abundances).
#'
#' @return species roles dataframe
#' @export
#'
#' @examples
species_roles <- function(structural_metrics,lambdas,sp.rates,abundances){

  all.sp <- unique(lambdas$sp)#sort(unique(sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]))
  years <- unique(lambdas$year)
  all.plots <- sort(unique(lambdas$plot))
  n.plots <- length(all.plots)
  
  # add vital rates to lambdas
  lambdas$germ.rate <- sp.rates$germination.rate[match(lambdas$sp,sp.rates$species.code)]
  lambdas$surv.rate <- sp.rates$seed.survival[match(lambdas$sp,sp.rates$species.code)]
  
  # transform lambdas from the negbin to r in a linear LV 
  growth_rate_equivalence <- function(lambda,g,s,sign = "negative"){
    if(sign == "negative"){
      r <- log((1-(1-g)*s)/g) - lambda
    }else{
      r <- log((1-(1-g)*s)/g) + lambda
    }
    r
  }
  
  # see comment in "structural_metrics" function
  lambdas$lambda.log <- log(lambdas$lambda)
  
  lambdas$gr <- growth_rate_equivalence(lambda = lambdas$lambda.log,
                                        g = lambdas$germ.rate,
                                        s = lambdas$surv.rate,sign = "positive")
  
  # hetfeas <- subset(feasres.clean,fit.type == "heterogeneous_both")
  
  abplot <- abundances %>% filter(species %in% all.sp) %>% group_by(year,plot,species) %>% summarise(num = sum(individuals))
  
  sproles <- expand.grid(year = years,plot = all.plots,sp = all.sp)
  sproles$pairwise <- FALSE
  sproles$indirect <- FALSE
  sproles$dominant <- FALSE
  sproles$transient <- FALSE
  
  # feasible and not feasible combinations are checked, to
  # assess transient and dominant sp.
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:n.plots){
      for(i.sp in 1:length(all.sp)){
        my.pos <- which(sproles$year == years[i.year] &
                          sproles$plot == i.plot &
                          sproles$sp == all.sp[i.sp])
        
        # feasibility values for this site and with this sp
        # present
        sitefeas <- subset(structural_metrics,year == years[i.year] & 
                             plot == i.plot)
        mycol <- which(names(sitefeas) == all.sp[i.sp])
        spfeas <- sitefeas[which(sitefeas[,mycol] == TRUE),]
        abpresent <- abplot$species[abplot$year == years[i.year] & 
                                      abplot$plot == i.plot & 
                                      abplot$num > 0]
        
        # if there are combinations with i.sp, check them
        # otherwise, if no combinations but the sp is present in abpresent,
        # for now we treat it as transient
        if(nrow(spfeas)>0){
          spfeas$presentsp <- rowSums(spfeas[,all.sp])
          # nas to zeroes
          spfeas$feasibility[which(is.na(spfeas$feasibility))] <- 0
          # pair combinations
          pairfeas <- subset(spfeas,presentsp == 2)
          # n-sp combinations
          netfeas <- subset(spfeas,presentsp > 2)
          # any feasible pair?
          if(sum(pairfeas$feasibility == 1) > 0){
            sproles$pairwise[my.pos] <- TRUE
          }
          # any feasible network?
          if(sum(netfeas$feasibility == 1) > 0){
            sproles$indirect[my.pos] <- TRUE
          }
          # if no pairs or network,
          # is it dominant in this connected component?
          if(sproles$pairwise[my.pos] == FALSE &
             sproles$indirect[my.pos] == FALSE){
            
            # with which species does it "share" the network?
            intspt <- colSums(spfeas[,all.sp])
            intsp <- names(intspt[intspt>0])
            
            # not this sp
            intsp <- intsp[which(intsp != all.sp[i.sp])]
            
            # compare growth rates in this plot/year
            mygr <- lambdas[lambdas$year == years[i.year] &
                              lambdas$plot == i.plot,c("sp","gr")]
            # growth rate of focal species
            myspgr <- mygr$gr[mygr$sp == all.sp[i.sp]]
            # growth rates of interacting sp
            mygr <- subset(mygr,!is.na(gr) & sp %in% intsp)
            # if no growth rate is higher, focal sp is dominant
            if(!is.na(myspgr)){
              if(sum(mygr$gr > myspgr) == 0){
                sproles$dominant[my.pos] <- TRUE
              }
            }
          }# if no pairs or indirect
          
          # if no pairs, indirect, or dominant, 
          # it is transient
          if(!sproles$pairwise[my.pos] & 
             !sproles$indirect[my.pos] &
             !sproles$dominant[my.pos]){
            sproles$transient[my.pos] <- TRUE
          }
          # if sp not in any combination but present in original data,
          # mark it as transient
        }else if(all.sp[i.sp] %in% abpresent){
          sproles$transient[my.pos] <- TRUE
        }# if-else sp present in this plot/year
      }# for i.sp
    }# for i.plot
  }# for i.year
  
  sproles
  # write.csv2(sproles,file = "./results/structural_metrics_sproles.csv",row.names = FALSE)
  
}