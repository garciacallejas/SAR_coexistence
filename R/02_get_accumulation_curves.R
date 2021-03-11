
# in this file, we combine info on richness and coexistence to
# plot accumulation curves

# the homogeneous coexistence curve assumes that only one community
# coexists per plot

# the heterogeneous coexistence curve assumes that more than one
# communities (modules) may coexist per plot, thus implying spatial
# segregation

# these are calculated in different ways. Here we calculate both curves
# and save them to disk.

# read data ---------------------------------------------------------------
library(tidyverse)

structural_metrics <- read.csv2(file = "./results/02_observed_structural_metrics.csv")

years <- sort(unique(structural_metrics$year))
all.plots <- sort(unique(structural_metrics$plot))
n.plots <- length(all.plots)
all.sp <- names(structural_metrics)[which(!names(structural_metrics) 
                                          %in% c("year",
                                                 "plot",
                                                 "feasibility",
                                                 "feasibility.domain",
                                                 "com.pair.differential",
                                                 "fit"))]

# only feasible communities
feasna <- rowSums(structural_metrics[,all.sp])
feasvalid <- structural_metrics[which(!is.na(feasna)),]
feasvalid <- subset(feasvalid, feasibility == 1)

feashet <- subset(feasvalid,fit == "heterogeneous_space")
feashom <- subset(feasvalid,fit == "homogeneous_space")

# 1 - non-additive homogeneous curve --------------------------------

feashomax <- NULL
for(i.year in 1:length(years)){
  for(i.plot in 1:n.plots){
    # for(i.type in 1:length(param.types)){
    feasdata <- feashom[which(feashom$year == years[i.year] &
                                feashom$plot == i.plot),]
    # feasdata <- subset(feasres,year == years[i.year] & 
    #                      plot == i.plot &
    #                      fit.type == param.types[i.type])
    # feasdata <- subset(feasdata,feasibility == 1)
    if(nrow(feasdata)>0){
      maxrich <- max(rowSums(feasdata[,all.sp]))#feasdata$n.sp
      frows <- which(rowSums(feasdata[,all.sp]) == maxrich)
      # feasdata <- feasdata[frows,]
      feashomax <- rbind(feashomax,feasdata[frows,])
    }
    # }# for i.type
  }# for i.plot
}# for i.year


# for each year, select combinatios of 1:n plots
# and calculate explicitly the number of coexisting sp
# in each combination
# note that a given plot may have two or more combinations
# that are feasible with the same number of species 
# but different identities
# homcoex <- NULL
homcoex <- list()

for(i.year in 1:length(years)){
  
  year.data <- subset(feashomax, year == years[i.year])
  yearcoex <- NULL
  
  for(n.plots in 1:max(all.plots)){
    nplotscoex <- NULL
    # select n plots
    allcombos <- combn(all.plots,n.plots)
    # how many combinations of n.plots
    ncombs <- ncol(allcombos)
    # for each of these combinations
    for(i.comb in 1:ncombs){
      my.plots <- allcombos[,i.comb]
      my.data <- subset(year.data, plot %in% my.plots)
      if(nrow(my.data)>0){
        combsperplot <- my.data %>% group_by(plot) %>% summarise(num = n())
        # translate to list, in order to use expand.grid
        combslist <- list()
        for(ip in 1:nrow(combsperplot)){
          combslist[[ip]] <- 1:combsperplot$num[ip]
        }
        # mycombs is a dataframe in which rows are 
        # combinations from each plot (e.g. the nth feasible community)
        # and columns are plots
        
        mycombs <- expand.grid(combslist)
        names(mycombs) <- combsperplot$plot
        
        # generate results dataframe
        combres <- expand.grid(#fit.type = param.types[i.type],
          year = years[i.year],
          n.plots = n.plots,
          plot.comb = i.comb,
          sp.comb = 1:nrow(mycombs),
          n.sp = 0,stringsAsFactors = FALSE)
        
        # count the combinations in each plot
        my.data <- my.data %>% 
          group_by(plot) %>% 
          mutate(cc = row_number())
        
        for(ipl in 1:nrow(mycombs)){
          
          # select the comb from the ipl-th row in mycombs
          # this needs some crafting, not intuitive
          
          # if(ncol(mycombs)>1){
          spcomb <- my.data %>%  
            group_by(plot) %>% 
            mutate(pid = group_indices()) %>%
            filter(cc == mycombs[ipl,pid]) %>%
            # and the easy part,
            # extract richness for this combination
            gather(key = "sp",value = "comb",all_of(all.sp)) %>%
            group_by(sp) %>%
            summarise(present = sum(comb)) %>%
            summarise(richness = sum(present > 0))
          
          combres$n.sp[ipl] <- as.numeric(spcomb)
        }
        
        combm <- combres %>% group_by(#fit.type,
          year,
          n.plots,
          plot.comb) %>% summarise(mean.sp = mean(n.sp),
                                   min.sp = min(n.sp),
                                   max.sp = max(n.sp))
        nplotscoex <- dplyr::bind_rows(nplotscoex,combm)
      }# if my.data > 0
    }# for i.comb
    yearcoex <- dplyr::bind_rows(yearcoex,nplotscoex)
  }# for n.plots
  # homcoex <- dplyr::bind_rows(homcoex,yearcoex)
  homcoex[[i.year]] <- yearcoex
  print(paste("year",i.year,"ok"))
}# for i.year

homcoexdf <- bind_rows(homcoex)

# 2 - additive heterogeneous curve --------------------------------------------

fun.na <- function(x){all(is.na(x))}

# set of coexisting sp per plot and year
coex.plot.year <- expand.grid(years,all.plots)
names(coex.plot.year) <- c("year","plot")
coex.plot.year[,all.sp] <- FALSE

for(i.year in 1:length(years)){
  for(i.plot in all.plots){
    
    # species in feasible combinations
    feas.comb <- subset(feashet, year == years[i.year] & plot == i.plot)
    
    # these combinations are already robust,
    # in that combos with all NAs and any negative values are excluded 
    
    if(nrow(feas.comb)>0){
      comb.sp <- colSums(feas.comb[,all.sp])
      stable.sp <- names(comb.sp[which(comb.sp > 0)])
      coex.plot.year[coex.plot.year$year == years[i.year] & 
                       coex.plot.year$plot == i.plot,stable.sp] <- TRUE
      
    }# if any stable combination
  }# for plot
}# for year

area.results <- NULL

for(i.num in 1:n.plots){
  my.combinations <- combn(1:9,i.num)
  for(i.comb in 1:ncol(my.combinations)){
    my.plots <- as.character(my.combinations[,i.comb])
    # concatenate into a single string
    my.plots.char <- paste(my.plots,collapse = "_")
    for(i.year in 1:length(years)){
      my.data <- subset(coex.plot.year, year == years[i.year] & plot %in% my.plots)
      
      # how many coexisting sp in this combination of plots
      # identity of sp does not matter
      my.comb.sp <- sum(colSums(my.data[,all.sp]) != 0)
      
      combres <- data.frame(fit.type = "heterogeneous",#param.types[i.type],
                            year = years[i.year],
                            n.plots = i.num,
                            plots = my.plots.char,
                            num.sp = my.comb.sp,stringsAsFactors = F)
      
      area.results <- bind_rows(area.results,combres)
    }# for i.year
  }# for i.comb
}# for i.num

# store curves ------------------------------------------------------------

write.csv2(area.results,"./results/03_heterogeneous_curve.csv",row.names = F)
write.csv2(homcoexdf,"./results/03_07_homogeneous_curve.csv",row.names = F)

