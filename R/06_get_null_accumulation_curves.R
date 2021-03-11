
# parallelized code for obtaining accumulation curves
# As with script 05, be aware that this may be computationally intensive.

# load packages -----------------------------------------------------------

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 10
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# read replicates ---------------------------------------------------------

# only need structural metrics

files <- list.files(path = "./results/null_replicates/")
str.files <- files[grep("null_structural_metrics",files)]

null_metrics <- NULL

for(i.file in 1:length(str.files)){
    my.file <- read.csv2(file = paste("./results/null_replicates",str.files[i.file],sep=""),
                         header = TRUE,
                         stringsAsFactors = FALSE)
  null_metrics <- dplyr::bind_rows(null_metrics,my.file)
}

id.char <- unique(paste(null_metrics$null_model,"_",null_metrics$year,"_",null_metrics$replicate,sep=""))

# calculate accumulation curves -------------------------------------------

null.plots <- sort(unique(null_metrics$plot))
n.null.plots <- length(null.plots)
null.sp <- names(null_metrics)[which(!names(null_metrics)
                                     %in% c("year",
                                            "plot",
                                            "feasibility",
                                            "feasibility.domain",
                                            "com.pair.differential",
                                            "replicate",
                                            "null_model"))]

# null curve --------------------------------------------------------------
# obtain the accumulation curve for every rep:year:model combination
# and store them

# this is parallelized
nullcoex <- foreach(i.id = 1:length(id.char), 
                    .combine=rbind, 
                    .packages = 'tidyverse') %dopar% {

  # recover year, model, and replicate from the ID

i.year <- substr(id.char[i.id],4,7)
i.model <- substr(id.char[i.id],1,2)
i.rep <- substr(id.char[i.id],9,nchar(id.char[i.id]))

  repdata <- null_metrics[null_metrics$replicate == i.rep &
                            null_metrics$year == i.year &
                            null_metrics$null_model == i.model, ]

  repmax <- NULL
  for(i.plot in 1:n.null.plots){
    feasdata <- repdata[which(repdata$plot == i.plot),]
    if(nrow(feasdata)>0){
      maxrich <- max(rowSums(feasdata[,null.sp]))
      frows <- which(rowSums(feasdata[,null.sp]) == maxrich)
      repmax <- rbind(repmax,feasdata[frows,])
    }
  }# for i.plot


  # select combinatios of 1:n plots
  # and calculate explicitly the number of coexisting sp
  # in each combination
  # note that a given plot may have two or more combinations
  # that are feasible with the same number of species
  # but different identities
  
  repcoex <- NULL

  for(i.plots in 1:n.null.plots){
    nplotscoex <- list()
    # select n plots
    allcombos <- combn(n.null.plots,i.plots)
    # how many combinations of n.plots
    ncombs <- ncol(allcombos)
    # for each of these combinations
    for(i.comb in 1:ncombs){
      my.plots <- allcombos[,i.comb]
      my.data <- subset(repmax, plot %in% my.plots)
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
        combres <- expand.grid(
          n.plots = i.plots,
          plot.comb = i.comb,
          sp.comb = 1:nrow(mycombs),
          n.sp = 0,stringsAsFactors = FALSE)

        # count the combinations in each plot
        my.data <- my.data %>%
          group_by(plot) %>%
          mutate(cc = row_number())

        for(ipl in 1:nrow(mycombs)){

          testcomb <- NULL
          for(ip in 1:ncol(mycombs)){
            testcomb <- rbind(testcomb,my.data[my.data$plot == names(mycombs)[ip] &
                                                 my.data$cc == mycombs[ipl,ip],])
          }

          spr <- testcomb %>% gather(key = "sp",value = "comb",all_of(null.sp)) %>%
            group_by(sp) %>%
            summarise(present = sum(comb)) %>%
            summarise(richness = sum(present > 0))

          combres$n.sp[ipl] <- as.numeric(spr)
        }

        combm <- combres %>% group_by(
          n.plots,
          plot.comb) %>% summarise(mean.sp = mean(n.sp),
                                   min.sp = min(n.sp),
                                   max.sp = max(n.sp))
        nplotscoex[[i.comb]] <- combm
      }# if my.data > 0
    }# for i.comb

    np <- bind_rows(nplotscoex)
    repcoex <- bind_rows(repcoex,np)
  }# for n.plots

  repcoex$year <- i.year
  repcoex$null_model <- i.model
  repcoex$replicate <- i.rep
  repcoex
#
}# for i.id

# # write curves to disk ----------------------------------------------------

write.csv2(nullcoex,
           file = paste(getwd(),"/results/07_null_accumulation_curves.csv",sep=""),
           row.names = FALSE)

# clear cores -------------------------------------------------------------

stopCluster(cl)
          