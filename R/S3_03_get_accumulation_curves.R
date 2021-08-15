
# obtain coexistence-area curves from a series of replicated communities

# load packages -----------------------------------------------------------
options(tidyverse.quiet = TRUE)

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# read replicates ---------------------------------------------------------

# some years do not have all species, in order to append them I need all columns
valid.sp <- c("BEMA","CETE","CHFU","CHMI","HOMA","LEMA","LYTR","MEEL","MESU", 
              "PAIN","PLCO","POMA","POMO","PUPA","SASO","SCLA","SOAS","SPRU","SUSP")

# only need structural metrics
files <- list.files(path = "./results/posterior_draws/")
str.files <- files[grep("S3_posterior_draws_structural_metrics",files)]

null_metrics <- NULL

for(i.file in 1:length(str.files)){
  my.file <- read.csv2(file = paste("./results/posterior_draws/",str.files[i.file],sep=""),
                       header = TRUE,
                       stringsAsFactors = FALSE)
  
  mysp <- names(my.file)[which(!names(my.file)
                               %in% c("year",
                                      "plot",
                                      "feasibility",
                                      "feasibility.domain",
                                      "com.pair.differential",
                                      "replicate",
                                      "fit.type"))]
  miss <- valid.sp[which(!valid.sp %in% mysp)]
  my.file[,miss] <- FALSE
  
  null_metrics <- dplyr::bind_rows(null_metrics,my.file)
}

id.char <- unique(paste(null_metrics$fit.type,"_",null_metrics$year,"_",null_metrics$replicate,sep=""))

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
                                            "fit.type"))]

# null curve --------------------------------------------------------------
# obtain the accumulation curve for every rep:year:model combination
# and store them

cat(format(Sys.time(),usetz = TRUE)," - STARTING with",workers,"cores")

# this is parallelized
nullcoex <- foreach(i.id = 1:length(id.char), 
                    .packages = 'tidyverse') %dopar% {
                      
                      # recover year, type, and replicate from the ID
                      
                      i.type <- substr(id.char[i.id],1,18)
                      i.year <- substr(id.char[i.id],20,23)
                      i.rep <- substr(id.char[i.id],25,nchar(id.char[i.id]))
                      
                      repdata <- null_metrics[which(null_metrics$replicate == i.rep &
                                                      null_metrics$year == i.year &
                                                      null_metrics$fit.type == i.type),]
                      
                      # the accumulation curves are computed differently
                      # for homogeneous and heterogeneous parameterizations
                      if(i.type == "heterogeneous_both"){
                        
                        # set of coexisting sp per plot and year
                        coex.plot.year <- data.frame(year = i.year, plot = 1:n.null.plots)
                        coex.plot.year[,null.sp] <- FALSE
                        
                        for(i.plot in null.plots){
                          
                          # species in feasible combinations
                          feas.comb <- subset(repdata, plot == i.plot)
                          
                          # these combinations are already robust,
                          # in that combos with all NAs and any negative values are excluded 
                          
                          if(nrow(feas.comb)>0){
                            comb.sp <- colSums(feas.comb[,null.sp],na.rm = TRUE)
                            stable.sp <- names(comb.sp[which(comb.sp > 0)])
                            coex.plot.year[coex.plot.year$plot == i.plot,stable.sp] <- TRUE
                            
                          }# if any stable combination
                        }# for plot
                        
                        repcoex <- NULL
                        
                        for(i.num in 1:n.null.plots){
                          my.combinations <- combn(1:9,i.num)
                          for(i.comb in 1:ncol(my.combinations)){
                            my.plots <- as.character(my.combinations[,i.comb])
                            # concatenate into a single string
                            my.plots.char <- paste(my.plots,collapse = "_")
                            # for(i.year in 1:length(years)){
                            my.data <- subset(coex.plot.year, plot %in% my.plots)
                              
                              # how many coexisting sp in this combination of plots
                              # identity of sp does not matter
                              my.comb.sp <- sum(colSums(my.data[,null.sp]) != 0)
                              
                              combres <- data.frame(n.plots = i.num,
                                                    plots = my.plots.char,
                                                    num.sp = my.comb.sp,stringsAsFactors = F)
                              
                              repcoex <- bind_rows(repcoex,combres)
                            # }# for i.year
                          }# for i.comb
                        }# for i.num
                        
                        repcoex$year <- i.year
                        repcoex$fit.type <- i.type
                        repcoex$replicate <- i.rep
                        # return
                        repcoex
                      }else{
                        
                        repmax <- NULL
                        for(i.plot in 1:n.null.plots){
                          feasdata <- repdata[which(repdata$plot == i.plot),]
                          feasdata[is.na(feasdata)] <- FALSE 
                          
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
                                
                                testcomb.list <- list()
                                for(ip in 1:ncol(mycombs)){
                                  testcomb.list[[ip]] <- my.data[my.data$plot == names(mycombs)[ip] &
                                                                       my.data$cc == mycombs[ipl,ip],]
                                }
                                
                                testcomb <- bind_rows(testcomb.list)
                                
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
                        repcoex$fit.type <- i.type
                        repcoex$replicate <- i.rep
                        # return
                        repcoex
                        
                      }# if-else heterogeneous_both or heterogeneous_time
                    }# for i.id

cat("\n",format(Sys.time(),usetz = TRUE)," - FINISHED with",workers,"cores")

# # write curves to disk ----------------------------------------------------
hom.curves <- NULL
het.curves <- NULL

for(i.id in 1:length(nullcoex)){
  if(nullcoex[[i.id]]$fit.type[1] == "heterogeneous_both"){
    het.curves <- bind_rows(het.curves,nullcoex[[i.id]])
  }else{
    hom.curves <- bind_rows(hom.curves,nullcoex[[i.id]])
  }
}

if(!is.null(hom.curves)){
  write.csv2(hom.curves,
             file = paste(getwd(),"results/posterior_draws/S3_posterior_draws_accumulation_curves_heterogeneous_time.csv",sep=""),
             row.names = FALSE)
}

if(!is.null(het.curves)){
  write.csv2(het.curves,
             file = paste(getwd(),"results/posterior_draws/S3_posterior_draws_accumulation_curves_heterogeneous_both.csv",sep=""),
             row.names = FALSE)
}
# clear cores -------------------------------------------------------------

stopCluster(cl)
