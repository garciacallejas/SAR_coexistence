
# obtain coexistence-area curves from a series of replicated communities
# this script is very demanding computationally
# we include it for reference: for the actual calculations, we adapted it
# to run it in a HPC at Cadiz University

# load packages -----------------------------------------------------------
options(tidyverse.quiet = TRUE)

library(foreach)
library(doParallel)
library(tidyverse)

# set number of cores -----------------------------------------------------

workers <- 8
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# read replicates ---------------------------------------------------------

# some years do not have all species, in order to append them I need all columns
valid.sp <- c("BEMA","CETE","CHFU","CHMI","HOMA","LEMA","LYTR","MEEL","MESU", 
              "PAIN","PLCO","POMA","POMO","PUPA","SASO","SCLA","SOAS","SPRU","SUSP")

# only need structural metrics
files <- list.files(path = "./results/rate_sensitivity_analysis/")
str.files <- files[grep("S4_modified_rates_structural_metrics",files)]

rates_metrics <- NULL

for(i.file in 1:length(str.files)){
  my.file <- read.csv2(file = paste("./results/rate_sensitivity_analysis/",str.files[i.file],sep=""),
                       header = TRUE,
                       stringsAsFactors = FALSE)
  
  mysp <- names(my.file)[which(!names(my.file)
                               %in% c("year",
                                      "plot",
                                      "feasibility",
                                      "feasibility.domain",
                                      "com.pair.differential",
                                      "replicate",
                                      "analysis",
                                      "fit.type"))]
  miss <- valid.sp[which(!valid.sp %in% mysp)]
  my.file[,miss] <- FALSE
  
  rates_metrics <- dplyr::bind_rows(rates_metrics,my.file)
}

# -------------------------------------------------------------------------
# NOTE
# we calculate increases and decreases in both rates for the 
# heterogeneous parameterization, 
# as well as random variation; but the possibility is open
# to calculate other modified curves. We only calculate these 
# due to the computational effort needed, and 
# because if these show no clear trend, it is expected that 
# the same calculation for the homogeneous curves
# will also show no clear trend.

rates_metrics <- subset(rates_metrics,fit.type == "heterogeneous_both" &
                          analysis %in% c("germ_increase","germ_decrease", "germ_random",
                                          "surv_increase","surv_decrease", "surv_random"))

# -------------------------------------------------------------------------
# NOTE: change the number of leading zeros if needed, depending on the number
# of replicates
id.char <- unique(paste(rates_metrics$fit.type,"_",
                        rates_metrics$year,"_",
                        formatC(rates_metrics$replicate, width = 2, format = "d", flag = "0"),"_",
                        rates_metrics$analysis,
                        sep=""))

# calculate accumulation curves -------------------------------------------

rates.plots <- sort(unique(rates_metrics$plot))
n.rates.plots <- length(rates.plots)
rates.sp <- names(rates_metrics)[which(!names(rates_metrics)
                                     %in% c("year",
                                            "plot",
                                            "feasibility",
                                            "feasibility.domain",
                                            "com.pair.differential",
                                            "replicate",
                                            "analysis",
                                            "fit.type"))]

# null curve --------------------------------------------------------------
# obtain the accumulation curve for every rep:year:type:analysis combination
# and store them

cat(format(Sys.time(),usetz = TRUE)," - STARTING with",workers,"cores")

# this is parallelized
ratescoex <- foreach(i.id = 1:length(id.char), 
                    .packages = 'tidyverse') %dopar% {
                      
                      # recover year, type, and replicate from the ID
                      
                      i.type <- substr(id.char[i.id],1,18)
                      i.year <- substr(id.char[i.id],20,23)
                      i.rep <- as.numeric(substr(id.char[i.id],25,26))
                      i.analysis <- substr(id.char[i.id],28,nchar(id.char[i.id]))
                      
                      repdata <- rates_metrics[which(rates_metrics$replicate == i.rep &
                                                       rates_metrics$analysis == i.analysis &
                                                       rates_metrics$year == i.year &
                                                       rates_metrics$fit.type == i.type),]
                      
                      # the accumulation curves are computed differently
                      # for homogeneous and heterogeneous parameterizations
                      if(i.type == "heterogeneous_both"){
                        
                        # set of coexisting sp per plot and year
                        coex.plot.year <- data.frame(year = i.year, plot = 1:n.rates.plots)
                        coex.plot.year[,rates.sp] <- FALSE
                        
                        for(i.plot in rates.plots){
                          
                          # species in feasible combinations
                          feas.comb <- subset(repdata, plot == i.plot)
                          
                          # these combinations are already robust,
                          # in that combos with all NAs and any negative values are excluded 
                          
                          if(nrow(feas.comb)>0){
                            comb.sp <- colSums(feas.comb[,rates.sp],na.rm = TRUE)
                            stable.sp <- names(comb.sp[which(comb.sp > 0)])
                            coex.plot.year[coex.plot.year$plot == i.plot,stable.sp] <- TRUE
                            
                          }# if any stable combination
                        }# for plot
                        
                        repcoex <- NULL
                        
                        for(i.num in 1:n.rates.plots){
                          my.combinations <- combn(1:9,i.num)
                          for(i.comb in 1:ncol(my.combinations)){
                            my.plots <- as.character(my.combinations[,i.comb])
                            # concatenate into a single string
                            my.plots.char <- paste(my.plots,collapse = "_")
                            # for(i.year in 1:length(years)){
                            my.data <- subset(coex.plot.year, plot %in% my.plots)
                              
                              # how many coexisting sp in this combination of plots
                              # identity of sp does not matter
                              my.comb.sp <- sum(colSums(my.data[,rates.sp]) != 0)
                              
                              combres <- data.frame(n.plots = i.num,
                                                    plots = my.plots.char,
                                                    num.sp = my.comb.sp,stringsAsFactors = F)
                              
                              repcoex <- bind_rows(repcoex,combres)
                            # }# for i.year
                          }# for i.comb
                        }# for i.num
                        
                        repcoex$year <- i.year
                        repcoex$fit.type <- i.type
                        repcoex$analysis <- i.analysis
                        repcoex$replicate <- i.rep
                        # return
                        repcoex
                      }else{
                        
                        repmax <- NULL
                        for(i.plot in 1:n.rates.plots){
                          feasdata <- repdata[which(repdata$plot == i.plot),]
                          feasdata[is.na(feasdata)] <- FALSE 
                          
                          if(nrow(feasdata)>0){
                            maxrich <- max(rowSums(feasdata[,rates.sp]))
                            frows <- which(rowSums(feasdata[,rates.sp]) == maxrich)
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
                        
                        for(i.plots in 1:n.rates.plots){
                          nplotscoex <- list()
                          # select n plots
                          allcombos <- combn(n.rates.plots,i.plots)
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
                                  testcomb.list[[ip]] <- subset(my.data,plot == names(mycombs)[ip] &
                                                                  cc == mycombs[ipl,ip])
                                }
                                
                                testcomb <- bind_rows(testcomb.list)
                                
                                spr <- testcomb %>% gather(key = "sp",value = "comb",all_of(rates.sp)) %>%
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
                        repcoex$analysis <- i.analysis
                        repcoex$replicate <- i.rep
                        
                        write_csv2(repcoex,file = paste("rate_sensitivity_analysis/S4_modified_rates_accumulation_curves_",
                                                        i.year,"_",
                                                        i.type,"_",
                                                        i.analysis,"_",
                                                        i.rep,".csv",sep=""))
                        
                        # return
                        repcoex
                        
                      }# if-else heterogeneous_both or heterogeneous_time
                    }# for i.id

cat("\n",format(Sys.time(),usetz = TRUE)," - FINISHED with",workers,"cores")

# # write curves to disk ----------------------------------------------------
hom.curves <- NULL
het.curves <- NULL

for(i.id in 1:length(ratescoex)){
  if(ratescoex[[i.id]]$fit.type[1] == "heterogeneous_both"){
    het.curves <- bind_rows(het.curves,ratescoex[[i.id]])
  }else{
    hom.curves <- bind_rows(hom.curves,ratescoex[[i.id]])
  }
}

if(!is.null(hom.curves)){
  write.csv2(hom.curves,
             file = "results/rate_sensitivity_analysis/S4_modified_rates_accumulation_curves_heterogeneous_time.csv",
             row.names = FALSE)
}

if(!is.null(het.curves)){
  write.csv2(het.curves,
             file = "results/rate_sensitivity_analysis/S4_modified_rates_accumulation_curves_heterogeneous_both.csv",
             row.names = FALSE)
}
# clear cores -------------------------------------------------------------

stopCluster(cl)
beepr::beep(3)
