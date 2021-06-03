
# plot the partition of coexisting species into different categories,
# namely species coexisting in pairs, networks, or species that do not coexist
# and are either dominant or transient.

# this is done for the heterogeneous space parameters
# Figure 3 of the manuscript

# heterogeneous curve partition -------------------------------------------
library(tidyverse)
sproles <- read.csv2(file = "./results/04_observed_species_roles.csv",header = T,stringsAsFactors = F)

years <- sort(unique(sproles$year))
all.plots <- sort(unique(sproles$plot))
n.plots <- length(all.plots)
all.sp <- sort(unique(sproles$sp))

# first, if a species enters from pairs, it does not enter from indirect
spr1 <- sproles
spr1$indirect2 <- spr1$indirect
spr1$indirect <- ifelse(spr1$pairwise == TRUE,FALSE,spr1$indirect2)
spr1$indirect2 <- NULL

# this should be equal to the richness values in rich.area
# spnum <- spr1 %>% group_by(year,plot) %>% summarise(num = sum(pairwise == TRUE | 
#                                                                 indirect == TRUE | 
#                                                                 dominant == TRUE | 
#                                                                 transient == TRUE))

spr2 <- gather(spr1,key = "sp.type",value = "value",-year,-plot,-sp)
spr3 <- spread(spr2,key = sp,value = value)

categories <- sort(c(unique(spr3$sp.type),"absent"))

# obtain, for each area, the
# percentage of sp that appear on a single category
# e.g. if sp A appears as pairwise in plot 1 and
# transient in plot 2, in area 1+2, it is pairwise
catacc <- NULL
for(i.year in 1:length(years)){
  for(i.num in 1:n.plots){
    my.combinations <- combn(1:9,i.num)
    for(i.comb in 1:ncol(my.combinations)){
      
      # to keep track of the category of each sp
      sp.cat <- character(length(all.sp))
      names(sp.cat) <- all.sp
      
      my.plots <- as.character(my.combinations[,i.comb])
      # concatenate into a single string
      my.plots.char <- paste(my.plots,collapse = "_")
      
      my.data <- subset(spr2,year == years[i.year] & 
                          plot %in% my.plots)
      
      # update sp.cat
      # tidy would be nice, 
      # but maybe too difficult
      for(i.sp in 1:length(sp.cat)){
        if(any(my.data$value[my.data$sp == names(sp.cat)[i.sp] &
                             my.data$sp.type == "pairwise"])){
          sp.cat[i.sp] <- "pairwise"
        }else if(any(my.data$value[my.data$sp == names(sp.cat)[i.sp] &
                                   my.data$sp.type == "indirect"])){
          sp.cat[i.sp] <- "indirect"
        }else if(any(my.data$value[my.data$sp == names(sp.cat)[i.sp] &
                                   my.data$sp.type == "dominant"])){
          sp.cat[i.sp] <- "dominant"
        }else if(any(my.data$value[my.data$sp == names(sp.cat)[i.sp] &
                                   my.data$sp.type == "transient"])){
          sp.cat[i.sp] <- "transient"
        }else{
          sp.cat[i.sp] <- "absent"
        }
      }
      
      freq <- table(sp.cat)
      
      comb.data <- expand.grid(year = years[i.year],
                               n.plots = i.num,
                               plots = my.plots.char,
                               category = categories,
                               num.sp = 0)
      
      for(i.cat in 1:length(categories)){
        comb.data$num.sp[comb.data$category == categories[i.cat]] <- 
          ifelse(is.na(freq[categories[i.cat]]),0,freq[categories[i.cat]])
      }
      
      catacc <- rbind(catacc,comb.data)
      
    }# for i.comb
  }# for i.num
}# for i.year

catmean <- catacc %>% 
  group_by(category,year,n.plots) %>% 
  summarise(num = mean(num.sp),
            median.num = median(num.sp),
            max.sp = max(num.sp),
            min.sp = min(num.sp))

# percentage of transient species
transient.ratio <- catmean %>%
  group_by(year,n.plots,category) %>%
  summarise(ratio = num/19)
avg.transient.ratio <- transient.ratio %>% 
  filter(category == "transient") %>% 
  # group_by(category) %>%
  summarise(avg.transient.ratio = mean(ratio))

# species roles plot ------------------------------------------------------

catmean$category <- factor(catmean$category,levels = c("transient","dominant","indirect","pairwise","absent"))

# rename to comply with the terminology of the main text
catmean$category <- dplyr::recode(catmean$category, indirect = "multispecies coexistence", pairwise = "pairwise coexistence")
catmean$area <- catmean$n.plots * 8.5^2

# custom palette
library(wesanderson)
custom_palette <- wes_palette("Zissou1", 4, type = "discrete")

catmean2 <- subset(catmean, category != "absent")

catbar <- ggplot(catmean2,aes(x = area, y = num))+
  geom_col(aes(fill = category))+
  facet_grid(.~year)+
  # scale_fill_manual(name="category",values = palette.values) +
  scale_fill_manual(name = "category", values = custom_palette)+
                       # direction = -1,
                       # option = "E") +
  scale_y_continuous(breaks=c(0,seq(1, 20, 2)))+
  
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(face = "bold"))+
  # theme(panel.grid.minor=element_blank())+
  ylab("average number of species")+
  xlab(bquote('area'~(m^2)))+
  NULL
catbar

ggsave(filename = paste("./images/Fig_3.pdf",sep=""),
       plot = catbar,
       device = cairo_pdf,
       width = 9,height = 3,dpi = 600)

# ggsave(filename = paste("./images/Fig_3.png",sep=""),
#        plot = catbar,
#        # device = cairo_pdf,
#        width = 9,height = 3,dpi = 600)
