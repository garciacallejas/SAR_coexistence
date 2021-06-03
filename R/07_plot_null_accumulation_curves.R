
# plot accumulation curves of null model data
# this is Fig 4 of the manuscript

# accumulation curves -----------------------------------------------------

# read data ---------------------------------------------------------------
library(tidyverse)

nullcoex <- read.csv2(file = "./results/07_null_accumulation_curves.csv",header = TRUE,stringsAsFactors = FALSE)
homcoexdf <- read.csv2("./results/03_07_homogeneous_curve.csv",header = T,stringsAsFactors = F)
rich.area <- read.csv2("./data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)

base.year <- 2015:2019

# plot --------------------------------------------------------------------
# include error estimates
std <- function(x) sd(x)/sqrt(length(x))

# richness

rich.mean <- rich.area %>%
  filter(year %in% base.year) %>%
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(richness),
            se.sp = std(richness),
            sd.sp = sd(richness),
            max.sp = max(richness),
            min.sp = min(richness))
# rich.mean$fit.type <- "richness"
rich.mean$area <- rich.mean$n.plots * 8.5^2

# homogeneous parameterization

homaddob <- homcoexdf %>% 
  filter(year %in% base.year) %>% 
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(mean.sp),
            se.sp = std(mean.sp),
            sd.sp = sd(mean.sp),
            maxm.sp = max(max.sp),
            minm.sp = min(min.sp))
names(homaddob)[c(6,7)] <- c("max.sp","min.sp")
homaddob$area <- homaddob$n.plots * 8.5^2

# null models

homaddnull <- nullcoex %>%
  group_by(null_model,replicate,n.plots,year) %>%
  summarise(coexisting.sp = mean(mean.sp),
            se.sp = std(mean.sp),
            sd.sp = mean(mean.sp),
            maxm.sp = max(max.sp),
            minm.sp = min(min.sp))
names(homaddnull)[c(8,9)] <- c("max.sp","min.sp")
homaddnull$area <- homaddnull$n.plots * 8.5^2
 
# include mean null curves
mean.curve <- homaddnull %>% 
  group_by(year,null_model,n.plots) %>% 
  summarise(coex = mean(coexisting.sp))
names(mean.curve)[4] <- "coexisting.sp"
mean.curve$area <- mean.curve$n.plots * 8.5^2

homaddnull1 <- subset(homaddnull,null_model == "KD")
homaddnull2 <- subset(homaddnull,null_model == "RA")
mean.curve1 <- subset(mean.curve,null_model == "KD")
mean.curve2 <- subset(mean.curve,null_model == "RA")

# plot
my.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette.values <- my.palette[c(2,6,4,3,7)]

null.color.1 <- "grey"#"#99DAFF"#"d3ebf8"
null.color.2 <- "grey"#"#66C7FF"
null.color.mean <- palette.values[2]

pd <- 20
ew <- 20
pointsize <- 1.5
barsize <- .5
linesize <- .9

null.coex.plot <- ggplot(homaddob, aes(x = area, y = coexisting.sp))+
  
  # null lines
  geom_line(data = homaddnull1,aes(x = area,y = coexisting.sp, group = replicate),size = .1,color = null.color.1,alpha = .5) +
  geom_line(data = homaddnull2,aes(x = area,y = coexisting.sp, group = replicate),size = .1,color = null.color.2,alpha = .5) +
  
  geom_line(data = mean.curve1,aes(x = area, y = coexisting.sp),color = null.color.mean,linetype = "dashed",size = .6) +
  geom_line(data = mean.curve2,aes(x = area, y = coexisting.sp),color = null.color.mean,linetype = "dotted",size = .6) +
  
  # richness
  geom_line(data = rich.mean,aes(x = area,y = coexisting.sp),color = palette.values[3],
            size = linesize) +
  geom_errorbar(data = rich.mean, aes(ymin = coexisting.sp-se.sp, ymax = coexisting.sp+se.sp),
                color = palette.values[3],
                size = barsize,
                # position = position_dodge(pd),
                width = ew)+
  geom_point(data = rich.mean,aes(x = area,y = coexisting.sp),
             fill = palette.values[3],
             shape = 21,
             size = pointsize) +
  
  # homogeneous curve
  geom_line(color = palette.values[2],
            size = linesize)+
  geom_errorbar(aes(ymin = coexisting.sp-se.sp, ymax = coexisting.sp+se.sp),
                color = palette.values[2],
                size = barsize,
                # position = position_dodge(pd),
                width = ew)+
  geom_point(fill = palette.values[2],
             shape = 21,
             size = pointsize)+

  scale_y_continuous(breaks=seq(0, 21, 2))+
  facet_grid(.~year)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(face = "bold"))+
  theme(panel.grid.minor=element_blank())+
  ylab("average number of species")+
  xlab(bquote('area'~(m^2)))+
  NULL
# null.coex.plot

ggsave(filename = paste("./images/Fig_4.pdf",sep=""),
       plot = null.coex.plot,
       device = cairo_pdf,
       width = 9,height = 3,dpi = 600)
# ggsave(filename = paste("./images/Fig_4.png",sep=""),
#        plot = null.coex.plot,
#        # device = cairo_pdf,
#        width = 9,height = 3,dpi = 600)
