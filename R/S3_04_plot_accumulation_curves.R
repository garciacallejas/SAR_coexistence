# plot accumulation curves from posterior draws of alpha and lambda coefficients
# this is Fig S3 of the manuscript

# read data ---------------------------------------------------------------
library(tidyverse)

hetrep <- read.csv2("./results/posterior_draws/S3_posterior_draws_accumulation_curves_heterogeneous_both.csv")
homrep <- read.csv2("./results/posterior_draws/S3_posterior_draws_accumulation_curves_heterogeneous_time.csv")

homcoexdf <- read.csv2("./results/03_07_homogeneous_curve.csv",header = T,stringsAsFactors = F)
hetcoexdf <- read.csv2("./results/03_heterogeneous_curve.csv",header = T,stringsAsFactors = F)
rich.area <- read.csv2("./data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)

base.year <- 2015:2019

# plot --------------------------------------------------------------------

# richness

rich.mean <- rich.area %>%
  filter(year %in% base.year) %>%
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(richness),
            max.sp = max(richness),
            min.sp = min(richness))
# rich.mean$fit.type <- "richness"
rich.mean$area <- rich.mean$n.plots * 8.5^2

# homogeneous parameterization

homaddob <- homcoexdf %>% 
  filter(year %in% base.year) %>% 
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(mean.sp),
            maxm.sp = max(max.sp),
            minm.sp = min(min.sp))
names(homaddob)[c(4,5)] <- c("max.sp","min.sp")
homaddob$area <- homaddob$n.plots * 8.5^2
homaddob$fit.type <- "homogeneous"

# heterogeneous parameterization

hetaddob <- hetcoexdf %>% 
  filter(year %in% base.year) %>%
  group_by(year,n.plots) %>% 
  summarise(coexisting.sp = mean(num.sp),
            max.sp = max(num.sp),
            min.sp = min(num.sp))
hetaddob$area <- hetaddob$n.plots * 8.5^2
hetaddob$fit.type <- "heterogeneous"

addob <- rbind(homaddob,hetaddob)

# null models

# heterogeneous
hetaddnull <- hetrep %>%
  filter(year %in% base.year) %>%
  group_by(fit.type,replicate,n.plots,year) %>%
  summarise(coexisting.sp = mean(num.sp),
            max.sp = max(num.sp),
            min.sp = min(num.sp))
# names(addnull)[c(6,7)] <- c("max.sp","min.sp")
hetaddnull$area <- hetaddnull$n.plots * 8.5^2
hetaddnull$fit.type <- "heterogeneous"

homaddnull <- homrep %>%
  filter(year %in% base.year) %>% 
  group_by(replicate,year,n.plots) %>%
  summarise(coexisting.sp = mean(mean.sp),
            maxm.sp = max(max.sp),
            minm.sp = min(min.sp))
names(homaddnull)[c(5,6)] <- c("max.sp","min.sp")
homaddnull$area <- homaddnull$n.plots * 8.5^2
homaddnull$fit.type <- "homogeneous"

addnull <- bind_rows(hetaddnull,homaddnull)

# plot
my.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette.values <- my.palette[c(2,6,4,3,7)]

null.color.1 <- "grey"#"#99DAFF"#"d3ebf8"
null.color.2 <- "grey"#"#66C7FF"

n1 <- ggplot(addob, aes(x = area, y = coexisting.sp))+
  
  # null lines
  geom_line(data = addnull,aes(x = area,y = coexisting.sp, group = replicate),
            size = .1,
            color = null.color.1,
            alpha = .8) +
  
  # mean curves
  geom_line(aes(color = fit.type),size = .7)+
  geom_point(aes(fill = fit.type),shape = 21,size = 2)+
  scale_color_manual(name = "",values = palette.values,labels = c("heterogeneous\ncoexistence",
                                                        "homogeneous\ncoexistence",
                                                        "richness")) +
  scale_fill_manual(name = "",values = palette.values,labels = c("heterogeneous\ncoexistence",
                                                       "homogeneous\ncoexistence",
                                                       "richness")) +
  scale_y_continuous(breaks=seq(0, 21, 2))+
  facet_grid(fit.type~year)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(face = "bold"))+
  theme(panel.grid.minor=element_blank())+
  theme(legend.position = "none") +
  ylab("average number of species")+
  xlab(bquote('area'~(m^2)))+
  NULL
# n1

# ggsave(filename = paste("./images/Fig_S3.pdf",sep=""),
#        plot = n1,
#        device = cairo_pdf,
#        width = 7,height = 3,dpi = 600)
# ggsave(filename = paste("./images/Fig_S3.png",sep=""),
#        plot = n1,
#        # device = cairo_pdf,
#        width = 7,height = 3,dpi = 600)
