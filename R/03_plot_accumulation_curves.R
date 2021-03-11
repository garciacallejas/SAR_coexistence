
# plot coexistence-area curves alongside
# species-area curves. 
# This is Fig. 2 of the manuscript.

# read data ---------------------------------------------------------------
library(tidyverse)

hom <- read.csv2("./results/03_07_homogeneous_curve.csv",header = T,stringsAsFactors = F)
het <- read.csv2("./results/03_heterogeneous_curve.csv",header = T,stringsAsFactors = F)
rich.area <- read.csv2("./data/01_03_05_07_cummulative_richness.csv",header = T,stringsAsFactors = F)

# plot curves --------------------------------------------------

# average by number of plots

rich.mean <- rich.area %>% 
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(richness),
            max.sp = max(richness),
            min.sp = min(richness))
rich.mean$fit.type <- "richness"

homadd <- hom %>% 
  group_by(year,n.plots) %>%
  summarise(coexisting.sp = mean(mean.sp),
            maxm.sp = max(max.sp),
            minm.sp = min(min.sp))
names(homadd)[c(4,5)] <- c("max.sp","min.sp")
homadd$fit.type <- "homogeneous"

hetadd <- het %>% 
  group_by(fit.type,year,n.plots) %>% 
  summarise(coexisting.sp = mean(num.sp),
            max.sp = max(num.sp),
            min.sp = min(num.sp))

# full dataset
addrcoex <- bind_rows(rich.mean,homadd,hetadd)
addrcoex$n.plots <- as.factor(addrcoex$n.plots)
addrcoex$area <- as.numeric(as.character(addrcoex$n.plots)) * 8.5^2

# plot
my.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette.values <- my.palette[c(2,6,4,3,7)]

add.coex.plot <- ggplot(addrcoex, aes(x = area, y = coexisting.sp))+
  geom_line(aes(group = fit.type, color = fit.type),size = 1.2)+
  geom_point(aes(group = fit.type, fill = fit.type),shape = 21,position = position_dodge(.2),size = 3)+
  scale_color_manual(values = palette.values,labels = c("heterogeneous\ncoexistence",
                                                                       "homogeneous\ncoexistence",
                                                                       "richness")) +
  scale_fill_manual(values = palette.values,labels = c("heterogeneous\ncoexistence",
                                                                      "homogeneous\ncoexistence",
                                                                      "richness")) +
  scale_y_continuous(breaks=seq(0, 21, 2))+
  facet_grid(.~year)+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.justification=c(1.05,-0.05), legend.position=c(1,0))+
  theme(legend.background = element_rect(fill=alpha('white', 0.4)),
        legend.margin = margin(0,0,0,0),
        legend.key = element_rect(size = 0.5))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(face = "bold"))+
  theme(panel.grid.minor=element_blank())+
  ylab("species")+
  xlab(bquote('area'~(m^2)))+
  NULL
add.coex.plot

ggsave(filename = paste("./images/Fig_2.pdf",sep=""),plot = add.coex.plot,device = cairo_pdf,width = 9,height = 3,dpi = 600)
