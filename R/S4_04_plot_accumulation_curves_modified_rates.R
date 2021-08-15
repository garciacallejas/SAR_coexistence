
# Supplementary Table S4 and Supplementary Figure S4

library(tidyverse)
library(colorblindr)

# -------------------------------------------------------------------------
hetcoexdf <- read.csv2("./results/03_heterogeneous_curve.csv",header = T,stringsAsFactors = F)
hetcoexdf$analysis <- "main"

mod.rates <- read.csv2("results/rate_sensitivity_analysis/S4_modified_rates_accumulation_curves_heterogeneous_both.csv")
mod.rates$analysis <- recode(mod.rates$analysis,germ_decrease = "lower germination rate",
                             germ_increase = "higher germination rate",
                             germ_random = "random variation in germination rate",
                             surv_decrease = "lower seed survival",
                             surv_increase = "higher seed survival",
                             surv_random = "random variation in seed survival"
                             )


# -------------------------------------------------------------------------
# random data is replicated 100 times. Move it out, get the averages, move it in
rand.rates <- subset(mod.rates, analysis %in% c("random variation in seed survival",
                                                "random variation in germination rate"))

rand.rates.avg <- rand.rates %>% 
  group_by(n.plots,plots,year,fit.type,analysis) %>%
  summarise(mean.num.sp = mean(num.sp)) %>%
  rename("num.sp" = mean.num.sp)

mod.rates.norand <- subset(mod.rates, !analysis %in% c("random variation in seed survival",
                                                       "random variation in germination rate"))  
mod.rates.avg <- bind_rows(mod.rates.norand,rand.rates.avg)

# -------------------------------------------------------------------------
all.curves <- bind_rows(hetcoexdf,mod.rates.avg)

std <- function(x) sd(x)/sqrt(length(x))
aggregated.data <- all.curves %>% 
  group_by(year,n.plots,analysis,fit.type) %>% 
  summarise(coexisting.sp = mean(num.sp),
            max.sp = max(num.sp),
            min.sp = min(num.sp),
            se.sp = std(num.sp))
aggregated.data$area <- aggregated.data$n.plots * 8.5^2
aggregated.data$n.plots <- as.factor(aggregated.data$n.plots)

# -------------------------------------------------------------------------
my.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette.values <- my.palette[c(2,6,4,3,7,1,5,8)]

pd <- 20
ew <- 50
pointsize <- 1.5
barsize <- .5
linesize <- .9

add.coex.plot <- ggplot(aggregated.data, aes(x = area, y = coexisting.sp, group = analysis))+
  geom_line(aes(group = analysis, color = analysis),
            position = position_dodge(pd),
            size = linesize)+
  
  geom_errorbar(aes(ymin = coexisting.sp-se.sp, ymax = coexisting.sp+se.sp,
                    color = analysis),
                size = barsize,
                position = position_dodge(pd),
                width = ew)+
  
  geom_point(aes(group = analysis, fill = analysis),
             shape = 21,
             position = position_dodge(pd),
             size = pointsize)+
  scale_color_OkabeIto() +
  scale_fill_OkabeIto() +
  # scale_color_manual(values = palette.values,labels = c("heterogeneous\ncoexistence",
  #                                                       "homogeneous\ncoexistence",
  #                                                       "richness")) +
  # scale_fill_manual(values = palette.values,labels = c("heterogeneous\ncoexistence",
  #                                                      "homogeneous\ncoexistence",
  #                                                      "richness")) +
  scale_y_continuous(breaks=seq(0, 21, 2))+
  facet_grid(.~year)+
  theme_bw()+
  # theme(legend.title=element_blank())+
  # theme(legend.justification=c(1.05,-0.05), legend.position=c(1,0))+
  # theme(legend.background = element_rect(fill=alpha('white', 0.4)),
  #       legend.margin = margin(0,0,0,0),
  #       legend.key = element_rect(size = 0.5))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(face = "bold"))+
  theme(panel.grid.minor=element_blank())+
  ylab("average number of species")+
  xlab(bquote('area'~(m^2)))+
  NULL
# add.coex.plot

# ggsave(filename = paste("./images/Fig_S4_v2.pdf",sep=""),
#       plot = add.coex.plot,
#       device = cairo_pdf,
#       width = 8,height = 4,dpi = 600)


# -------------------------------------------------------------------------
# The curves look identical in most cases, except for 2019, higher germ rate
# so, it is interesting to show that there are no underlying errors
# by showing the raw number of feasible combinations in each analysis.
# This shows that, in fact, they are not identical, but slightly different combinations
# can give rise to the same curves
# -------------------------------------------------------------------------
# main data
feas.data <- read.csv2("results/02_observed_structural_metrics.csv")

# damn all the name changes I've been making to these categories!
feas.het <- subset(feas.data,fit == "heterogeneous_space")
names(feas.het)[which(names(feas.het) == "fit")] <- "fit.type"

feas.het$feasibility.domain <- NULL
feas.het$com.pair.differential <- NULL

files <- list.files(path = "./results/rate_sensitivity_analysis/",
                    pattern = "S4_modified_rates_structural_metrics",
                    full.names = TRUE)
rates.files <- files[grep("heterogeneous_both",files)]
rates.data <- rates.files %>%
  map(read.csv2) %>%
  reduce(rbind)

# to be honest, this is all badly coded and sloppy. I don't blame you 
# if you can hardly follow it... it happens to me. I've just spent way
# too much time with this study, my brain barely functions anymore.
# I don't think I will be able to look at a SAR in the next 5? years xD

rates.data <- subset(rates.data, analysis %in% c("germ_increase","germ_decrease",
                                                 "germ_random","surv_random",
                                                 "surv_increase","surv_decrease"))
rates.data$analysis <- recode(rates.data$analysis,germ_decrease = "lower germination rate",
                              germ_increase = "higher germination rate",
                              germ_random = "random variation in germination rate",
                              surv_decrease = "lower seed survival",
                              surv_increase = "higher seed survival",
                              surv_random = "random variation in seed survival"
                              )
rates.data$fit.type <- recode(rates.data$fit.type, heterogeneous_both = "heterogeneous_space")

# -------------------------------------------------------------------------
# random have replicates
# as above, move them out, average, move back in
rates.rand <- subset(rates.data,analysis %in% c("random variation in seed survival",
                                                "random variation in germination rate"))
rates.rand.avg <- rates.rand %>%
  group_by(BEMA,CETE,CHFU,CHMI,HOMA,LEMA,LYTR,MEEL,MESU,PAIN,
           PLCO,POMA,POMO,PUPA,SASO,SCLA,SOAS,SPRU,SUSP,
           year,plot,analysis,fit.type) %>%
  summarise(mean.feasibility = mean(feasibility)) %>%
  rename("feasibility" = mean.feasibility)

rates.no.rand <- subset(rates.data, !analysis %in% c("random variation in seed survival",
                                                "random variation in germination rate"))
rates.no.rand$replicate <- NULL
rates.no.rand$feasibility.domain <- NULL
rates.no.rand$com.pair.differential <- NULL
rates.avg <- bind_rows(rates.rand.avg,rates.no.rand)

rates.w <- rates.avg %>% pivot_wider(names_from = analysis,
                                       values_from = feasibility)

all.data <- left_join(feas.het,rates.w)
names(all.data)[which(names(all.data) == "feasibility")] <- "main"

# -------------------------------------------------------------------------
# find how many feasible combos per year and analysis type
aggr.data <- all.data %>%
  pivot_longer(cols = c(main,`higher germination rate`,
                        `lower germination rate`,
                        `random variation in germination rate`,
                        `higher seed survival`,
                        `lower seed survival`,
                        `random variation in seed survival`
                        ),
               names_to = "analysis",values_to = "feasibility") %>%
  group_by(year,analysis) %>%
  summarise(num.feasible.combos = as.integer(sum(feasibility,na.rm = TRUE))) %>%
  pivot_wider(names_from = year,values_from = num.feasible.combos)

# this is supplementary table 4
# print(xtable::xtable(aggr.data),
#       floating=FALSE,
#       latex.environments=NULL,
#       booktabs=TRUE,
#       include.rownames = FALSE)




