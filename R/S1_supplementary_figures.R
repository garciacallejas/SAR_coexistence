
# -------------------------------------------------------------------------
# Figure S1

library(tidyverse)
library(colorblindr)

het.alphas <- read.csv2(file = "data/01_alpha_heterogeneous.csv",header = TRUE,stringsAsFactors = FALSE)
hom.alphas <- read.csv2(file = "data/01_05_alpha_homogeneous.csv",header = TRUE,stringsAsFactors = FALSE)

years <- sort(unique(het.alphas$year))
plots <- sort(unique(het.alphas$plot))
focals <- sort(unique(het.alphas$focal))

hom.diag.dom <- expand.grid(year = years,plot = plots, focal = focals, model = "homogeneous", intra = NA_real_, mean.inter = NA_real_, sum.inter = NA_real_,stringsAsFactors = FALSE)
het.diag.dom <- expand.grid(year = years,plot = plots, focal = focals, model = "heterogeneous",intra = NA_real_, mean.inter = NA_real_, sum.inter = NA_real_,stringsAsFactors = FALSE)

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.focal in 1:length(focals)){
      hom.row <- subset(hom.alphas, year == years[i.year] & plot == plots[i.plot] & focal == focals[i.focal])
      het.row <- subset(het.alphas, year == years[i.year] & plot == plots[i.plot] & focal == focals[i.focal])
      
      hom.pos <- which(hom.diag.dom$year == years[i.year] & hom.diag.dom$plot == plots[i.plot] & hom.diag.dom$focal == focals[i.focal])
      het.pos <- which(het.diag.dom$year == years[i.year] & het.diag.dom$plot == plots[i.plot] & het.diag.dom$focal == focals[i.focal])
      
      hom.diag.dom$intra[hom.pos] <- hom.row$magnitude[hom.row$neighbour == focals[i.focal]]
      hom.diag.dom$mean.inter[hom.pos] <- mean(hom.row$magnitude[which(!is.na(hom.row$magnitude) & hom.row$neighbour != focals[i.focal])])
      hom.diag.dom$sum.inter[hom.pos] <- sum(hom.row$magnitude[which(!is.na(hom.row$magnitude) & hom.row$neighbour != focals[i.focal])])
      
      het.diag.dom$intra[het.pos] <- het.row$magnitude[het.row$neighbour == focals[i.focal]]
      het.diag.dom$mean.inter[het.pos] <- mean(het.row$magnitude[which(!is.na(het.row$magnitude) & het.row$neighbour != focals[i.focal])])
      het.diag.dom$sum.inter[het.pos] <- sum(het.row$magnitude[which(!is.na(het.row$magnitude) & het.row$neighbour != focals[i.focal])])
      
    }# i.focal
  }# i.plot
}# i.year

interactions.observed <- bind_rows(hom.diag.dom,het.diag.dom)
# write.csv2(diag.all.models,file = "./results/S2_interactions_summary.csv",row.names = FALSE)
# interactions.observed
interactions.observed$year <- as.factor(interactions.observed$year)
interactions.observed$intra.inter.ratio <- interactions.observed$intra/interactions.observed$mean.inter

# n
# interactions.observed %>%
#   group_by(model,year) %>%
#   filter(!is.na(intra)) %>%
#   summarise(num.intra = n())

intra.plot <- ggplot(interactions.observed, aes(x = year, y = intra, fill = model))+
  geom_boxplot(alpha = 0.60,outlier.size = .7,outlier.colour = "darkgrey") +
  scale_fill_OkabeIto(order = c(1,5)) +
  guides(fill = FALSE) +
  theme_bw()+
  ggtitle("a)")+
  labs(y = "Intraspecific\ninteraction strength")+
  NULL
yintra = boxplot.stats(interactions.observed$intra)$stats[c(1, 5)]
intra.scale.plot <- intra.plot + coord_cartesian(ylim = yintra*3.5)
intra.scale.plot

# n
# interactions.observed %>%
#   group_by(model,year) %>%
#   filter(!is.na(mean.inter)) %>%
#   summarise(mean.inter = n())

inter.plot <- ggplot(interactions.observed, aes(x = year, y = mean.inter, fill = model))+
  geom_boxplot(alpha = 0.60,outlier.size = .7,outlier.colour = "darkgrey") +
  scale_fill_OkabeIto(order = c(1,5)) +
  guides(fill = FALSE) +
  theme_bw()+
  ggtitle("b)")+
  labs(y = "Average interspecific\ninteraction strength")+
  NULL
yinter = boxplot.stats(interactions.observed$mean.inter)$stats[c(1, 5)]
inter.scale.plot <- inter.plot + coord_cartesian(ylim = yinter*2)
inter.scale.plot

# n
# interactions.observed %>%
#   group_by(model,year) %>%
#   filter(!is.na(intra.inter.ratio)) %>%
#   summarise(num.ratio = n())

ratio.plot <- ggplot(interactions.observed, aes(x = year, y = intra.inter.ratio, fill = model))+
  geom_boxplot(alpha = 0.60,outlier.size = .7,outlier.colour = "darkgrey") +
  scale_fill_OkabeIto(order = c(1,5)) +
  guides(fill = FALSE) +
  theme_bw()+
  ggtitle("c)")+
  labs(y = "Average ratio\nintra/interspecific strength")+
  NULL
yratio = boxplot.stats(interactions.observed$intra.inter.ratio)$stats[c(1, 5)]
ratio.scale.plot <- ratio.plot + coord_cartesian(ylim = yratio*1.5)
ratio.scale.plot

library(patchwork)
is.plot <- intra.scale.plot + inter.scale.plot + ratio.scale.plot
# ggsave(filename = "./images/Fig_S1.pdf",is.plot,device = cairo_pdf,width = 11,height = 4,dpi = 600)

# -------------------------------------------------------------------------
# Figure S2

intra.inter.pair <- list()
param <- c("heterogeneous","homogeneous")
all.sp <- focals

for(i.param in 1:length(param)){
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      if(param[i.param] == "heterogeneous"){
        my.data <- subset(het.alphas,
                          year == years[i.year] & plot == plots[i.plot])
      }else{
        my.data <- subset(hom.alphas,
                          year == years[i.year] & plot == plots[i.plot])
      }
      
      ii <- expand.grid(year = years[i.year],
                        plot = plots[i.plot],
                        parameterization = param[i.param],
                        sp1 = all.sp,
                        sp2 = all.sp,
                        intra.sp1 = NA_real_,
                        intra.sp2 = NA_real_,
                        inter = NA_real_)
      # no need for sp1 = sp2
      ii <- subset(ii, sp1 != sp2)
      
      # all intraspecific coefs for this year-plot
      all.intras <- data.frame(sp = all.sp,intra = NA_real_)
      for(i.sp in 1:length(all.sp)){
        all.intras$intra[all.intras$sp == all.sp[i.sp]] <- 
          my.data$magnitude[my.data$focal == my.data$neighbour & 
                              my.data$focal == all.sp[i.sp]]
      }# for i.sp
      
      # fill up all coefs for each pair
      for(i.sp in 1:length(all.sp)){
        for(j.sp in 1:length(all.sp)){
          if(i.sp != j.sp){
            ii$intra.sp1[ii$sp1 == all.sp[i.sp] &
                           ii$sp2 == all.sp[j.sp]] <- 
              all.intras$intra[all.intras$sp == all.sp[i.sp]]
            
            ii$intra.sp2[ii$sp1 == all.sp[i.sp] &
                           ii$sp2 == all.sp[j.sp]] <- 
              all.intras$intra[all.intras$sp == all.sp[j.sp]]
            
            ii$inter[ii$sp1 == all.sp[i.sp] &
                       ii$sp2 == all.sp[j.sp]] <- 
              my.data$magnitude[my.data$focal == all.sp[i.sp] & 
                                  my.data$neighbour == all.sp[j.sp]]
            
          }
        }# for j.sp
      }# for i.sp
      
      intra.inter.pair[[length(intra.inter.pair) +1]] <- ii
      
    }# for i.plot
  }# for i.year
}# for i.param

iip <- bind_rows(intra.inter.pair)
iip <- bind_rows(intra.inter.pair)

# distribution of intra-inter

iip$sp2.minus.inter <- iip$intra.sp2 - iip$inter

summary.ii <- iip %>% 
  filter(!is.na(intra.sp1) & !is.na(intra.sp2) & !is.na(inter)) %>%
  group_by(year,plot,parameterization) %>%
  mutate(num = n()) %>%
  group_by(year,plot,parameterization) %>%
  summarise(higher.intra.sp1 = sum(intra.sp1 > inter)/num,
            higher.intra.sp2 = sum(intra.sp2 > inter)/num,
            higher.intra.both = sum(intra.sp1 > inter & intra.sp2 > inter)/num)

summary.ii$year <- as.factor(summary.ii$year)

summary.ii.het <- summary.ii[summary.ii$parameterization == "heterogeneous",]

# n
summary.ii.het %>%
  group_by(year) %>%
  summarise(num = n()) 


s2 <- ggplot(summary.ii.het, aes(x = year, y = higher.intra.both)) + 
  # geom_point(fill = "grey", shape = 21, size = .9, position = position_dodge(.2)) + 
  geom_boxplot(fill = "lightgrey",outlier.colour = "grey",outlier.size = .8) +
  # geom_hline(yintercept = .5, linetype = "dashed", color = "darkgrey") +
  # facet_grid(parameterization~.)+
  theme_bw()+
  # theme(panel.grid.minor=element_blank())+
  # theme(strip.background = element_blank())+
  # theme(strip.text.y = element_text(face = "bold"))+
  ylab("frequency of higher intraspecific strength") +
  NULL

# ggsave("results/images/Fig_S2.pdf",
#        plot = s2,
#        device = cairo_pdf,
#        width = 8,height = 4,dpi = 600)




