
# sensitivity analysis for germination and survival rates
# generate modified rates and store them in a csv file

library(tidyverse)

sp.rates <- read.csv2(file = "data/01_05_plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)

# tidy the table
obs.rates <- sp.rates %>% select(-species.name) %>%
  pivot_longer(!species.code,names_to = "rate",values_to = "value") %>%
  mutate(analysis = "main") %>%
  select(species.code,analysis,rate,value)

# mean and sd of the sampling distribution
rate.dist <- obs.rates %>% group_by(rate) %>%
  summarise(mean.rate = mean(value),sd.rate = sd(value))

# generate new samples for each rate (survival, germination, both) and three scenarios:
# random variation, increase, decrease

# -------------------------------------------------------------------------
# increase in seed.survival
si.rates <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "surv_increase",
            rate = "seed.survival",
            mod.value = value + rate.dist$sd.rate[rate.dist$rate == "seed.survival"])
names(si.rates)[4] <- "value"

# just in case
si.rates$value[si.rates$value > 1] <- 1
# attach the other, non-modified rate
si.rates <- bind_rows(si.rates,obs.rates[obs.rates$rate == "germination.rate",])
si.rates$analysis <- "surv_increase"
si.rates$replicate <- 1

# -------------------------------------------------------------------------
# decrease in seed.survival
sd.rates <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "surv_decrease",
            rate = "seed.survival",
            mod.value = value - rate.dist$sd.rate[rate.dist$rate == "seed.survival"])
names(sd.rates)[4] <- "value"

# just in case
sd.rates$value[sd.rates$value < 0] <- 0.01
# attach the other, non-modified rate
sd.rates <- bind_rows(sd.rates,obs.rates[obs.rates$rate == "germination.rate",])
sd.rates$analysis <- "surv_decrease"
sd.rates$replicate <- 1

# -------------------------------------------------------------------------
# random variation in seed.survival
sr.rates <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "surv_random",
            rate = "seed.survival",
            mod.value = value + runif(1,-rate.dist$sd.rate[rate.dist$rate == "seed.survival"],
                                      rate.dist$sd.rate[rate.dist$rate == "seed.survival"]),
            replicate = 1:100)
names(sr.rates)[4] <- "value"

# just in case
sr.rates$value[sr.rates$value > 1] <- 1
sr.rates$value[sr.rates$value < 0] <- 0.01
# attach the other, non-modified rate
# in this case, replicated 100 times
rep.germ <- obs.rates[obs.rates$rate == "germination.rate",] %>% 
  slice(rep(row_number(), 100)) %>%
  mutate(replicate = sort(rep(1:100,nrow(sp.rates))))
sr.rates <- bind_rows(sr.rates,rep.germ)
sr.rates$analysis <- "surv_random"

# -------------------------------------------------------------------------
# increase in germination.rate
gi.rates <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "germ_increase",
            rate = "germination.rate",
            mod.value = value + rate.dist$sd.rate[rate.dist$rate == "germination.rate"])
names(gi.rates)[4] <- "value"

# just in case
gi.rates$value[gi.rates$value > 1] <- 1
# attach the other, non-modified rate
gi.rates <- bind_rows(gi.rates,obs.rates[obs.rates$rate == "seed.survival",])
gi.rates$analysis <- "germ_increase"
gi.rates$replicate <- 1

# -------------------------------------------------------------------------
# decrease in germination.rate
gd.rates <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "germ_decrease",
            rate = "germination.rate",
            mod.value = value - rate.dist$sd.rate[rate.dist$rate == "germination.rate"])
names(gd.rates)[4] <- "value"

# just in case
gd.rates$value[gd.rates$value < 0] <- 0.01
# attach the other, non-modified rate
gd.rates <- bind_rows(gd.rates,obs.rates[obs.rates$rate == "seed.survival",])
gd.rates$analysis <- "germ_decrease"
gd.rates$replicate <- 1

# -------------------------------------------------------------------------
# random variation in germination.rate
gr.rates <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "germ_random",
            rate = "germination.rate",
            mod.value = value + runif(100,-rate.dist$sd.rate[rate.dist$rate == "germination.rate"],
                                      rate.dist$sd.rate[rate.dist$rate == "germination.rate"]),
            replicate = 1:100)
names(gr.rates)[4] <- "value"

# just in case
gr.rates$value[gr.rates$value > 1] <- 1
gr.rates$value[gr.rates$value < 0] <- 0.01
# attach the other, non-modified rate
# in this case, replicated 100 times
rep.surv <- obs.rates[obs.rates$rate == "seed.survival",] %>% 
  slice(rep(row_number(), 100)) %>%
  mutate(replicate = sort(rep(1:100,nrow(sp.rates))))
gr.rates <- bind_rows(gr.rates,rep.surv)
gr.rates$analysis <- "germ_random"

# -------------------------------------------------------------------------
# increase in both
bi.rates.g <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_increase",
            rate = "germination.rate",
            mod.value = value + rate.dist$sd.rate[rate.dist$rate == "germination.rate"])
names(bi.rates.g)[4] <- "value"
# just in case
bi.rates.g$value[bi.rates.g$value > 1] <- 1

bi.rates.s <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_increase",
            rate = "seed.survival",
            mod.value = value + rate.dist$sd.rate[rate.dist$rate == "seed.survival"])
names(bi.rates.s)[4] <- "value"
# just in case
bi.rates.s$value[bi.rates.s$value > 1] <- 1

bi.rates <- bind_rows(bi.rates.g,bi.rates.s)
bi.rates$replicate <- 1

# -------------------------------------------------------------------------
# decrease in both
bd.rates.g <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_decrease",
            rate = "germination.rate",
            mod.value = value - rate.dist$sd.rate[rate.dist$rate == "germination.rate"])
names(bd.rates.g)[4] <- "value"
# just in case
bd.rates.g$value[bd.rates.g$value < 0] <- 0.01

bd.rates.s <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_decrease",
            rate = "seed.survival",
            mod.value = value - rate.dist$sd.rate[rate.dist$rate == "seed.survival"])
names(bd.rates.s)[4] <- "value"
# just in case
bd.rates.s$value[bd.rates.s$value < 0] <- 0.01

bd.rates <- bind_rows(bd.rates.g,bd.rates.s)
bd.rates$replicate <- 1

# -------------------------------------------------------------------------
# random variation in both
br.rates.g <- obs.rates %>% 
  filter(rate == "germination.rate") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_random",
            rate = "germination.rate",
            mod.value = value + runif(100,-rate.dist$sd.rate[rate.dist$rate == "germination.rate"],
                                      rate.dist$sd.rate[rate.dist$rate == "germination.rate"]))
names(br.rates.g)[4] <- "value"

# just in case
br.rates.g$value[br.rates.g$value > 1] <- 1
br.rates.g$value[br.rates.g$value < 0] <- 0.01

br.rates.s <- obs.rates %>% 
  filter(rate == "seed.survival") %>%
  group_by(species.code) %>%
  summarise(analysis = "both_random",
            rate = "seed.survival",
            mod.value = value + runif(100,-rate.dist$sd.rate[rate.dist$rate == "seed.survival"],
                                      rate.dist$sd.rate[rate.dist$rate == "seed.survival"]))
names(br.rates.s)[4] <- "value"

# just in case
br.rates.s$value[br.rates.s$value > 1] <- 1
br.rates.s$value[br.rates.s$value < 0] <- 0.01

br.rates <- bind_rows(br.rates.g,br.rates.s)

# -------------------------------------------------------------------------
# put together all replicates
sim.rates <- bind_rows(#sim.rates,
                       si.rates,
                       sd.rates,
                       sr.rates,
                       gi.rates,
                       gd.rates,
                       gr.rates,
                       bi.rates,
                       bd.rates,
                       br.rates) %>%
  arrange(analysis,replicate,species.code,rate)

write.csv2(sim.rates,
           file = "data/rate_sensitivity_analysis/S4_01_modified_vital_rates.csv",
           row.names = FALSE)




