

# Table S1 ----------------------------------------------------------------

# install.packages("mmSAR", repos="http://R-Forge.R-project.org")
library(mmSAR)

sar <- read.csv2("./results/S3_accumulation_curves_data.csv",header = TRUE,stringsAsFactors = FALSE)

years <- sort(unique(sar$year))
types <- sort(unique(sar$fit.type))
areas <- sort(unique(sar$area))

# transform to list of lists

data.list <- list()
i.data <- 1
for(i.year in 1:length(years)){
  for(i.type in 1:length(types)){
    sarname <- paste(types[i.type],"_",years[i.year],sep="")
    sardata <- sar$coexisting.sp[sar$year == years[i.year] &
                                   sar$fit.type == types[i.type]]
    sardf <- data.frame(a = areas,s = sardata)
    
    sarlist <- list(name = sarname,data = sardf)
    
    data.list[[i.data]] <- sarlist
    i.data <- i.data + 1
  }# i.type
}# i.year

# calculate params for each curve
modelList <- c("power","expo")
data("power")
data("expo")

# this is not deterministic, i.e. outcomes may vary slightly,
# but overall, power-law seems a better fit for most curves
i.data <- 1
for(i.year in 1:length(years)){
  for(i.type in 1:length(types)){
    sarfit <- NULL
    try(sarfit <- mmSAR::multiSAR(modelList = modelList,data = data.list[[i.data]], crit = "Info"))
    if(!is.null(sarfit)){
      if(sarfit$DeltaIC[1] == 0){
        mod <- "power"
        dif <- sarfit$DeltaIC[2]
      }else{
        mod <- "expo"
        dif <- sarfit$DeltaIC[1]
      }
      
      # cat("year:",years[i.year],",type:",types[i.type],"- best fit:",mod,"- dif:",dif,"\n")
    }
    i.data <- i.data + 1
  }
}

# create results dataframe
# resdf <- expand.grid(year = years,fit.type = types,model = "power-law",c = 0,z = 0,stringsAsFactors = FALSE)
resdf <- expand.grid(year = years,fit.type = types,
                     parameter = c("c","z"), 
                     estimate = NA, 
                     se = NA,
                     t_value = NA,
                     p_value = NA,
                     lower.bound = NA,
                     upper.bound = NA,
                     stringsAsFactors = FALSE)

# extract parameters from power-law fits
i.data <- 1
for(i.year in 1:length(years)){
  for(i.type in 1:length(types)){
    sarfit <- NULL
    try(sarfit <- mmSAR::rssoptim(power,data.list[[i.data]],verb = FALSE))
    if(!is.null(sarfit)){
      pos.c <- which(resdf$year == years[i.year] &
                       resdf$fit.type == types[i.type] & 
                       resdf$parameter == "c")
      pos.z <- which(resdf$year == years[i.year] &
                       resdf$fit.type == types[i.type] & 
                       resdf$parameter == "z")
      resdf$estimate[pos.c] <- sarfit$par["c"]
      resdf$estimate[pos.z] <- sarfit$par["z"]
      
      resdf$se[pos.c] <- sarfit$sigConf["c",2]
      resdf$se[pos.z] <- sarfit$sigConf["z",2]      
      
      resdf$t_value[pos.c] <- sarfit$sigConf["c",3]
      resdf$t_value[pos.z] <- sarfit$sigConf["z",3]
      
      resdf$p_value[pos.c] <- sarfit$sigConf["c",4]
      resdf$p_value[pos.z] <- sarfit$sigConf["z",4]
      
      resdf$lower.bound[pos.c] <- sarfit$sigConf["c",5]
      resdf$lower.bound[pos.z] <- sarfit$sigConf["z",5]
      
      resdf$upper.bound[pos.c] <- sarfit$sigConf["c",2]
      resdf$upper.bound[pos.z] <- sarfit$sigConf["z",6]
      
      # cat("year:",years[i.year],",type:",types[i.type],"- best fit:",mod,"- dif:",dif,"\n")
    }
    i.data <- i.data + 1
  }
}

# write to disk
# write.csv2(resdf,"./results/S3_curve_power_law_parameterization.csv",row.names = FALSE)

# Table S3 ----------------------------------------------------------------

abund <- read.csv2("data/01_05_abundances.csv",header = TRUE,stringsAsFactors = FALSE)

ab <- abund %>% 
  group_by(year,species) %>%
  summarise(abund = sum(individuals))
names(ab)[2] <- "sp"

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

spr2 <- gather(spr1,key = "sp.type",value = "value",-year,-plot,-sp)
spr <- spr2 %>% filter(value == TRUE) %>%
  group_by(year,sp,sp.type) %>%
  summarise(freq = sum(value))

sprab <- left_join(spr,ab)  

# frequency as transient
sprtr <- spread(sprab,key = sp.type,value = freq,fill = 0) %>%
  mutate(total = dominant+indirect+pairwise+transient) %>%
  mutate(rel.transient = transient/total)

trall <- sprtr %>% group_by(sp) %>% summarise(all.t = sum(transient),
                                              all.total = sum(total), 
                                              rel.t = sum(transient)/sum(total))
avg.tr <- mean(trall$rel.t)

# plot --------------------------------------------------------------------

abtr.plot <- ggplot(sprtr,aes(x = rel.transient,y = log(abund)))+
  geom_point()+
  geom_smooth(method = "lm")+
  NULL
# abtr.plot

# glm analysis ------------------------------------------------------------

library(lme4)
library(DHARMa)

# transient in a given plot/year (0-1) as a function of its abundance
sprab$transient <- ifelse(sprab$sp.type == "transient",1,0)
sprab.model <- subset(sprab,sp.type != "dominant")

# binomial glm
bglm <- glm(transient ~ log(abund),data = sprab.model,family = "binomial")

# model summary and checks
# Table S3 is taken from here
# summary(bglm)
# 
# resb <- DHARMa::simulateResiduals(bglm)
# plot(resb)
# testResiduals(resb)
# plotResiduals(resb)

# ggplot(sprab.model, aes(x=log(abund), y=transient)) + geom_point() + 
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange = TRUE, se=FALSE)


