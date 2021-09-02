
library(tidyverse)

# initial data ------------------------------------------------------------

local.stab.het <- read.csv2(file = "results/local_stability_analysis/S5_local_stability_heterogeneous.csv")
local.stab.hom <- read.csv2(file = "results/local_stability_analysis/S5_local_stability_homogeneous.csv")

# -------------------------------------------------------------------------
# join both datasets

local.stab.het$fit <- "heterogeneous_space"
local.stab.hom$fit <- "homogeneous_space"

local.stab <- bind_rows(local.stab.het,local.stab.hom)
local.stab <- subset(local.stab,!is.na(max_eigenvalue))
sum(local.stab$max_eigenvalue <= 1)

# -------------------------------------------------------------------------
# it does not make sense to analyze any more: there are only 2 locally stable combos
# out of ~500k

# Further studies may look into this, mathematically and ecologically.

