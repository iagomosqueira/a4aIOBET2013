

#library(devtools)
#install_github("FLa4a", "colinpmillar")

# select a setup option

source("simpleSetup.R")
#source("BETSetup.R")

nsim <- 100
biomass <- TRUE

# run simulations
source("simulation.R")

# summarise sims
sims <- do.call(rbind, simout)

print(bwplot(sdiff ~ parname | type, data = sims, scales = list(relation = "free")))


