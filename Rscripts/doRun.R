

#library(devtools)
#install_github("FLa4a", "colinpmillar")

# select a setup option

#source("simpleSetup.R")
source("BETSetup.R")

nsim <- 1000
biomass <- FALSE

# run simulations
source("simulation.R")

head(sims)

#print(bwplot(sdiff ~ parname | type, data = sims, scales = list(relation = "free")))

# plot catchability estimates
qest <- array(Xq %*% matrix(subset(sims, type == "qMod") $ est, length(bq)), dim = c(nages, nyears, nsim))[,1,]
qsim <- matrix(Xq %*% bq, nages, nyears)[,1]

plot(ages, qsim, type = "l", ylim = range(qest, qsim))
apply(qest, 2, lines, x = ages, col = grey(0.8))



