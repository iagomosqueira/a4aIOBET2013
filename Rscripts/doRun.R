

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

load("01.22.Wednesday.RData")

require(FLa4a)

par(mfrow = c(2,2))

# plot catchability estimates
qest <- exp(array(Xq %*% matrix(subset(sims, type == "qMod") $ est, length(bq)), dim = c(nages, nyears, nsim))[,1,])
qsim <- exp(matrix(Xq %*% bq, nages, nyears)[,1])

plot(ages, qsim, type = "l", ylim = range(qest, qsim))
apply(qest, 2, lines, x = ages, col = grey(0.8))
lines(ages, qsim)


# plot recruitment estimates
rest <- exp(matrix(Xr %*% matrix(subset(sims, type == "rMod") $ est, length(br)), nyears, nsim))
rsim <- exp(c(Xr %*% br))

plot(years, rsim, type = "l", ylim = range(rest, rsim))
apply(rest, 2, lines, x = years, col = grey(0.8))
lines(years, rsim)

# plot fbar estimates
fest <- array(Xf %*% matrix(subset(sims, type == "fMod") $ est, length(bf)), dim = c(nages, nyears, nsim))
fsim <- matrix(Xf %*% bf, nages, nyears)

fbarsim <- colMeans(exp(fsim[fbar[1]:fbar[2],]))
fbarest <- apply(exp(fest[fbar[1]:fbar[2],,]), 2:3, mean)

plot(years, fbarsim, type = "l", ylim = range(fbarest, fbarsim))
apply(fbarest, 2, lines, x = years, col = grey(0.8))
lines(years, fbarsim)

# plot fage estimates for last year
fest <- exp(array(Xf %*% matrix(subset(sims, type == "fMod") $ est, length(bf)), dim = c(nages, nyears, nsim))[,nyears,])
fsim <- exp(matrix(Xf %*% bf, nages, nyears)[,nyears])

plot(ages, fsim, type = "l", ylim = range(fest, fsim))
apply(fest, 2, lines, x = ages, col = grey(0.8))
lines(ages, fsim)


