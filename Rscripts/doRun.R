

#library(devtools)
#install_github("FLa4a", "colinpmillar")

library(FLa4a)
source("setupSimulation.R")
source("simulation.R")



# select a setup option

load('../data/bet_iotc_2011.RData')
range(bet)[c("minfbar","maxfbar")] <- c(4,9)
units(harvest(bet)) <- "f"
load('../data/betSS3.RData')
bet <- window(bet, 1952, 2011)
harvest[harvest == 0] <- min(harvest[harvest > 0])
harvest(bet) <- c(harvest)
stock.n(bet) <- c(stock.n[,-ncol(stock.n)])

biomass <- TRUE

fmodel <- ~ s(age, k = 4) + s(year, k = 20)
srmodel <- ~ factor(year)
n1model <- ~ factor(age)
qmodel <- list(~ 1)

sim <- setupStock(fmodel, srmodel, n1model, bet)
sim <- c(sim, setupIndices(sim $ stock, biomass, qtype = "flat"))

args <- c(sim, list(nsim = 10))
# overwrite qmodel to fit with something different than was used to sim
args $ qmodel <- qmodel

# run simulations
sims <- do.call(doSimulation, args)

head(sims)

#print(bwplot(sdiff ~ parname | type, data = sims, scales = list(relation = "free")))


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


