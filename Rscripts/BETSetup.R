
require(FLCore)
load('../data/bet_iotc_2011.RData')

# CHANGE Fully selected ages
range(bet)[c("minfbar","maxfbar")] <- c(4,9)
units(harvest(bet)) <- "f"

dms <- dims(bet)

years <- dms$minyear:dms$maxyear 
ages <- dms$min:dms$max

nages <- length(ages)
nyears <- length(years)

# for now make up some F values and N values
harvest(bet) <- c(outer(c(scale(c(.1,.1,.2,.2,.3,.35,.4, 0.45, .3, .2), center=FALSE)), 
                      c(rep(0.01, 32), seq(0.01, 0.15, length = 15), seq(0.15, 0.1, length = 15))))

stock.n(bet)[1,] <- 100000
stock.n(bet)[-1,1] <- 1000

qest <- rep(1e-5 * c(.1,.1,.9,1,1,.9,.5, 0.3, .3, .2), nyears)

# model data
data <- expand.grid(age = ages, year = years)

# model formulas
fmodel <- ~ s(age, k = 4) + s(year, k = 10)
qmodel <- list(~ s(age, k = 4))
srmodel <- ~ factor(year)
n1model <- ~ factor(age)

# design matrices
require(FLa4a)
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
Xf <- getX(fmodel, data)
Xq <- getX(qmodel[[1]], data)
Xr <- getX(srmodel, subset(data, age == ages[1]))
Xn1 <- getX(n1model, subset(data, age > ages[1] & year == years[1]))
options(opts)

# model parameters
bf <- c(solve(t(Xf) %*% Xf) %*% t(Xf) %*% log(c(harvest(bet))))
bq <- c(solve(t(Xq) %*% Xq) %*% t(Xq) %*% log(qest))
br <- c(solve(t(Xr) %*% Xr) %*% t(Xr) %*% log(c(rec(bet))))
bn1 <- c(solve(t(Xn1) %*% Xn1) %*% t(Xn1) %*% log(c(stock.n(bet)[-1,1])))
lsigc <- log(0.1)
lsigi <- log(0.1)

# m, wts and maturity
M <- matrix(c(m(bet)[,1]), nages, nyears)
Wt <- matrix(c(stock.wt(bet)[,1]), nages, nyears)
Mat <- matrix(c(mat(bet)[,1]), nages, nyears)
pz <- 0.6

